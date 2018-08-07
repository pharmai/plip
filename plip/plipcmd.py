#! /usr/bin/env python
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
plipcmd.py - Main script for PLIP command line execution.
"""

# Compatibility
from __future__ import print_function
from __future__ import absolute_import

# Own modules
try:
    from plip.modules.preparation import *
    from plip.modules.plipremote import VisualizerData
    from plip.modules.report import StructureReport,__version__
    from plip.modules import config
    from plip.modules.mp import parallel_fn
    from plip.modules.webservices import check_pdb_status, fetch_pdb
except ImportError:
    from modules.preparation import *
    from modules.plipremote import VisualizerData
    from modules.report import StructureReport,__version__
    from modules import config
    from modules.mp import parallel_fn
    from modules.webservices import check_pdb_status, fetch_pdb

# Python standard library
import sys
import argparse
from argparse import ArgumentParser
import time
import multiprocessing
import json

# External libraries
import lxml.etree as et

descript = "Protein-Ligand Interaction Profiler (PLIP) v%s " \
           "is a command-line based tool to analyze interactions in a protein-ligand complex. " \
           "If you are using PLIP in your work, please cite: " \
           "Salentin,S. et al. PLIP: fully automated protein-ligand interaction profiler. " \
           "Nucl. Acids Res. (1 July 2015) 43 (W1): W443-W447. doi: 10.1093/nar/gkv315" % __version__


def threshold_limiter(aparser, arg):
    arg = float(arg)
    if arg <= 0:
        aparser.error("All thresholds have to be values larger than zero.")
    return arg





def process_pdb(pdbfile, outpath, as_string=False, outputprefix='report'):
    """Analysis of a single PDB file. Can generate textual reports XML, PyMOL session files and images as output."""
    if not as_string:
        startmessage = '\nStarting analysis of %s\n' % pdbfile.split('/')[-1]
    else:
        startmessage = "Starting analysis from stdin.\n"
    write_message(startmessage)
    write_message('='*len(startmessage)+'\n')
    mol = PDBComplex()
    mol.output_path = outpath
    mol.load_pdb(pdbfile, as_string=as_string)
    # #@todo Offers possibility for filter function from command line (by ligand chain, position, hetid)
    for ligand in mol.ligands:
        mol.characterize_complex(ligand)

    create_folder_if_not_exists(outpath)

    # Generate the report files
    streport = StructureReport(mol, outputprefix=outputprefix)

    config.MAXTHREADS = min(config.MAXTHREADS, len(mol.interaction_sets))


    ######################################
    # PyMOL Visualization (parallelized) #
    ######################################

    if config.PYMOL or config.PICS:
        try:
            from plip.modules.visualize import visualize_in_pymol
        except ImportError:
            from modules.visualize import visualize_in_pymol
        complexes = [VisualizerData(mol, site) for site in sorted(mol.interaction_sets)
                     if not len(mol.interaction_sets[site].interacting_res) == 0]
        if config.MAXTHREADS > 1:
            write_message('\nGenerating visualizations in parallel on %i cores ...' % config.MAXTHREADS)
            parfn = parallel_fn(visualize_in_pymol)
            parfn(complexes, processes=config.MAXTHREADS)
        else:
            [visualize_in_pymol(plcomplex) for plcomplex in complexes]

    if config.XML:  # Generate report in xml format
        streport.write_xml(as_string=config.STDOUT)

    if config.TXT:  # Generate report in txt (rst) format
        streport.write_txt(as_string=config.STDOUT)

def download_structure(inputpdbid):
    """Given a PDB ID, downloads the corresponding PDB structure.
    Checks for validity of ID and handles error while downloading.
    Returns the path of the downloaded file."""
    try:
        if len(inputpdbid) != 4 or extract_pdbid(inputpdbid.lower()) == 'UnknownProtein':
            sysexit(3, 'Invalid PDB ID (Wrong format)\n')
        pdbfile, pdbid = fetch_pdb(inputpdbid.lower())
        pdbpath = tilde_expansion('%s/%s.pdb' % (config.BASEPATH.rstrip('/'), pdbid))
        create_folder_if_not_exists(config.BASEPATH)
        with open(pdbpath, 'w') as g:
            g.write(pdbfile)
        write_message('file downloaded as %s\n\n' % pdbpath)
        return pdbpath, pdbid

    except ValueError:  # Invalid PDB ID, cannot fetch from RCBS server
        sysexit(3, 'Invalid PDB ID (Entry does not exist)\n')

def remove_duplicates(slist):
    """Checks input lists for duplicates and returns
    a list with unique entries"""
    unique = list(set(slist))
    difference = len(slist) - len(unique)
    if difference == 1:
        write_message("Removed one duplicate entry from input list.\n")
    if difference > 1:
        write_message("Removed %i duplicate entries from input list.\n" % difference)
    return unique


def main(inputstructs, inputpdbids):
    """Main function. Calls functions for processing, report generation and visualization."""
    pdbid, pdbpath = None, None
    # #@todo For multiprocessing, implement better stacktracing for errors
    # Print title and version
    title = "* Protein-Ligand Interaction Profiler v%s *" % __version__
    write_message('\n' + '*' * len(title) + '\n')
    write_message(title)
    write_message('\n' + '*' * len(title) + '\n\n')
    outputprefix = config.OUTPUTFILENAME

    if inputstructs is not None:  # Process PDB file(s)
        num_structures = len(inputstructs)
        inputstructs = remove_duplicates(inputstructs)
        read_from_stdin = False
        for inputstruct in inputstructs:
            if inputstruct == '-':
                inputstruct = sys.stdin.read()
                read_from_stdin = True
                if config.RAWSTRING:
                    if sys.version_info < (3,):
                        inputstruct = bytes(inputstruct).decode('unicode_escape')
                    else:
                        inputstruct = bytes(inputstruct, 'utf8').decode('unicode_escape')
            else:
                if os.path.getsize(inputstruct) == 0:
                    sysexit(2, 'Empty PDB file\n')  # Exit if input file is empty
                if num_structures > 1:
                    basename = inputstruct.split('.')[-2].split('/')[-1]
                    config.OUTPATH = '/'.join([config.BASEPATH, basename])
                    outputprefix = 'report'
            process_pdb(inputstruct, config.OUTPATH, as_string=read_from_stdin, outputprefix=outputprefix)
    else:  # Try to fetch the current PDB structure(s) directly from the RCBS server
        num_pdbids = len(inputpdbids)
        inputpdbids =remove_duplicates(inputpdbids)
        for inputpdbid in inputpdbids:
            pdbpath, pdbid = download_structure(inputpdbid)
            if num_pdbids > 1:
                config.OUTPATH = '/'.join([config.BASEPATH, pdbid[1:3].upper(), pdbid.upper()])
                outputprefix = 'report'
            process_pdb(pdbpath, config.OUTPATH, outputprefix=outputprefix)

    if (pdbid is not None or inputstructs is not None) and config.BASEPATH is not None:
        if config.BASEPATH in ['.', './']:
            write_message('\nFinished analysis. Find the result files in the working directory.\n\n')
        else:
            write_message('\nFinished analysis. Find the result files in %s\n\n' % config.BASEPATH)


    ##############################
    # Parse command line arguments
    ##############################
def main_init():
    """Parse command line arguments and start main script for analysis."""
    parser = ArgumentParser(prog="PLIP", description=descript)
    pdbstructure = parser.add_mutually_exclusive_group(required=True)  # Needs either PDB ID or file
    pdbstructure.add_argument("-f", "--file", dest="input", nargs="+", help="Set input file, '-' reads from stdin") # '-' as file name reads from stdin
    pdbstructure.add_argument("-i", "--input", dest="pdbid", nargs="+")
    outputgroup = parser.add_mutually_exclusive_group(required=False)  # Needs either outpath or stdout
    outputgroup.add_argument("-o", "--out", dest="outpath", default="./")
    outputgroup.add_argument("-O", "--stdout", dest="stdout", action="store_true", default=False, help="Write to stdout instead of file")
    parser.add_argument("--rawstring", dest="use_raw_string", default=False, action="store_true", help="Use Python raw strings for stdout and stdin")
    parser.add_argument("-v", "--verbose", dest="verbose", default=False, help="Set verbose mode", action="store_true")
    parser.add_argument("-p", "--pics", dest="pics", default=False, help="Additional pictures", action="store_true")
    parser.add_argument("-x", "--xml", dest="xml", default=False, help="Generate report file in XML format",
                        action="store_true")
    parser.add_argument("-t", "--txt", dest="txt", default=False, help="Generate report file in TXT (RST) format",
                        action="store_true")
    parser.add_argument("-y", "--pymol", dest="pymol", default=False, help="Additional PyMOL session files",
                        action="store_true")
    parser.add_argument("--maxthreads", dest="maxthreads", default=multiprocessing.cpu_count(),
                        help="Set maximum number of main threads (number of binding sites processed simultaneously)."
                             "If not set, PLIP uses all available CPUs if possible.",
                        type=int)
    parser.add_argument("--breakcomposite", dest="breakcomposite", default=False,
                        help="Don't combine ligand fragments with covalent bonds but treat them as single ligands for the analysis.",
                        action="store_true")
    parser.add_argument("--altlocation", dest="altlocation", default=False,
                        help="Also consider alternate locations for atoms (e.g. alternate conformations).",
                        action="store_true")
    parser.add_argument("--debug", dest="debug", default=False,
                        help="Turn on DEBUG mode with extended log.",
                        action="store_true")
    parser.add_argument("--nofix", dest="nofix", default=False,
                        help="Turns off fixing of PDB files.",
                        action="store_true")
    parser.add_argument("--nofixfile", dest="nofixfile", default=False,
                        help="Turns off writing files for fixed PDB files.",
                        action="store_true")
    parser.add_argument("--nopdbcanmap", dest="nopdbcanmap", default=False,
                        help="Turns off calculation of mapping between canonical and PDB atom order for ligands.",
                        action="store_true")
    parser.add_argument("--dnareceptor", dest="dnareceptor", default=False,
                        help="Uses the DNA instead of the protein as a receptor for interactions.",
                        action="store_true")
    parser.add_argument("--name", dest="outputfilename", default="report",
                        help="Set a filename for the report TXT and XML files. Will only work when processing single structures.")
    ligandtype = parser.add_mutually_exclusive_group()  # Either peptide/inter or intra mode
    ligandtype.add_argument("--peptides", "--inter", dest="peptides", default=[],
                        help="Allows to define one or multiple chains as peptide ligands or to detect inter-chain contacts",
                        nargs="+")
    ligandtype.add_argument("--intra", dest="intra", help="Allows to define one chain to analyze intra-chain contacts.")
    parser.add_argument("--keepmod", dest="keepmod", default=False,
                        help="Keep modified residues as ligands",
                        action="store_true")
    # Optional threshold arguments, not shown in help
    thr = namedtuple('threshold', 'name type')
    thresholds = [thr(name='aromatic_planarity', type='angle'),
                  thr(name='hydroph_dist_max', type='distance'), thr(name='hbond_dist_max', type='distance'),
                  thr(name='hbond_don_angle_min', type='angle'), thr(name='pistack_dist_max', type='distance'),
                  thr(name='pistack_ang_dev', type='other'), thr(name='pistack_offset_max', type='distance'),
                  thr(name='pication_dist_max', type='distance'), thr(name='saltbridge_dist_max', type='distance'),
                  thr(name='halogen_dist_max', type='distance'), thr(name='halogen_acc_angle', type='angle'),
                  thr(name='halogen_don_angle', type='angle'), thr(name='halogen_angle_dev', type='other'),
                  thr(name='water_bridge_mindist', type='distance'), thr(name='water_bridge_maxdist', type='distance'),
                  thr(name='water_bridge_omega_min', type='angle'), thr(name='water_bridge_omega_max', type='angle'),
                  thr(name='water_bridge_theta_min', type='angle')]
    for t in thresholds:
        parser.add_argument('--%s' % t.name, dest=t.name, type=lambda val: threshold_limiter(parser, val),
                            help=argparse.SUPPRESS)

    arguments = parser.parse_args()
    config.VERBOSE = True if (arguments.verbose or arguments.debug) else False
    config.DEBUG = True if arguments.debug else False
    config.MAXTHREADS = arguments.maxthreads
    config.XML = arguments.xml
    config.TXT = arguments.txt
    config.PICS = arguments.pics
    config.PYMOL = arguments.pymol
    config.STDOUT = arguments.stdout
    config.RAWSTRING = arguments.use_raw_string
    config.OUTPATH = arguments.outpath
    config.OUTPATH = tilde_expansion("".join([config.OUTPATH, '/'])
                                     if not config.OUTPATH.endswith('/') else config.OUTPATH)
    config.BASEPATH = config.OUTPATH  # Used for batch processing
    config.BREAKCOMPOSITE = arguments.breakcomposite
    config.ALTLOC = arguments.altlocation
    config.PEPTIDES = arguments.peptides
    config.INTRA = arguments.intra
    config.NOFIX = arguments.nofix
    config.NOFIXFILE = arguments.nofixfile
    config.NOPDBCANMAP = arguments.nopdbcanmap
    config.KEEPMOD = arguments.keepmod
    config.DNARECEPTOR = arguments.dnareceptor
    config.OUTPUTFILENAME = arguments.outputfilename
    # Make sure we have pymol with --pics and --pymol
    if config.PICS or config.PYMOL:
        try:
            import pymol
        except ImportError:
            write_message("PyMOL is required for --pics and --pymol.\n", mtype='error')
            raise
    # Assign values to global thresholds
    for t in thresholds:
        tvalue = getattr(arguments, t.name)
        if tvalue is not None:
            if t.type == 'angle' and not 0 < tvalue < 180:  # Check value for angle thresholds
                parser.error("Threshold for angles need to have values within 0 and 180.")
            if t.type == 'distance':
                if tvalue > 10:  # Check value for angle thresholds
                    parser.error("Threshold for distances must not be larger than 10 Angstrom.")
                elif tvalue > config.BS_DIST + 1:  # Dynamically adapt the search space for binding site residues
                    config.BS_DIST = tvalue + 1
            setattr(config, t.name.upper(), tvalue)
    # Check additional conditions for interdependent thresholds
    if not config.HALOGEN_ACC_ANGLE > config.HALOGEN_ANGLE_DEV:
        parser.error("The halogen acceptor angle has to be larger than the halogen angle deviation.")
    if not config.HALOGEN_DON_ANGLE > config.HALOGEN_ANGLE_DEV:
        parser.error("The halogen donor angle has to be larger than the halogen angle deviation.")
    if not config.WATER_BRIDGE_MINDIST < config.WATER_BRIDGE_MAXDIST:
        parser.error("The water bridge minimum distance has to be smaller than the water bridge maximum distance.")
    if not config.WATER_BRIDGE_OMEGA_MIN < config.WATER_BRIDGE_OMEGA_MAX:
        parser.error("The water bridge omega minimum angle has to be smaller than the water bridge omega maximum angle")
    expanded_path = tilde_expansion(arguments.input) if arguments.input is not None else None
    main(expanded_path, arguments.pdbid)  # Start main script
if __name__ == '__main__':
  main_init()
