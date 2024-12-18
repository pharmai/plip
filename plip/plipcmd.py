#! /usr/bin/env python
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
plipcmd.py - Main script for PLIP command line execution.
"""

# system imports
import argparse
import logging
import multiprocessing
import os
import sys
import ast
from argparse import ArgumentParser
from collections import namedtuple

from plip.basic import config, logger

logger = logger.get_logger()

from plip.basic.config import __version__
from plip.basic.parallel import parallel_fn
from plip.basic.remote import VisualizerData
from plip.exchange.report import StructureReport
from plip.exchange.webservices import fetch_pdb
from plip.structure.preparation import create_folder_if_not_exists, extract_pdbid
from plip.structure.preparation import tilde_expansion, PDBComplex

description = f"The Protein-Ligand Interaction Profiler (PLIP) Version {__version__} " \
              "is a command-line based tool to analyze interactions in a protein-ligand complex. " \
              "If you are using PLIP in your work, please cite: " \
              f"{config.__citation_information__} " \
              f"Supported and maintained by: {config.__maintainer__}"


def threshold_limiter(aparser, arg):
    arg = float(arg)
    if arg <= 0:
        aparser.error("All thresholds have to be values larger than zero.")
    return arg

def residue_list(input_string):
    """Parse mix of residue numbers and ranges passed with the --residues flag into one list"""
    result = []
    for part in input_string.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            result.extend(range(start, end + 1))
        else:
            result.append(int(part))
    return result

def process_pdb(pdbfile, outpath, as_string=False, outputprefix='report'):
    """Analysis of a single PDB file with optional chain filtering."""
    if not as_string:
        pdb_file_name = pdbfile.split('/')[-1]
        startmessage = f'starting analysis of {pdb_file_name}'
    else:
        startmessage = 'starting analysis from STDIN'
    logger.info(startmessage)
    mol = PDBComplex()
    mol.output_path = outpath
    mol.load_pdb(pdbfile, as_string=as_string)
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
        from plip.visualization.visualize import visualize_in_pymol
        complexes = [VisualizerData(mol, site) for site in sorted(mol.interaction_sets)
                     if not len(mol.interaction_sets[site].interacting_res) == 0]
        if config.MAXTHREADS > 1:
            logger.info(f'generating visualizations in parallel on {config.MAXTHREADS} cores')
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
            logger.error(f'invalid PDB-ID (wrong format): {inputpdbid}')
            sys.exit(1)
        pdbfile, pdbid = fetch_pdb(inputpdbid.lower())
        pdbpath = tilde_expansion('%s/%s.pdb' % (config.BASEPATH.rstrip('/'), pdbid))
        create_folder_if_not_exists(config.BASEPATH)
        with open(pdbpath, 'w') as g:
            g.write(pdbfile)
        logger.info(f'file downloaded as {pdbpath}')
        return pdbpath, pdbid

    except ValueError:  # Invalid PDB ID, cannot fetch from RCBS server
        logger.error(f'PDB-ID does not exist: {inputpdbid}')
        sys.exit(1)


def remove_duplicates(slist):
    """Checks input lists for duplicates and returns
    a list with unique entries"""
    unique = list(set(slist))
    difference = len(slist) - len(unique)
    if difference == 1:
        logger.info('removed one duplicate entry from input list')
    if difference > 1:
        logger.info(f'Removed {difference} duplicate entries from input list')
    return unique


def run_analysis(inputstructs, inputpdbids, chains=None):
    """Main function. Calls functions for processing, report generation and visualization."""
    pdbid, pdbpath = None, None
    # @todo For multiprocessing, implement better stacktracing for errors
    # Print title and version
    logger.info(f'Protein-Ligand Interaction Profiler (PLIP) {__version__}')
    logger.info(f'brought to you by: {config.__maintainer__}')
    logger.info(f'please cite: {config.__citation_information__}')
    output_prefix = config.OUTPUTFILENAME

    if inputstructs is not None:  # Process PDB file(s)
        num_structures = len(inputstructs)  # @question: how can it become more than one file? The tilde_expansion function does not consider this case.
        inputstructs = remove_duplicates(inputstructs)
        read_from_stdin = False
        for inputstruct in inputstructs:
            if inputstruct == '-':  # @expl: when user gives '-' as input, pdb file is read from stdin
                inputstruct = sys.stdin.read()
                read_from_stdin = True
                if config.RAWSTRING:
                    if sys.version_info < (3,):
                        inputstruct = bytes(inputstruct).decode('unicode_escape')  # @expl: in Python2, the bytes object is just a string.
                    else:
                        inputstruct = bytes(inputstruct, 'utf8').decode('unicode_escape')
            else:
                if os.path.getsize(inputstruct) == 0:
                    logger.error('empty PDB file')
                    sys.exit(1)
                if num_structures > 1:
                    basename = inputstruct.split('.')[-2].split('/')[-1]
                    config.OUTPATH = '/'.join([config.BASEPATH, basename])
                    output_prefix = 'report'
            process_pdb(inputstruct, config.OUTPATH, as_string=read_from_stdin, outputprefix=output_prefix)
    else:  # Try to fetch the current PDB structure(s) directly from the RCBS server
        num_pdbids = len(inputpdbids)
        inputpdbids = remove_duplicates(inputpdbids)
        for inputpdbid in inputpdbids:
            pdbpath, pdbid = download_structure(inputpdbid)
            if num_pdbids > 1:
                config.OUTPATH = '/'.join([config.BASEPATH, pdbid[1:3].upper(), pdbid.upper()])
                output_prefix = 'report'
            process_pdb(pdbpath, config.OUTPATH, outputprefix=output_prefix)

    if (pdbid is not None or inputstructs is not None) and config.BASEPATH is not None:
        if config.BASEPATH in ['.', './']:
            logger.info('finished analysis, find the result files in the working directory')
        else:
            logger.info(f'finished analysis, find the result files in {config.BASEPATH}')


def main():
    """Parse command line arguments and start main script for analysis."""
    parser = ArgumentParser(prog="PLIP", description=description)
    pdbstructure = parser.add_mutually_exclusive_group(required=True)  # Needs either PDB ID or file
    # '-' as file name reads from stdin
    pdbstructure.add_argument("-f", "--file", dest="input", nargs="+", help="Set input file, '-' reads from stdin")
    pdbstructure.add_argument("-i", "--input", dest="pdbid", nargs="+")
    outputgroup = parser.add_mutually_exclusive_group(required=False)  # Needs either outpath or stdout
    outputgroup.add_argument("-o", "--out", dest="outpath", default="./")
    outputgroup.add_argument("-O", "--stdout", dest="stdout", action="store_true", default=False,
                             help="Write to stdout instead of file")
    parser.add_argument("--rawstring", dest="use_raw_string", default=False, action="store_true",
                        help="Use Python raw strings for stdin")
    parser.add_argument("-v", "--verbose", dest="verbose", default=False, help="Turn on verbose mode",
                        action="store_true")
    parser.add_argument("-q", "--quiet", dest="quiet", default=False, help="Turn on quiet mode", action="store_true")
    parser.add_argument("-s", "--silent", dest="silent", default=False, help="Turn on silent mode", action="store_true")
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
                        help="Treat nucleic acids as part of the receptor structure (together with any present protein) instead of as a ligand.",
                        action="store_true")
    parser.add_argument("--name", dest="outputfilename", default="report",
                        help="Set a filename for the report TXT and XML files. Will only work when processing single structures.")
    ligandtype = parser.add_mutually_exclusive_group()  # Either peptide/inter or intra mode
    ligandtype.add_argument("--peptides", "--inter", dest="peptides", default=[],
                            help="Allows to define one or multiple chains as peptide ligands or to detect inter-chain contacts",
                            nargs="+")
    ligandtype.add_argument("--intra", dest="intra", help="Allows to define one chain to analyze intra-chain contacts.")
    parser.add_argument("--residues", dest="residues", default=[], nargs="+",
                        help="""Allows to specify which residues of the chain(s) should be considered as peptide ligands.
                        Give single residues (separated with comma) or ranges (with dash) or both, for several chains separate selections with one space""")
    parser.add_argument("--keepmod", dest="keepmod", default=False,
                        help="Keep modified residues as ligands",
                        action="store_true")
    parser.add_argument("--nohydro", dest="nohydro", default=False,
                        help="Do not add polar hydrogens in case your structure already contains hydrogens.",
                        action="store_true")
    parser.add_argument("--model", dest="model", default=1, type=int,
                        help="Model number to be used for multi-model structures.")
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

    # Add argument to define receptor and ligand chains
    parser.add_argument("--chains", dest="chains", type=str,
                        help="Specify chains as receptor/ligand groups, e.g., '[['A'], ['B']]'. "
                             "Use format [['A'], ['B', 'C']] to define A as receptor, and B, C as ligands.")


    arguments = parser.parse_args()
    # make sure, residues is only used together with --inter (could be expanded to --intra in the future)
    if arguments.residues and not (arguments.peptides or arguments.intra):
        parser.error("The --residues option requires specification of a chain with --inter or --peptide")
    if arguments.residues and len(arguments.residues)!=len(arguments.peptides):
        parser.error("Please provide residue numbers or ranges for each chain specified. Separate selections with a single space.")
    # configure log levels
    config.VERBOSE = True if arguments.verbose else False
    config.QUIET = True if arguments.quiet else False
    config.SILENT = True if arguments.silent else False
    if config.VERBOSE:
        logger.setLevel(logging.DEBUG)
    elif config.QUIET:
        logger.setLevel(logging.WARN)
    elif config.SILENT:
        logger.setLevel(logging.CRITICAL)
    else:
        logger.setLevel(config.DEFAULT_LOG_LEVEL)
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
    config.RESIDUES = dict(zip(arguments.peptides, map(residue_list, arguments.residues)))
    config.INTRA = arguments.intra
    config.NOFIX = arguments.nofix
    config.NOFIXFILE = arguments.nofixfile
    config.NOPDBCANMAP = bool(arguments.nopdbcanmap or config.INTRA or config.PEPTIDES)
    config.KEEPMOD = arguments.keepmod
    config.DNARECEPTOR = arguments.dnareceptor
    config.OUTPUTFILENAME = arguments.outputfilename
    config.NOHYDRO = arguments.nohydro
    config.MODEL = arguments.model

    try:
        # add inner quotes for python backend
        if not arguments.chains:
            config.CHAINS = None
        else:
            import re
            quoted_input = re.sub(r'(?<!["\'])\b([a-zA-Z0-9_]+)\b(?!["\'])', r'"\1"', arguments.chains)
            config.CHAINS = ast.literal_eval(quoted_input)
        print(config.CHAINS)
        if config.CHAINS and not all(isinstance(c, list) for c in config.CHAINS):
            raise ValueError("Chains should be specified as a list of lists, e.g., '[[A], [B, C]]'.")
    except (ValueError, SyntaxError):
        parser.error("The --chains option must be in the format '[[A], [B, C]]'.")


    # Make sure we have pymol with --pics and --pymol
    if config.PICS or config.PYMOL:
        try:
            import pymol
        except ImportError:
            logger.error('PyMOL is required for the --pics and --pymol option')
            sys.exit(1)
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
    run_analysis(expanded_path, arguments.pdbid)  # Start main script


if __name__ == '__main__':
    main()
