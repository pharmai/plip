"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
"""

# Compatibility
from __future__ import print_function

# Own modules
from modules.preparation import *
from modules.visualize import visualize_in_pymol
from modules.report import TextOutput

# Python standard library
import sys
from argparse import ArgumentParser
import urllib2
import time
import multiprocessing

# External libraries
import lxml.etree as et

__version__ = '0.9.5'
descript = "Protein-Ligand Interaction Profiler (PLIP) v%s " \
           "is a command-line based tool to analyze interactions in a protein-ligand complex." % __version__


def sysexit(code, msg):
    sys.stderr.write(msg)
    sys.exit(code)


def check_pdb_status(pdbid):
    """Returns the status of a file in the PDB"""
    url = 'http://www.rcsb.org/pdb/rest/idStatus?structureId=%s' % pdbid
    xmlf = urllib2.urlopen(url)
    xml = et.parse(xmlf)
    xmlf.close()
    status = None
    current_pdbid = pdbid
    for df in xml.xpath('//record'):
        status = df.attrib['status']
        if status == 'OBSOLETE':
            current_pdbid = df.attrib['replacedBy']
    return [status, current_pdbid.lower()]


def fetch_pdb(pdbid, verbose_mode):
    """Get the newest entry from the RCSB server for the given PDB ID. Exits with '1' if PDB ID is invalid."""
    pdbid = pdbid.lower()
    if verbose_mode:
        sys.stdout.write('Checking status of PDB ID %s ... ' % pdbid)
    status = check_pdb_status(pdbid)
    if verbose_mode:
        if status[0] == 'OBSOLETE':
            sys.stdout.write('entry is obsolete, getting %s instead.\n' % status[1])
        elif status[0] == 'UNKOWN':
            pass
        else:
            sys.stdout.write('entry is up to date.\n')
    if status[0] == 'UNKNOWN':
        raise(ValueError, 'Invalid PDB ID')
    if verbose_mode:
        sys.stdout.write('Downloading file from PDB ... ')
    pdburl = 'http://www.rcsb.org/pdb/files/%s.pdb' % status[1]
    return [urllib2.urlopen(pdburl).read(), status[1]]


def process_pdb(pdbfile, outpath, is_zipped=False, text=False, verbose=False, pics=False):
    """Analysis of a single PDB file. Can generate textual reports and PyMOL session files as output."""
    tmpmol = PDBComplex()
    tmpmol.output_path = outpath
    tmpmol.load_pdb(pdbfile, is_zipped)
    report = et.Element('report')
    plipversion = et.SubElement(report, 'plipversion')
    plipversion.text = __version__
    pdbid = et.SubElement(report, 'pdbid')
    pdbid.text = tmpmol.pymol_name.upper()
    textlines = []
    textlines.append('Prediction of noncovalent interactions for PDB structure %s' % tmpmol.pymol_name.upper())
    textlines.append("="*len(textlines[0]))
    textlines.append('Created on %s using PLIP v%s\n' % (time.strftime("%Y/%m/%d"), __version__))
    if verbose:
        num_ligs = len([site for site in tmpmol.interaction_sets if not tmpmol.interaction_sets[site].no_interactions])
        if num_ligs == 1:
            sys.stdout.write("Analyzing %s with one ligand.\n" % tmpmol.pymol_name)
        elif num_ligs > 1:
            sys.stdout.write("Analyzing %s with %i ligands.\n" % (tmpmol.pymol_name, num_ligs))
        else:
            sys.stdout.write("%s contains no ligands.\n" % tmpmol.pymol_name)

    ###########################################################################################
    # Generate XML- and rST-formatted reports for each binding site#
    ###########################################################################################
    threads = []
    #@todo Improve thread handling
    for i, site in enumerate(tmpmol.interaction_sets):
        s = tmpmol.interaction_sets[site]
        bindingsite = TextOutput(s).generate_xml()
        bindingsite.set('id', str(i+1))
        if not s.no_interactions:
            bindingsite.set('has_interactions', 'True')
        else:
            bindingsite.set('has_interactions', 'False')
        report.insert(i+1, bindingsite)
        for itype in TextOutput(s).generate_rst():
            textlines.append(itype)
        if not s.no_interactions:
            if verbose:
                sys.stdout.write("  @ %s\n" % site)
            p = multiprocessing.Process(target=visualize_in_pymol, args=(tmpmol, site, False, pics))
            threads.append(p)
        else:
            textlines.append('No interactions detected.')
        sys.stdout = sys.__stdout__  # Change back to original stdout, gets changed when PyMOL has been used before
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    tree = et.ElementTree(report)
    create_folder_if_not_exists(tilde_expansion(outpath))
    tree.write('%s/report.xml' % tilde_expansion(outpath), pretty_print=True, xml_declaration=True)

    ################################
    # Write rST to the output file #
    ################################

    if text:
        with open('%s/report.rst' % tilde_expansion(outpath), 'w') as f:
            for textline in textlines:
                f.write(textline+'\n')


def main(args):
    """Main function. Processes command line arguments and processes PDB files in single or batch mode."""
    pdbid, outp = None, None
    if args.verbose:
        title = "* Protein-Ligand Interaction Profiler v%s *" % __version__
        sys.stdout.write('\n'+'*'*len(title)+'\n')
        sys.stdout.write(title)
        sys.stdout.write('\n'+'*'*len(title)+'\n\n')

    outp = "".join([args.outpath, '/']) if not args.outpath.endswith('/') else args.outpath
    if args.input is not None:  # Process PDB file
        if os.path.getsize(args.input) == 0:
            sysexit(2, 'Error: Empty PDB file')
        process_pdb(args.input, outp, text=args.txt, verbose=args.verbose, pics=args.pics)
    else:  # Try to fetch the current PDB structure directly from the RCBS server
        try:
            pdbfile, pdbid = fetch_pdb(args.pdbid.lower(), verbose_mode=args.verbose)
            pdbpath = '%s/%s.pdb' % (args.outpath, pdbid)
            if args.verbose:
                sys.stdout.write('file downloaded as %s\n\n' % pdbpath)
            create_folder_if_not_exists(args.outpath)
            with open(tilde_expansion(pdbpath), 'w') as g:
                g.write(pdbfile)
            process_pdb(tilde_expansion(pdbpath), tilde_expansion(outp), text=args.txt, verbose=args.verbose,
                        pics=args.pics)
        except ValueError:
            sysexit(3, 'Error: Invalid PDB ID')
    if pdbid is not None and outp is not None:
        if outp in ['.', './']:
            outp = 'the working directory.'
        sys.stdout.write('\nFinished with analysis of %s. Find the result files in %s\n\n' % (pdbid, outp))

if __name__ == '__main__':
    ##############################
    # Parse command line arguments
    ##############################
    parser = ArgumentParser(prog="PLIP", description=descript)
    pdbstructure = parser.add_mutually_exclusive_group(required=True)
    pdbstructure.add_argument("-f", "--file", dest="input")
    pdbstructure.add_argument("-i", "--input", dest="pdbid")
    parser.add_argument("-o", "-out", dest="outpath", default="./")
    parser.add_argument("-v", "--verbose", dest="verbose", default=False, help="Set verbose mode", action="store_true")
    parser.add_argument("-p", "--pics", dest="pics", default=False, help="Additional pictures", action="store_true")
    parser.add_argument("-t", "--text", dest="txt", default=False, help="Additional rST output for reports",
                        action="store_true")
    arguments = parser.parse_args()
    main(arguments)