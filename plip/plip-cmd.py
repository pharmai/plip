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

# External libraries
import lxml.etree as et

__version__ = '0.9.3'
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


def fetch_pdb(pdbid):
    """Get the newest entry from the RCSB server for the given PDB ID. Exits with '1' if PDB ID is invalid."""
    pdbid = pdbid.lower()
    status = check_pdb_status(pdbid)
    if status[0] == 'UNKNOWN':
        raise(ValueError, 'Invalid PDB ID')
    pdburl = 'http://www.rcsb.org/pdb/files/%s.pdb' % status[1]
    return [urllib2.urlopen(pdburl).read(), status[1]]


def process_pdb(pdbfile, outpath, is_zipped=False, text=False, verbose=False, pics=False):
    """Analysis of a single PDB file. Can generate textual reports and PyMOL session files as output."""
    tmpmol = PDBComplex()
    tmpmol.output_path = outpath
    tmpmol.load_pdb(pdbfile, is_zipped)
    report = et.Element('report')
    pdbid = et.SubElement(report, 'pdbid')
    pdbid.text = tmpmol.pymol_name.upper()
    textlines = []
    textlines.append('Prediction of noncovalent interactions for PDB structure %s' % tmpmol.pymol_name.upper())
    textlines.append("="*len(textlines[0]))
    textlines.append('Created on %s using PLIP v%s\n' % (time.strftime("%Y/%m/%d"), __version__))
    if verbose:
        num_ligs = len(tmpmol.interaction_sets)
        if num_ligs == 1:
            print("\nAnalyzing %s with one ligand." % tmpmol.pymol_name)
        elif num_ligs > 1:
            print("\nAnalyzing %s with %i ligands." % (tmpmol.pymol_name, num_ligs))
        else:
            print("\n%s contains no ligands." % tmpmol.pymol_name)

    ###########################################################################################
    # Generate XML-formatted reports for each binding site and send to stdout as one XML file #
    ###########################################################################################
    for i, site in enumerate(tmpmol.interaction_sets):
        s = tmpmol.interaction_sets[site]
        if verbose:
            print("  @ %s" % site)
        bindingsite = TextOutput(s).generate_xml()
        bindingsite.set('id', str(i+1))
        report.insert(i+1, bindingsite)
        visualize_in_pymol(tmpmol, site, show=False, pics=pics)
        sys.stdout = sys.__stdout__  # Change back to original stdout, gets changed when PyMOL has been used before
    tree = et.ElementTree(report)
    create_folder_if_not_exists(tilde_expansion(outpath))
    tree.write('%s/report.xml' % tilde_expansion(outpath), pretty_print=True, xml_declaration=True)
    #@todo Fully implement textual output as restructured text
    if text:
        with open('%s/report.rst' % tilde_expansion(outpath), 'w') as f:
            for textline in textlines:
                f.write(textline+'\n')


def main(args):
    """Main function. Processes command line arguments and processes PDB files in single or batch mode."""
    if args.verbose:
        print("*** Protein-Ligand Interaction Profiler v%s ***" % __version__)
    outp = "".join([args.outpath, '/']) if not args.outpath.endswith('/') else args.outpath
    if args.input is not None:  # Process PDB file
        if os.path.getsize(args.input) == 0:
            sysexit(2, 'Error: Empty PDB file')
        process_pdb(args.input, outp, text=args.txt, verbose=args.verbose, pics=args.pics)
    else:  # Try to fetch the current PDB structure directly from the RCBS server
        try:
            pdbfile, pdbid = fetch_pdb(args.pdbid.lower())
            pdbpath = '%s/%s.pdb' % (args.outpath, pdbid)
            create_folder_if_not_exists(args.outpath)
            with open(tilde_expansion(pdbpath), 'w') as g:
                g.write(pdbfile)
            process_pdb(tilde_expansion(pdbpath), tilde_expansion(outp), text=args.txt, verbose=args.verbose,
                        pics=args.pics)
        except ValueError:
            sysexit(3, 'Error: Invalid PDB ID')

if __name__ == '__main__':
    #@todo Fully implement flat text output for reports
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
    parser.add_argument("-t", "--text", dest="txt", default=False, help="Additional flat text output for reports",
                        action="store_true")
    arguments = parser.parse_args()
    main(arguments)