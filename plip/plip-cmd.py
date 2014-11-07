"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
Script for webserver usage.
"""

# Own modules
from modules.preparation import *
from modules.visualize import visualize_in_pymol
from modules.report import TextOutput

# Python standard library
import sys
from argparse import ArgumentParser
import urllib2

# External libraries
import lxml.etree as et


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


def process_pdb(pdbfile, outpath, is_zipped=False):
    """Analysis of a single PDB file. Can generate textual reports and PyMOL session files as output."""
    tmpmol = PDBComplex()
    tmpmol.output_path = outpath
    tmpmol.load_pdb(pdbfile, is_zipped)
    report = et.Element('report')
    pdbid = et.SubElement(report, 'pdbid')
    pdbid.text = tmpmol.pymol_name.upper()

    ###########################################################################################
    # Generate XML-formatted reports for each binding site and send to stdout as one XML file #
    ###########################################################################################
    for i, site in enumerate(tmpmol.interaction_sets):
        s = tmpmol.interaction_sets[site]
        bindingsite = TextOutput(s).generate_xml()
        bindingsite.set('id', str(i+1))
        report.insert(i+1, bindingsite)
        visualize_in_pymol(tmpmol, site, show=False, pics=True)
    sys.stdout = sys.__stdout__  # Change back to original stdout, gets changed when PyMOL has been used before
    tree = et.ElementTree(report)
    create_folder_if_not_exists(tilde_expansion(outpath))
    tree.write('%s/%s.xml' % (tilde_expansion(outpath), tmpmol.pymol_name), pretty_print=True, xml_declaration=True)


def main(args):
    """Main function. Processes command line arguments and processes PDB files in single or batch mode."""

    outp = "".join([args.outpath, '/']) if not args.outpath.endswith('/') else args.outpath
    if args.input is not None:  # Process PDB file
        if os.path.getsize(args.input) == 0:
            sysexit(2, 'Error: Empty PDB file')
        process_pdb(args.input, outp)
    else:  # Try to fetch the current PDB structure directly from the RCBS server
        try:
            pdbfile, pdbid = fetch_pdb(args.pdbid.lower())
            pdbpath = '%s/%s.pdb' % (args.outpath, pdbid)
            create_folder_if_not_exists(args.outpath)
            with open(tilde_expansion(pdbpath), 'w') as g:
                g.write(pdbfile)
            process_pdb(tilde_expansion(pdbpath), tilde_expansion(outp))
        except ValueError:
            sysexit(3, 'Error: Invalid PDB ID')

if __name__ == '__main__':

    ##############################
    # Parse command line arguments
    ##############################
    parser = ArgumentParser(prog="PLIP")
    pdbstructure = parser.add_mutually_exclusive_group(required=True)
    pdbstructure.add_argument("-f", dest="input")
    pdbstructure.add_argument("-i", dest="pdbid")
    parser.add_argument("-o", dest="outpath", default="/tmp/")
    arguments = parser.parse_args()
    main(arguments)