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


class PDBRepo(object):
    """Copy from Joachim's repo"""

    def __init__(self, pdbpath='~/lib/pdb', obsoletefile='~/lib/pdbsuppl/obsolete.dat'):
        # read mapping of obsolete PDB ids
        path = os.path.expanduser(obsoletefile)
        f = open(path, 'r')

        self.succ = dict()  # dictionary mapping a pdb id to its successor

        for line in f:
            line = line.strip().split()
            try:
                self.succ[line[2].lower()] = line[3].lower()
            except IndexError:
                pass
        f.close()
        self.pdbpath = pdbpath

    def getpdbid(self, pdbid):
        return pdbid.lower()

    def getlatestpdbid(self, pdbid):
        """Checks the list of obsolete PDB IDs and returns the successor."""
        pdbid = self.getpdbid(pdbid)
        if pdbid in self.succ:
            pdbid = self.succ[pdbid]
        return pdbid

    def getpdbfile(self, pdbid):
        """Returns the path to the file of a given pdb id."""
        pdbid = self.getlatestpdbid(pdbid)
        pdbdir = os.path.expanduser(self.pdbpath)
        pdbfile = os.path.join(pdbdir, pdbid[1:3], 'pdb'+pdbid+'.ent.gz')
        if not os.path.exists(pdbfile):
            raise IOError("No such file: '"+pdbfile+"'")
        return pdbfile


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
    tree.write('%s/%s.xml' % (outpath, tmpmol.pymol_name), pretty_print=True, xml_declaration=True)


def main(args):
    """Main function. Processes command line arguments and processes PDB files in single or batch mode."""

    outp = "".join([args.outpath, '/']) if not args.outpath.endswith('/') else args.outpath
    if args.input is not None:  # Process PDB file
        #@todo Implement checking function for uploaded PDB files here. e.g. os.path.getsize(f) == 0
        process_pdb(args.input, outp)
    else:  # Fetch from biodata and process PDB file if corresponding PDB can be found
        try:
            pdbfile = pdbrepo.getpdbfile(args.pdbid)
            process_pdb(pdbfile, outp, is_zipped=True)  # All files in the repo are ent.gz archives
        except IOError:
            sys.stderr.write("Error: PDB file for PDB ID %s not found in the database in %s\n"
                             % (args.pdbid, pdbrepo.pdbpath))
            sys.exit()

if __name__ == '__main__':

    ##############################
    # Parse command line arguments
    ##############################
    parser = ArgumentParser(prog="PLIP Webserver")
    pdbstructure = parser.add_mutually_exclusive_group(required=True)
    pdbstructure.add_argument("-f", dest="input")
    pdbstructure.add_argument("-i", dest="pdbid")
    parser.add_argument("-o", dest="outpath", default="/tmp/")
    parser.add_argument("-pdbrepo", dest="pdbrepo", default="/biodata/biodb/ABG/databases/pdb")
    parser.add_argument("-obsolete", dest="obsolete", default="~/lib/pdbsuppl/obsolete.dat")
    arguments = parser.parse_args()

    pdbrepo = PDBRepo(pdbpath=arguments.pdbrepo, obsoletefile=arguments.obsolete)
    main(arguments)