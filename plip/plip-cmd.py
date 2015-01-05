"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
Module for visualization in PyMOL.
Copyright (C) 2014  Sebastian Salentin

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
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

__version__ = '1.0.0'
descript = "Protein-Ligand Interaction Profiler (PLIP) v%s " \
           "is a command-line based tool to analyze interactions in a protein-ligand complex." % __version__


def sysexit(code, msg):
    """Exit using an custom error message and error code."""
    sys.stderr.write(msg)
    sys.exit(code)


def check_pdb_status(pdbid):
    """Returns the status and up-to-date entry in the PDB for a given PDB ID"""
    url = 'http://www.rcsb.org/pdb/rest/idStatus?structureId=%s' % pdbid
    xmlf = urllib2.urlopen(url)
    xml = et.parse(xmlf)
    xmlf.close()
    status = None
    current_pdbid = pdbid
    for df in xml.xpath('//record'):
        status = df.attrib['status']  # Status of an entry can be either 'UNKWOWN', 'OBSOLETE', or 'CURRENT'
        if status == 'OBSOLETE':
            current_pdbid = df.attrib['replacedBy']  # Contains the up-to-date PDB ID for obsolete entries
    return [status, current_pdbid.lower()]


def fetch_pdb(pdbid, verbose_mode):
    """Get the newest entry from the RCSB server for the given PDB ID. Exits with '1' if PDB ID is invalid."""
    pdbid = pdbid.lower()
    if verbose_mode:
        sys.stdout.write('Checking status of PDB ID %s ... ' % pdbid)
    state, current_entry = check_pdb_status(pdbid)  # Get state and current PDB ID

    if verbose_mode:
        if state == 'OBSOLETE':
            sys.stdout.write('entry is obsolete, getting %s instead.\n' % current_entry)
        elif state == 'CURRENT':
            sys.stdout.write('entry is up to date.\n')
    if state == 'UNKNOWN':
        sysexit(3, 'Invalid PDB ID')
    if verbose_mode:
        sys.stdout.write('Downloading file from PDB ... ')
    pdburl = 'http://www.rcsb.org/pdb/files/%s.pdb' % current_entry  # Get URL for current entry
    pdbfile = None
    try:
        pdbfile = urllib2.urlopen(pdburl).read()
    except urllib2.HTTPError:
        sysexit(5, "Error: No file in PDB format available from wwPDB for the given PDB ID.")
    return [pdbfile, current_entry]


def process_pdb(pdbfile, outpath, xml=False, verbose_mode=False, pics=False, pymol=False, maxthreads=None):
    """Analysis of a single PDB file. Can generate textual reports XML, PyMOL session files and images as output."""
    mol = PDBComplex()
    mol.output_path = outpath
    mol.load_pdb(pdbfile)

    # Begin constructing the XML tree
    report = et.Element('report')
    plipversion = et.SubElement(report, 'plipversion')
    plipversion.text = __version__
    pdbid = et.SubElement(report, 'pdbid')
    pdbid.text = mol.pymol_name.upper()

    # Write header of rST file
    textlines = ['Prediction of noncovalent interactions for PDB structure %s' % mol.pymol_name.upper(), ]
    textlines.append("="*len(textlines[0]))
    textlines.append('Created on %s using PLIP v%s\n' % (time.strftime("%Y/%m/%d"), __version__))

    if verbose_mode:
        num_ligs = len([site for site in mol.interaction_sets if not mol.interaction_sets[site].no_interactions])
        if num_ligs == 1:
            sys.stdout.write("Analyzing %s with one ligand.\n" % mol.pymol_name)
        elif num_ligs > 1:
            sys.stdout.write("Analyzing %s with %i ligands.\n" % (mol.pymol_name, num_ligs))
        else:
            sys.stdout.write("%s contains no ligands.\n" % mol.pymol_name)

    ######################################################################################################
    # Generate XML- and rST-formatted reports for each binding site and initialize visualization threads #
    ######################################################################################################

    threads = []
    running_threads = []

    for i, site in enumerate(mol.interaction_sets):
        s = mol.interaction_sets[site]
        bindingsite = TextOutput(s).generate_xml()
        bindingsite.set('id', str(i+1))
        bindingsite.set('has_interactions', 'False')
        report.insert(i+1, bindingsite)
        for itype in TextOutput(s).generate_rst():
            textlines.append(itype)
        if not s.no_interactions:
            bindingsite.set('has_interactions', 'True')
            if verbose_mode:
                sys.stdout.write("  @ %s\n" % site)

            if pymol:
                # Initialize thread for PyMOL visualization and add to list of threads to be processed
                p = multiprocessing.Process(target=visualize_in_pymol, args=(mol, site, False, pics))
                threads.append(p)

        else:
            textlines.append('No interactions detected.')
        sys.stdout = sys.__stdout__  # Change back to original stdout, gets changed when PyMOL has been used before

    ##############################################
    # Use multithreading for PyMOL visualization #
    ##############################################

    if maxthreads is None:  # Use as many threads as there are processor cores (should be a safe value)
        maxthreads = multiprocessing.cpu_count()
    else:
        maxthreads = max(2, maxthreads)
    while len(threads) != 0:
        # Add threads as long as there are free threads available
        for i, t in enumerate(threads):
            if len(running_threads) <= maxthreads-2:  # One can be still added, one used for the main process
                t.start()
                running_threads.append(t)
                threads.pop(i)
        # Check if running threads have finished and delete them
        for j, t in enumerate(running_threads):
            if not t.is_alive():
                running_threads.pop(j)

    ###########################################
    # Write final rST and XML to output files #
    ###########################################

    tree = et.ElementTree(report)
    create_folder_if_not_exists(tilde_expansion(outpath))
    if xml:
        tree.write('%s/report.xml' % tilde_expansion(outpath), pretty_print=True, xml_declaration=True)

    with open('%s/report.rst.txt' % tilde_expansion(outpath), 'w') as f:
        [f.write(textline+'\n') for textline in textlines]


def main(args):
    """Main function. Calls functions for processing, report generation and visualization."""
    pdbid, outp = None, None
    outp = "".join([args.outpath, '/']) if not args.outpath.endswith('/') else args.outpath

    if args.verbose:
        # Print title and version
        title = "* Protein-Ligand Interaction Profiler v%s *" % __version__
        sys.stdout.write('\n'+'*'*len(title)+'\n')
        sys.stdout.write(title)
        sys.stdout.write('\n'+'*'*len(title)+'\n\n')

    if args.input is not None:  # Process PDB file
        if os.path.getsize(args.input) == 0:
            sysexit(2, 'Error: Empty PDB file')  # Exit if input file is empty
        process_pdb(args.input, outp, xml=args.xml, verbose_mode=args.verbose, pics=args.pics, pymol=args.pymol,
                    maxthreads=int(args.maxthreads))
    else:  # Try to fetch the current PDB structure directly from the RCBS server
        try:
            pdbfile, pdbid = fetch_pdb(args.pdbid.lower(), verbose_mode=args.verbose)
            pdbpath = '%s/%s.pdb' % (args.outpath, pdbid)
            create_folder_if_not_exists(args.outpath)

            if args.verbose:
                sys.stdout.write('file downloaded as %s\n\n' % pdbpath)

            with open(tilde_expansion(pdbpath), 'w') as g:
                g.write(pdbfile)
            process_pdb(tilde_expansion(pdbpath), tilde_expansion(outp), xml=args.xml, verbose_mode=args.verbose,
                        pics=args.pics, pymol=args.pymol, maxthreads=int(args.maxthreads))
        except ValueError:  # Invalid PDB ID, cannot fetch from RCBS server
            sysexit(3, 'Error: Invalid PDB ID')
    if pdbid is not None and outp is not None:
        if outp in ['.', './']:
            outp = 'the working directory.'
        if args.verbose:
            sys.stdout.write('\nNearly finished with analysis of %s. Find the result files in %s\n\n' % (pdbid, outp))

if __name__ == '__main__':

    ##############################
    # Parse command line arguments
    ##############################

    parser = ArgumentParser(prog="PLIP", description=descript)
    pdbstructure = parser.add_mutually_exclusive_group(required=True)  # Needs either PDB ID or file
    pdbstructure.add_argument("-f", "--file", dest="input")
    pdbstructure.add_argument("-i", "--input", dest="pdbid")
    parser.add_argument("-o", "-out", dest="outpath", default="./")
    parser.add_argument("-v", "--verbose", dest="verbose", default=False, help="Set verbose mode", action="store_true")
    parser.add_argument("-p", "--pics", dest="pics", default=False, help="Additional pictures", action="store_true")
    parser.add_argument("-x", "--xml", dest="xml", default=False, help="Additional XML output for reports",
                        action="store_true")
    parser.add_argument("-y", "--pymol", dest="pymol", default=False, help="Additional PyMOL session files",
                        action="store_true")
    parser.add_argument("--maxthreads", dest="maxthreads", default=1,
                        help="Set maximum number of main threads (number of binding sites processed simultaneously)",
                        type=int)
    arguments = parser.parse_args()

    main(arguments)  # Start main script