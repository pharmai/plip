"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
webservices.py - Connect to various webservices to retrieve data
"""

# Python standard library
from __future__ import absolute_import

try: # Python 3
    from urllib.request import urlopen
    from urllib.error import HTTPError
except ImportError: # Fallback Python 2.x
    from urllib2 import urlopen, HTTPError

# Own modules
from .supplemental import write_message, sysexit

# External libraries
import lxml.etree as et



def check_pdb_status(pdbid):
    """Returns the status and up-to-date entry in the PDB for a given PDB ID"""
    url = 'http://www.rcsb.org/pdb/rest/idStatus?structureId=%s' % pdbid
    xmlf = urlopen(url)
    xml = et.parse(xmlf)
    xmlf.close()
    status = None
    current_pdbid = pdbid
    for df in xml.xpath('//record'):
        status = df.attrib['status']  # Status of an entry can be either 'UNKWOWN', 'OBSOLETE', or 'CURRENT'
        if status == 'OBSOLETE':
            current_pdbid = df.attrib['replacedBy']  # Contains the up-to-date PDB ID for obsolete entries
    return [status, current_pdbid.lower()]


def fetch_pdb(pdbid):
    """Get the newest entry from the RCSB server for the given PDB ID. Exits with '1' if PDB ID is invalid."""
    pdbid = pdbid.lower()
    write_message('\nChecking status of PDB ID %s ... ' % pdbid)
    state, current_entry = check_pdb_status(pdbid)  # Get state and current PDB ID

    if state == 'OBSOLETE':
        write_message('entry is obsolete, getting %s instead.\n' % current_entry)
    elif state == 'CURRENT':
        write_message('entry is up to date.\n')
    elif state == 'UNKNOWN':
        sysexit(3, 'Invalid PDB ID (Entry does not exist on PDB server)\n')
    write_message('Downloading file from PDB ... ')
    pdburl = 'http://www.rcsb.org/pdb/files/%s.pdb' % current_entry  # Get URL for current entry
    try:
        pdbfile = urlopen(pdburl).read().decode()
        # If no PDB file is available, a text is now shown with "We're sorry, but ..."
        # Could previously be distinguished by an HTTP error
        if 'sorry' in pdbfile:
            sysexit(5, "No file in PDB format available from wwPDB for the given PDB ID.\n")
    except HTTPError:
        sysexit(5, "No file in PDB format available from wwPDB for the given PDB ID.\n")
    return [pdbfile, current_entry]
