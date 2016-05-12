"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
webservices.py - Connect to various webservices to retrieve data
Copyright 2014-2016 Sebastian Salentin

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

# Own modules
from supplemental import write_message, sysexit

# Python standard libary
import urllib2

# External libraries
import lxml.etree as et



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
        pdbfile = urllib2.urlopen(pdburl).read()
        # If no PDB file is available, a text is now shown with "We're sorry, but ..."
        # Could previously be distinguished by an HTTP error
        if 'sorry' in pdbfile:
            sysexit(5, "No file in PDB format available from wwPDB for the given PDB ID.\n")
    except urllib2.HTTPError:
        sysexit(5, "No file in PDB format available from wwPDB for the given PDB ID.\n")
    return [pdbfile, current_entry]
