import sys
from urllib.error import HTTPError
from urllib.request import urlopen

import lxml.etree as et

from plip.basic import logger

logger = logger.get_logger()


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
    # logger.info(f'checking status of PDB-ID {pdbid}')
    # @todo re-implement state check with ew RCSB API, see https://www.rcsb.org/news?year=2020&article=5eb18ccfd62245129947212a&feature=true
    # state, current_entry = check_pdb_status(pdbid)  # Get state and current PDB ID
    #
    # if state == 'OBSOLETE':
    #     logger.info(f'entry is obsolete, getting {current_entry} instead')
    # elif state == 'CURRENT':
    #     logger.info('entry is up-to-date')
    # elif state == 'UNKNOWN':
    #     logger.error('invalid PDB-ID (entry does not exist on PDB server)')
    #     sys.exit(1)
    logger.info('downloading file from PDB')
    # get URL for current entry
    # @todo needs update to react properly on response codes of RCSB servers
    pdburl = f'https://files.rcsb.org/download/{pdbid}.pdb'
    try:
        pdbfile = urlopen(pdburl).read().decode()
        # If no PDB file is available, a text is now shown with "We're sorry, but ..."
        # Could previously be distinguished by an HTTP error
        if 'sorry' in pdbfile:
            logger.error('no file in PDB format available from wwPDB for the given PDB ID.')
            sys.exit(1)
    except HTTPError:
        logger.error('no file in PDB format available from wwPDB for the given PDB ID')
        sys.exit(1)
    return [pdbfile, pdbid]
