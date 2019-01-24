# coding=utf-8
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
test_remote_services.py - Unit Tests for remote services.
"""


import unittest
from plip.modules.webservices import check_pdb_status


class TestPDB(unittest.TestCase):
    """Test PDB Web Service methods"""

    def test_pdb_entry_status(self):
        # 1a0v is an obsolete entry and is replaced by 1y46
        status, current_pdbid = check_pdb_status('1a0v')
        self.assertEqual(status, 'OBSOLETE')
        self.assertEqual(current_pdbid, '1y46')

        # 1vsn is an current entry
        status, current_pdbid = check_pdb_status('1vsn')
        self.assertEqual(status, 'CURRENT')
        self.assertEqual(current_pdbid, '1vsn')

        # xxxx is not an PDB entry
        status, current_pdbid = check_pdb_status('xxxx')
        self.assertEqual(status, 'UNKNOWN')
        self.assertEqual(current_pdbid, 'xxxx')
