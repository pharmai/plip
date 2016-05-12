# coding=utf-8
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
test_remote_services.py - Unit Tests for remote services.
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


        #print "Latest PDB ---", "passed" if get_latest_pdb("1a0v") == "1y46" else "not passed"
#print "Resolution ---", "passed" if get_resolution("4HHB") == 1.74 else "not passed"
