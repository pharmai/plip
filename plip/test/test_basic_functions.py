# coding=utf-8
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
test_basic_functions.py - Unit Tests for basic functionality.
Copyright 2014-2015 Sebastian Salentin

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
from plip.modules.preparation import PDBComplex


class TestMapping(unittest.TestCase):
    """Test"""

    def test_ids(self):
        """Test if the atom IDs are correctly mapped from internal to original PDB."""
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1vsn.pdb')
        bsid = 'NFT:A:283'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        for contact in s.hydrophobic_contacts:
            if contact.restype == 'ALA' and contact.resnr == 133:
                self.assertEqual(contact.ligatom_orig_idx, 1636)
                self.assertEqual(contact.bsatom_orig_idx, 994)
            if contact.restype == 'ASP' and contact.resnr == 61:
                self.assertEqual(contact.ligatom_orig_idx, 1639)
                self.assertEqual(contact.bsatom_orig_idx, 448)
        for contact in s.hbonds_ldon + s.hbonds_pdon:
            if contact.restype == 'GLN' and contact.resnr == 19:
                self.assertEqual(contact.a_orig_idx, 1649)
                self.assertEqual(contact.d_orig_idx, 153)
            if contact.restype == 'CYS' and contact.resnr == 25:
                self.assertEqual(contact.a_orig_idx, 1649)
                self.assertEqual(contact.d_orig_idx, 183)
            if contact.restype == 'ASN' and contact.resnr == 158:
                self.assertEqual(contact.d_orig_idx, 1629)
                self.assertEqual(contact.a_orig_idx, 1199)
        for contact in s.halogen_bonds:
            if contact.restype == 'TYR' and contact.resnr == 67:
                self.assertEqual(contact.don.x_orig_idx, 1627)
                self.assertEqual(contact.acc.o_orig_idx, 485)
            if contact.restype == 'LEU' and contact.resnr == 157:
                self.assertEqual(contact.don.x_orig_idx, 1628)
                self.assertEqual(contact.acc.o_orig_idx, 1191)
