# coding=utf-8
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
test_literature_validated.py - Unit Tests for literature-validated cases.
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


class MetalCoordinationTest(unittest.TestCase):
    """Checks predictions against literature-validated interactions for metal coordination."""

    ###############################################
    # Literature-validated cases from publication #
    ###############################################

    def test_test(self):
        """Binding of anti-Alzheimer drug E2020 to acetylcholinesterase from Torpedo californica (1eve)
        Reference: Chakrabarti et al. Geometry of nonbonded interactions involving planar groups in proteins. (2007)
        """
        # #@todo Modify this test and implement other unit tests
        # #@todo Add tests for coordination numbers and geometry
        # #@todo E.g. 1GMW (Cu)
        # #@todo e.g. 101M (Heme, Fe); 4HRO, 1A01 (both Fe)
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1eve.pdb')
        s = tmpmol.interaction_sets['E20-A-2001']
        # Aromatic stacking with Trp84 and Trp279
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({84, 279}.issubset(pistackres))
        # Pi-Cation interaction of Phe330 with ligand
        pication = {pication.resnr for pication in s.pication_paro}
        self.assertTrue({330}.issubset(pication))
