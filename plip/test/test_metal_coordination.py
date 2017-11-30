# coding=utf-8
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
test_metal_coordination.py - Unit Tests for Metal Coordination.
"""


import unittest
from plip.modules.preparation import PDBComplex


class MetalCoordinationTest(unittest.TestCase):
    """Checks predictions against literature-validated interactions for metal coordination."""

    ###############################################
    # Literature-validated cases from publication #
    ###############################################

    def test_1rmd(self):
        """Zinc binding sites in RAG1 dimerization domain (1rmd)
        Reference: Harding. The architecture of metal coordination groups in proteins. (2004), Fig. 1a
        """

        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1rmd.pdb')
        bsid = 'ZN:A:119'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Coordination by three cysteines and one histidine of the protein
        metalres = [mres.restype for mres in s.metal_complexes]
        self.assertEqual(metalres.count('CYS'), 3)
        self.assertEqual(metalres.count('HIS'), 1)
        # Zn atom with tetrahedral geometry (coordination number 4)
        self.assertEqual(s.metal_complexes[0].coordination_num, 4)
        self.assertEqual(s.metal_complexes[0].geometry, 'tetrahedral')

    def test_1rla(self):
        """Rat liver arginase, a binuclear manganese metalloenzyme (1rmd)
        Reference: Harding. The architecture of metal coordination groups in proteins. (2004), Fig. 1b
        """

        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1rla.pdb')
        bsid = 'MN:A:500'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Coordination by one histidine, three aspartic acid residues, and one water molecule
        metalres = [mres.restype for mres in s.metal_complexes]
        self.assertEqual(metalres.count('HIS'), 1)
        self.assertEqual(metalres.count('ASP'), 3)
        self.assertEqual(metalres.count('HOH'), 1)
        # Mn atom with square pyramidal geometry (coordination number 5)
        self.assertEqual(s.metal_complexes[0].coordination_num, 5)
        self.assertEqual(s.metal_complexes[0].geometry, 'square.pyramidal')

    def test_1het(self):
        """Liver alcohol deshydrogenase (1het)
        Reference: Harding. The architecture of metal coordination groups in proteins. (2004), Fig. 2
        """

        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1het.pdb')
        bsid = 'ZN:A:401'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Coordination by four cysteines
        metalres = [mres.restype + str(mres.resnr) for mres in s.metal_complexes]
        self.assertEqual(set(metalres), {'CYS97', 'CYS100', 'CYS103', 'CYS111'})
        # Zn atom with tetrahedral geometry (coordination number 4)
        self.assertEqual(s.metal_complexes[0].coordination_num, 4)
        self.assertEqual(s.metal_complexes[0].geometry, 'tetrahedral')

    def test_1vfy(self):
        """Phosphatidylinositol-3-phosphate binding FYVE domain of VPS27P protein (1vfy)
        Reference: Harding. The architecture of metal coordination groups in proteins. (2004), Fig. 5
        """

        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1vfy.pdb')
        bsid = 'ZN:A:300'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Coordination by four cysteines
        metalres = [mres.restype for mres in s.metal_complexes]
        self.assertEqual(set(metalres), {'CYS'})
        # Zn atom with tetrahedral geometry (coordination number 4)
        self.assertEqual(s.metal_complexes[0].coordination_num, 4)
        self.assertEqual(s.metal_complexes[0].geometry, 'tetrahedral')

    def test_2pvb(self):
        """Pike parvalbumin binding calcium (2pvb)
        Reference: Harding. The architecture of metal coordination groups in proteins. (2004), Fig. 6
        """

        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/2pvb.pdb')
        bsid = 'CA:A:110'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Ca atom with square pyramidal geometry (coordination number 5)
        self.assertEqual(s.metal_complexes[0].coordination_num, 5)
        self.assertEqual(s.metal_complexes[0].geometry, 'square.pyramidal')

    def test_2q8q(self):
        """Crystal Structure of S. aureus IsdE complexed with heme (2q8q)
        Reference: Grigg et al. Heme coordination by Staphylococcus aureus IsdE. (2007)
        """

        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/2q8q.pdb')
        bsid = 'HEM:A:300'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Coordination by four nitrogens of heme itself and one additional histidine from the protein
        metalres = [mres.restype for mres in s.metal_complexes]
        self.assertEqual(metalres.count('HEM'), 4)
        self.assertEqual(metalres.count('HIS'), 1)
        # Fe atom with square pyramidal geometry (coordination number 5)
        self.assertEqual(s.metal_complexes[0].coordination_num, 5)
        self.assertEqual(s.metal_complexes[0].geometry, 'square.pyramidal')
