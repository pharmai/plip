# coding=utf-8
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
test_basic_functions.py - Unit Tests for basic functionality.
"""

# Python Standard Library
import unittest
import numpy
import random

# Own modules
from plip.modules.preparation import PDBComplex
from plip.modules.supplemental import euclidean3d, vector, vecangle, projection
from plip.modules.supplemental import normalize_vector, cluster_doubles, centroid


class TestLigandSupport(unittest.TestCase):
    """Test for support of different ligands"""

    def test_dna_rna(self):
        """Test if DNA and RNA is correctly processed as ligands"""
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1tf6.pdb')
        # DNA ligand four times consisting of 31 parts (composite)
        self.assertEqual([len(ligand.members) for ligand in tmpmol.ligands].count(31), 4)
        for ligset in [set((x[0] for x in ligand.members)) for ligand in tmpmol.ligands]:
            if len(ligset) == 4:
                # DNA only contains four bases
                self.assertEqual(ligset, set(['DG', 'DC', 'DA', 'DT']))


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


class GeometryTest(unittest.TestCase):
    """Tests for geometrical calculations in PLIP"""

    def vector_magnitude(self, v):
        return numpy.sqrt(sum(x**2 for x in v))

    # noinspection PyUnusedLocal
    def setUp(self):
        """Generate random data for the tests"""
        # Generate two random n-dimensional float vectors, with -100 <= n <= 100 and values 0 <= i <= 1
        dim = random.randint(1, 100)
        self.rnd_vec = [random.uniform(-100, 100) for i in range(dim)]

    def test_euclidean(self):
        """Tests for mathematics.euclidean"""
        # Are the results correct?
        self.assertEqual(euclidean3d([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]), 0)
        self.assertEqual(euclidean3d([2.0, 3.0, 4.0], [2.0, 3.0, 4.0]), 0)
        self.assertEqual(euclidean3d([4.0, 5.0, 6.0], [4.0, 5.0, 8.0]), 2.0)
        # Does the function take vectors or tuples as an input? What about integers?
        self.assertEqual(euclidean3d((4.0, 5.0, 6.0), [4.0, 5.0, 8.0]), 2.0)
        self.assertEqual(euclidean3d((4.0, 5.0, 6.0), (4.0, 5.0, 8.0)), 2.0)
        self.assertEqual(euclidean3d((4, 5, 6), (4.0, 5.0, 8.0)), 2.0)
        # Is the output a float?
        self.assertIsInstance(euclidean3d([2.0, 3.0, 4.0], [2.0, 3.0, 4.0]), float)

    def test_vector(self):
        """Tests for mathematics.vector"""
        # Are the results correct?
        self.assertEqual(list(vector([1, 1, 1], [0, 1, 0])), [-1, 0, -1])
        self.assertEqual(list(vector([0, 0, 10], [0, 0, 4])), [0, 0, -6])
        # Do I get an Numpy Array?
        self.assertIsInstance(vector([1, 1, 1], [0, 1, 0]), numpy.ndarray)
        # Do I get 'None' if the points have different dimensions?
        self.assertEqual(vector([1, 1, 1], [0, 1, 0, 1]), None)

    def test_vecangle(self):
        """Tests for mathematics.vecangle"""
        # Are the results correct?
        self.assertEqual(vecangle([3, 4], [-8, 6], deg=False), numpy.radians(90.0))
        self.assertEqual(vecangle([3, 4], [-8, 6]), 90.0)
        self.assertAlmostEqual(vecangle([-1, -1], [1, 1], deg=False), numpy.pi)
        # Correct if both vectors are equal?
        self.assertEqual(vecangle([3, 3], [3, 3]), 0.0)

    def test_centroid(self):
        """Tests for mathematics.centroid"""
        # Are the results correct?
        self.assertEqual(centroid([[0, 0, 0], [2, 2, 2]]), [1.0, 1.0, 1.0])
        self.assertEqual(centroid([[-5, 1, 2], [10, 2, 2]]), [2.5, 1.5, 2.0])

    def test_normalize_vector(self):
        """Tests for mathematics.normalize_vector"""
        # Are the results correct?
        self.assertAlmostEqual(self.vector_magnitude(normalize_vector(self.rnd_vec)), 1)

    def test_projection(self):
        """Tests for mathematics.projection"""
        # Are the results correct?
        self.assertEqual(projection([-1, 0, 0], [3, 3, 3], [1, 1, 1]), [3, 1, 1])

    def test_cluster_doubles(self):
        """Tests for mathematics.cluster_doubles"""
        # Are the results correct?
        self.assertEqual(set(cluster_doubles([(1, 3), (4, 1), (5, 6), (7, 5)])), {(1, 3, 4), (5, 6, 7)})
