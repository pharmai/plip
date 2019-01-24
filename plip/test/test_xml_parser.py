# coding=utf-8
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
test_xml_parser.py - Unit Tests for XML Parser.
"""


import unittest
from plip.modules.plipxml import PLIPXML


class XMLParserTest(unittest.TestCase):
    """Checks if the XML parser is working correctly"""

    def setUp(self):
        self.px = PLIPXML('./xml/1vsn.report.xml')
        self.bsite = self.px.bsites['NFT:A:283']
        self.smiles = 'CC(C)CC(NC(c1ccc(cc1)c1ccc(cc1)S(N)(=O)=O)C(F)(F)F)C(=O)NCC=N'

    def test_general_information(self):
        """Test if general information is correctly parsed."""
        self.assertEqual(self.px.version, '1.4.2')
        self.assertEqual(self.px.pdbid, '1VSN')
        self.assertFalse(self.px.fixed)
        self.assertEqual(self.px.filename, '1vsn.pdb')
        self.assertEqual(self.px.excluded, [])

    def test_bsite_information(self):
        """Test if the binding site information is correctly parsed."""
        self.assertEqual(self.bsite.pdbid, '1VSN')
        self.assertEqual(self.bsite.uniqueid, '1VSN:NFT:A:283')
        self.assertEqual(self.bsite.hetid, 'NFT')
        self.assertEqual(self.bsite.longname, 'NFT')
        self.assertEqual(self.bsite.ligtype, 'SMALLMOLECULE')
        self.assertEqual(self.bsite.smiles, self.smiles)
        self.assertEqual(self.bsite.members, ['NFT:A:283'])
        self.assertFalse(self.bsite.composite)

        # ligand properties
        self.assertEqual(self.bsite.heavy_atoms, 33)
        self.assertEqual(self.bsite.hbd, 5)
        self.assertEqual(self.bsite.unpaired_hbd, 0)
        self.assertEqual(self.bsite.hba, 7)
        self.assertEqual(self.bsite.unpaired_hba, 2)
        self.assertEqual(self.bsite.hal, 3)
        self.assertEqual(self.bsite.unpaired_hal, 1)
        self.assertEqual(self.bsite.rings, 2)
        self.assertEqual(self.bsite.rotatable_bonds, 12)
        self.assertAlmostEqual(self.bsite.molweight, 484, 0)
        self.assertAlmostEqual(self.bsite.logp, 6, 0)

        # Atom mappings (non-exhaustive test)
        lmap = self.bsite.mappings['pdb_to_smiles']
        self.assertEqual(lmap[1625], 24)
        self.assertEqual(lmap[1649], 33)
        self.assertEqual(lmap[1617], 14)

        # Binding site residues
        self.assertEqual(len(self.bsite.bs_res), 35)

        # Interacting chains
        self.assertEqual(self.bsite.interacting_chains, ['A'])

        # Has Interactions?
        self.assertTrue(self.bsite.has_interactions, True)

    def test_interactions(self):
        """Test if interaction information is correctly parsed."""

        # Hydrophobic Contacts
        self.assertEqual(len(self.bsite.hydrophobics), 4)
        hydrophobic1 = self.bsite.hydrophobics[0]
        self.assertEqual(hydrophobic1.dist, 3.67)
        self.assertEqual(hydrophobic1.resnr, 61)
        self.assertEqual(hydrophobic1.restype, 'ASP')
        self.assertEqual(hydrophobic1.reschain, 'A')
        self.assertEqual(hydrophobic1.ligcarbonidx, 1639)
        self.assertEqual(hydrophobic1.protcarbonidx, 448)
        self.assertEqual(hydrophobic1.ligcoo, (-7.395, 24.225, 6.614))
        self.assertEqual(hydrophobic1.protcoo, (-6.900, 21.561, 9.090))

        # Hydrogen Bonds
        self.assertEqual(len(self.bsite.hbonds), 6)
        hbond1 = self.bsite.hbonds[0]
        self.assertEqual(hbond1.resnr, 19)
        self.assertEqual(hbond1.restype, 'GLN')
        self.assertEqual(hbond1.reschain, 'A')
        self.assertTrue(hbond1.sidechain)
        self.assertEqual(hbond1.dist_h_a, 2.16)
        self.assertEqual(hbond1.dist_d_a, 3.11)
        self.assertEqual(hbond1.don_angle, 160.05)
        self.assertTrue(hbond1.protisdon)
        self.assertEqual(hbond1.donoridx, 153)
        self.assertEqual(hbond1.donortype, 'Nam')
        self.assertEqual(hbond1.acceptoridx, 1649)
        self.assertEqual(hbond1.acceptortype, 'N2')
        self.assertEqual(hbond1.ligcoo, (2.820, 18.145, 6.806))
        self.assertEqual(hbond1.protcoo, (3.976, 15.409, 7.712))

        # Water Bridges
        self.assertEqual(len(self.bsite.wbridges), 1)
        wbridge1 = self.bsite.wbridges[0]
        self.assertEqual(wbridge1.resnr, 159)
        self.assertEqual(wbridge1.restype, 'HIS')
        self.assertEqual(wbridge1.reschain, 'A')
        self.assertEqual(wbridge1.dist_a_w, 3.67)
        self.assertEqual(wbridge1.dist_d_w, 3.13)
        self.assertEqual(wbridge1.don_angle, 126.73)
        self.assertEqual(wbridge1.water_angle, 116.36)
        self.assertTrue(wbridge1.protisdon)
        self.assertEqual(wbridge1.donor_idx, 1210)
        self.assertEqual(wbridge1.donortype, 'Nar')
        self.assertEqual(wbridge1.acceptor_idx, 1649)
        self.assertEqual(wbridge1.acceptortype, 'N2')
        self.assertEqual(wbridge1.ligcoo, (2.820, 18.145, 6.806))
        self.assertEqual(wbridge1.protcoo, (6.401, 19.307, 4.971))
        self.assertEqual(wbridge1.watercoo, (3.860, 18.563, 3.309))

        # Salt Bridges
        self.assertEqual(len(self.bsite.sbridges), 0)

        # Pi stacking
        self.assertEqual(len(self.bsite.pi_stacks), 0)

        # Pi cation interactions
        self.assertEqual(len(self.bsite.pi_cations), 0)

        # Halogen Bonds
        self.assertEqual(len(self.bsite.halogens), 2)
        hal1 = self.bsite.halogens[0]
        self.assertEqual(hal1.resnr, 67)
        self.assertEqual(hal1.restype, 'TYR')
        self.assertEqual(hal1.reschain, 'A')
        self.assertTrue(hal1.sidechain)
        self.assertEqual(hal1.dist, 3.37)
        self.assertEqual(hal1.don_angle, 156.70)
        self.assertEqual(hal1.acc_angle, 100.53)
        self.assertEqual(hal1.don_idx, 1627)
        self.assertEqual(hal1.donortype, 'F')
        self.assertEqual(hal1.acc_idx, 485)
        self.assertEqual(hal1.acceptortype, 'O3')
        self.assertEqual(hal1.ligcoo, (-1.862, 29.303, 4.507))
        self.assertEqual(hal1.protcoo, (-1.005, 26.276, 3.287))

        # Metal complexes
        self.assertEqual(len(self.bsite.metal_complexes), 0)
