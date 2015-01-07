# coding=utf-8
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
Module for testing the predictions with interactions from literature.
Copyright (C) 2014  Sebastian Salentin

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
"""


import unittest
from plip.modules.preparation import PDBComplex


class LiteratureValidatedTest(unittest.TestCase):
    """Checks predictions against literature-validated interactions"""

    def test_1eve(self):
        """Binding of anti-Alzheimer drug E2020 to acetylcholinesterase from Torpedo californica (1eve)
        Reference: Chakrabarti et al. Geometry of nonbonded interactions involving planar groups in proteins. (2007)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1eve.pdb')
        s = tmpmol.interaction_sets['E20-A-2001']
        # Aromatic stacking with Trp84 and Trp279
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({84, 279}.issubset(pistackres))
        # Pi-Cation interaction of Phe330 with ligand
        pication = {pication.resnr for pication in s.pication_paro}
        self.assertTrue({330}.issubset(pication))

    def test_1h2t(self):
        """Binding of methylated guanosine to heterodimeric nuclear-cap binding complex (1h2t)
        Reference: Chakrabarti et al. Geometry of nonbonded interactions involving planar groups in proteins. (2007)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1h2t.pdb')
        s = tmpmol.interaction_sets['7MG-Z-1152']
        # Sandwiched pi-stacking involving Tyr20 and Tyr43
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({20, 43}.issubset(pistackres))
        # Hydrogen bond with R112
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({112}.issubset(hbonds))
        # Salt bridge with D116
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_pneg}
        self.assertTrue({116}.issubset(saltb))

    def test_3pxf(self):
        """Binding of ANS to CDK2 (3pxf)
        Reference: Betzi et al. Discovery of a potential allosteric ligand binding site in CDK2 (2012)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/3pxf.pdb')
        s = tmpmol.interaction_sets['2AN-A-305']
        # Hydrogen bonding of Asp145 and Phe146
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({145, 146}.issubset(hbonds))
        # Salt bridge by Lys33 to sulfonate group
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_lneg}
        self.assertTrue({33}.issubset(saltb))
        # Naphtalene positioned between Leu55 and Lys56, indicating hydrophobic interactions
        hydroph = {hydroph.resnr for hydroph in s.hydrophobic_contacts}
        self.assertTrue({55, 56}.issubset(hydroph))

        s = tmpmol.interaction_sets['2AN-A-304']
        # Salt bridges to sulfonate group by Lys56 and His71
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_lneg}
        self.assertTrue({56, 71}.issubset(saltb))
        # Napthalene with hydrophobic interactions to Ile52 and Leu76
        hydroph = {hydroph.resnr for hydroph in s.hydrophobic_contacts}
        self.assertTrue({52, 76}.issubset(hydroph))

    def test_2reg(self):
        """Binding of choline to ChoX (2reg)
        Reference: Oswald et al. Crystal structures of the choline/acetylcholine substrate-binding protein ChoX
        from Sinorhizobium meliloti in the liganded and unliganded-closed states. (2008)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/2reg.pdb')
        s = tmpmol.interaction_sets['CHT-A-1']
        # Cation-pi interactions with Trp43, Trp90, Trp205, and Tyr119
        picat = {pication.resnr for pication in s.pication_paro}
        self.assertEqual({43, 90, 205, 119}, picat)
        # Saltbridge to Asp45
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_pneg}
        self.assertEqual({45}, saltb)

    def test_1osn(self):
        """Binding of VZV-tk to BVDU-MP (2reg)
        Reference: Bird et al. Crystal structures of Varicella Zoster Virus Thyrimidine Kinase. (2003)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1osn.pdb')
        s = tmpmol.interaction_sets['BVP-A-500']
        # Sandwiched pi-stacking involving Phe93 and Phe139
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({93, 139}.issubset(pistackres))
        # Hydrogen bonding of Gln90
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({90}.issubset(hbonds))

    def test_2w0s(self):
        """Binding of Vacc-TK to TDP (2w0s)
        Reference: Caillat et al. Crystal structure of poxvirus thymidylate kinase: An unexpected dimerization
        has implications for antiviral therapy (2008)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/2w0s.pdb')
        s = tmpmol.interaction_sets['BVP-B-1207']
        # Hydrogen bonding of Tyr101 and Arg72
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({101}.issubset(hbonds))
        # Halogen Bonding of Asn65
        halogens = {halogen.resnr for halogen in s.halogen_bonds}
        self.assertTrue({65}.issubset(halogens))
        # pi-stacking interaction with Phe68
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({68}.issubset(pistackres))
        # Saltbridge to Arg41 and Arg93
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_lneg}
        self.assertTrue({41, 93}.issubset(saltb))

    def test_1vsn(self):
        """Binding of NFT to Cathepsin K (1vsn)
        Reference: Li et al. Identification of a potent and selective non-basic cathepsin K inhibitor. (2006)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1vsn.pdb')
        s = tmpmol.interaction_sets['NFT-A-283']
        # Hydrogen bonding to Gly66
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({66}.issubset(hbonds))

    def test_3shy(self):
        """Binding of PDE5 to halogenated 1 (3shy)
        Reference: Xu et al. Utilization of halogen bond in lead optimization: a case study of rational design
        of potent phosphodiesterase type 5 (PDE5) inhibitors. (2006)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/3shy.pdb')
        s = tmpmol.interaction_sets['5FO-A-1']
        # Halogen Bonding of Tyr612
        halogens = {halogen.resnr for halogen in s.halogen_bonds}
        self.assertTrue({612}.issubset(halogens))
        # Hydrogen bonding to Gln817
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({817}.issubset(hbonds))
        # pi-stacking interaction with Phe820
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({820}.issubset(pistackres))

    def test_1p5e(self):
        """Binding of (1p5e)
        Reference: . (2006)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/3shy.pdb')
        s = tmpmol.interaction_sets['5FO-A-1']
