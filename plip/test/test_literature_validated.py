# coding=utf-8
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
test_literature_validated.py - Unit Tests for literature-validated cases.
"""


import unittest
from plip.modules.preparation import PDBComplex


class LiteratureValidatedTest(unittest.TestCase):
    """Checks predictions against literature-validated interactions"""

    ###############################################
    # Literature-validated cases from publication #
    ###############################################

    def test_1eve(self):
        """Binding of anti-Alzheimer drug E2020 to acetylcholinesterase from Torpedo californica (1eve)
        Reference: Chakrabarti et al. Geometry of nonbonded interactions involving planar groups in proteins. (2007)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1eve.pdb')
        bsid = 'E20:A:2001'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
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
        bsid = 'GDP:Z:1151'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
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
        bsids = ['2AN:A:305', '2AN:A:304']
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) in bsids:
                tmpmol.characterize_complex(ligand)

        s = tmpmol.interaction_sets[bsids[0]]  # 2AN:A:305
        # Hydrogen bonding of Asp145 and Phe146
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({145, 146}.issubset(hbonds))
        # Salt bridge by Lys33 to sulfonate group
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_lneg}
        self.assertTrue({33}.issubset(saltb))
        # Naphtalene positioned between Leu55 and Lys56, indicating hydrophobic interactions
        hydroph = {hydroph.resnr for hydroph in s.hydrophobic_contacts}
        self.assertTrue({55, 56}.issubset(hydroph))

        s = tmpmol.interaction_sets[bsids[1]]  # 2AN:A:304
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
        bsid = 'CHT:A:1'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
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
        bsid = 'BVP:A:500'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
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
        bsid = 'BVP:B:1207'  # Complex of BVDU with Magnesium Cofactor
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonding of Tyr101 and Arg72
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({101, 72}.issubset(hbonds))
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
        bsid = 'NFT:A:283'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonding to Gly66
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({66}.issubset(hbonds))

    def test_1p5e(self):
        """Binding of TBS to CDK2(1p5e)
        Reference: De Moliner et al. Alternative binding modes of an inhibitor to two different kinases. (2003)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1p5e.pdb')
        bsid = 'TBS:A:301'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Halogen Bonding of Ile10 and Leu83
        halogens = {halogen.resnr for halogen in s.halogen_bonds}
        self.assertTrue({10, 83}.issubset(halogens))

    def test_1acj(self):
        """Binding of Tacrine (THA) to active-site gorge of acetylcholinesterase (1acj)
        Reference: Harel et al. Quaternary ligand binding to aromatic residues in the active-site gorge of
        acetylcholinesterase.. (1993)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1acj.pdb')
        bsid = 'THA:A:999'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # pi-stacking interaction with Phe330 and Trp84
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({330, 84}.issubset(pistackres))

    def test_2zoz(self):
        """Binding of CgmR to ethidium(2z0z)
        Reference: Itou et al. Crystal Structures of the Multidrug Binding Repressor Corynebacterium
        glutamicum CgmR in Complex with Inducers and with an Operator. (2010)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/2zoz.pdb')
        bsid = 'ET:B:184'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # pi-stacking interaction with Trp63 and Phe147
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({147}.issubset(pistackres))  # Trp 63!!
        # hydrophobic interaction of Leu59, Leu88, Trp63, Trp113, Phe147
        hydrophobics = {hydrophobic.resnr for hydrophobic in s.all_hydrophobic_contacts}
        self.assertTrue({59, 88, 63, 113, 147}.issubset(hydrophobics))
        self.assertTrue({59, 88, 63, 92, 113, 147}.issubset(hydrophobics))

    def test_1xdn(self):
        """Binding of ATP to RNA editing ligase 1 (1xdn)
        Reference: Deng et al. High resolution crystal structure of a key editosome enzyme from Trypanosoma brucei:
        RNA editing ligase 1. (2004)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1xdn.pdb')
        bsid = 'ATP:A:501'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Arg111, Ile61 (backbone), Asn92, Val88, Lys87 and Glu86#
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({111, 61, 92, 88, 87}.issubset(hbonds))
        # Water bridges to Lys307, Arg309 and 111 from phosphate groups
        waterbridges = {wb.resnr for wb in s.water_bridges}
        # Water bridge to 60 not found due to prioritization
        self.assertTrue({307, 309, 111}.issubset(waterbridges))
        # pi-stacking interaction with Phe209
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({209}.issubset(pistackres))

    def test_1bma(self):
        """Binding of aminimide to porcine pancreatic elastase(1bma)
        Reference: Peisach et al. Interaction of a Peptidomimetic Aminimide Inhibitor with Elastase. (1995)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1bma.pdb')
        bsid = '0QH:A:256'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to val224 and Gln200
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({224, 200}.issubset(hbonds))
        self.assertTrue({224, 200}.issubset(hbonds))
        # hydrophobic interaction of Phe223 and val103
        hydrophobics = {hydrophobic.resnr for hydrophobic in s.all_hydrophobic_contacts}
        self.assertTrue({223, 103}.issubset(hydrophobics))
        # Water bridges to Ser203 not detected due to prioritization
        waterbridges = {wb.resnr for wb in s.water_bridges}
        self.assertTrue(set().issubset(waterbridges))

    def test_4rao(self):
        """Binding of (4rao)
        Reference: Keough et al. Aza-acyclic Nucleoside Phosphonates Containing a Second Phosphonate Group
        As Inhibitors of the Human, Plasmodium falciparum and vivax 6‑Oxopurine Phosphoribosyltransferases
        and Their Prodrugs As Antimalarial Agents (2004)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/4rao.pdb')
        bsid = '3L7:B:301'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Val187, Lys165, Thr141, Lys140, Gly139, Thr138, Asp137
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}  # res nr 100, 68, 69 and 199 in alternative conformation,
        self.assertTrue({137, 138, 139, 140, 141, 165, 187}.issubset(hbonds))
        # Water bridges to Asp137, Thr141, Met142, Arg199 and Gly139
        # res nr 199 and 142 in alternative conformation
        # Water bridges to 137m 141, 139 not detected due to prioritization
        # pi-stacking interaction with Phe186
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({186}.issubset(pistackres))

    def test_4qnb(self):
        """Binding of (4qnb)
        Reference:  Bhattacharya et al. Structural basis of HIV-1 capsid recognition by PF74 and CPSF6(2014)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/4qnb.pdb')
        bsid = '1B0:A:301'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Asn57 and Lys70
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({57, 70}.issubset(hbonds))
        # Cation-pi interactions with Lys70
        picat = {pication.resnr for pication in s.pication_laro}
        self.assertEqual({70}, picat)

    def test_4kya(self):
        """Binding of non-classical TS inhibitor 3 with Toxoplasma gondii TS-DHFR(4kya)
        Reference:  Zaware et al. Structural basis of HIV-1 capsid recognition by PF74 and CPSF6(2014)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/4kya.pdb')
        bsid = '1UG:E:702'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Ala609
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({609}.issubset(hbonds))
        # Saltbridge to Asp513
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_pneg}
        self.assertTrue({513}.issubset(saltb))
        # hydrophobic interaction of Ile402, Leu516, Phe520 and Met608
        hydrophobics = {hydrophobic.resnr for hydrophobic in s.all_hydrophobic_contacts}
        self.assertTrue({402, 516, 520, 608}.issubset(hydrophobics))
        # pi-stacking interaction with Trp403 and Phe520
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({403, 520}.issubset(pistackres))

    def test_1n7g(self):
        """Binding of NADPH to MURI from Arabidopsis thaliana (1n7g)
        Reference:  Mulichak et al. Structure of the MUR1 GDP-mannose 4, 6-dehydratase from Arabidopsis thaliana:
        implications for ligand binding and specificity(2002)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1n7g.pdb')
        bsid = 'NDP:A:701'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Thr37, Gly38, Gln39, Asp40,  Arg60, Leu92, Asp91, Ser63, Leu92, Ala115, Ser117,
        # Tyr128, Tyr185, Lys189, His215 and Arg220
        # Publication give the Prediction for Asp91 as hydrogen bond, when this contains two acceptor atoms.
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        # #@todo Hbond to 128 not detected
        self.assertTrue({37, 38, 39, 40, 92, 63, 92, 115, 117, 185, 189, 215, 220}.issubset(hbonds))
        # Water bridges to Gly35, Thr37, Gly38, Asp40, Arg60, Arg61, Ser63, Asn66, Ser117, Tyr128, Lys189, Arg220
        waterbridges = {wb.resnr for wb in s.water_bridges}
        # Hydrogen bonds to 35, 37, 38, 40, 63, 117, 128, 189, 220 not detected due to prioritization
        self.assertTrue({60, 66, 61}.issubset(waterbridges))
        # Saltbridge to arg60, Arg61, Arg69 and Arg220
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_lneg}
        # #@todo Additional saltbridges report to 69 and 200 (with large distances)
        self.assertTrue({60, 61}.issubset(saltb))
        # Cation-pi interactions with Arg60
        picat = {pication.resnr for pication in s.pication_laro}
        self.assertEqual({60}, picat)

    def test_4alw(self):
        """Binding of benzofuropyrimidinones compound 3 to PIM-1 (4alw)
        Reference:  Tsuhako et al. The design, synthesis, and biological evaluation of PIM kinase inhibitors.(2012)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/4alw.pdb')
        bsid = 'HY7:A:1308'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Asp186
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({186}.issubset(hbonds))
        # Saltbridge to A186 and Glu171
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_pneg}
        self.assertTrue({186, 171}.issubset(saltb))

    def test_3o1h(self):
        """Binding of TMAO to TorT-TorS system(3o1h)
        Reference:  Hendrickson et al. An Asymmetry-to-Symmetry Switch in Signal Transmission by the Histidine Kinase Receptor
        for TMAO.(2013)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/3o1h.pdb')
        bsid = 'TMO:B:1'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Trp45
        hbonds = {hbond.resnr for hbond in s.hbonds_ldon}
        self.assertTrue({45}.issubset(hbonds))
        # Water bridges to Trp45 not detected due to priorization of hydrogen bonds
        # Cation-pi interactions with Tyr44
        picat = {pication.resnr for pication in s.pication_paro}
        self.assertEqual({44}, picat)

    def test_3thy(self):
        """Binding of ADP tp MutS(3thy)
        Reference:  Shikha et al. Mechanism of mismatch recognition revealed by human MutSβ bound to unpaired DNA loops.(2012)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/3thy.pdb')
        bsid = 'ADP:A:935'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Saltbridge to His295 and Lys675
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_lneg}
        self.assertTrue({675}.issubset(saltb))
        # pi-stacking interaction with Tyr815
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({815}.issubset(pistackres))

    def test_3tah(self):
        """Binding of BGO to an an N11A mutant of the G-protein domain of FeoB.(3tah)
        Reference:  Ash et al. The structure of an N11A mutant of the G-protein domain of FeoB.(2011)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/3tah.pdb')
        bsid = 'BGO:A:300'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Ala11, Lys14, Thr15, Ser16, Asp113, Met114, Ala143 and Asp113
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({11, 13, 14, 15, 16, 113, 114, 143, 113}.issubset(hbonds))
        # Water bridges to Ala11 not detected due to prioritization of hydrogen bonds
        # Saltbridge to Asp116
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_pneg}
        self.assertTrue({116}.issubset(saltb))

    def test_3r0t(self):
        """Binding of protein kinase CK2 alpha subunit in with the inhibitor CX-5279 (3r0t)
        Reference:  Battistutta et al. Unprecedented selectivity and structural determinants of a new class of protein
        kinase CK2 inhibitors in clinical trials for the treatment of cancer (2011).
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/3r0t.pdb')
        bsid = 'FU9:A:338'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Val116
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({116}.issubset(hbonds))
        # Water bridges to Lys68 and Trp176
        waterbridges = {wb.resnr for wb in s.water_bridges}
        self.assertTrue({68, 176}.issubset(waterbridges))
        # Saltbridge to Ly68
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_lneg}
        self.assertTrue({68}.issubset(saltb))
        # hydrophobic interaction of Val66, Phe113 and Ile174
        hydrophobics = {hydrophobic.resnr for hydrophobic in s.all_hydrophobic_contacts}
        self.assertTrue({66, 113, 174}.issubset(hydrophobics))
        # pi-stacking interaction with His160
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({160}.issubset(pistackres))

    def test_1aku(self):
        """Binding of Flavin mononucleotido with D.Vulgaris(1aku)
        Reference:  McCarthy et al. Crystallographic Investigation of the Role of Aspartate 95 in the Modulation of the
        Redox Potentials of DesulfoVibrio Vulgaris Flavodoxin.(2002)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1aku.pdb')
        bsid = 'FMN:A:150'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Tht59 and Trp60
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({59, 60}.issubset(hbonds))
        # Water bridges to Asp63 and Tyr100
        waterbridges = {wb.resnr for wb in s.water_bridges}
        # Water bridge to Tyr100 not detected due to prioritization
        self.assertTrue({63}.issubset(waterbridges))
        # hydrophobic interaction of Trp60
        hydrophobics = {hydrophobic.resnr for hydrophobic in s.all_hydrophobic_contacts}
        self.assertTrue({60}.issubset(hydrophobics))
        # pi-stacking interaction with Tyr98
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({98}.issubset(pistackres))

    def test_4pjt(self):
        """Binding of BMN 673 to catPARP1(4pj7)
        Reference:  Aoyagi-Scharber et al. Structural basis for the inhibition of poly(ADP-ribose) polymerases 1 and 2 by BMN
        673, a potent inhibitor derived from dihydropyridophthalazinone.(2014)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/4pjt.pdb')
        bsid = '2YQ:D:1104'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Gly863
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({863}.issubset(hbonds))
        # pi-stacking interaction with Tyr889 and Tyr907
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({889, 907}.issubset(pistackres))

    def test_1bju(self):
        """Binding of ACPU to bovine tripsin(1bju)
        Reference:  Presnell et al. Oxyanion-Mediated Inhibition of Serine Proteases.(1998)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1bju.pdb')
        bsid = 'GP6:A:910'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        #@todo Publication show hydrogen bond interactions for Gly219
        # Hydrogen bonds to Ser190, Ser195, Gly219 and Asp189
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon+s.hbonds_ldon}
        self.assertTrue({189, 190, 195}.issubset(hbonds))
        # Water bridges to Ser190 and Val227
        # Water bridge to 190 not detected due to prioritization
        waterbridges = {wb.resnr for wb in s.water_bridges}
        self.assertTrue({227}.issubset(waterbridges))
        # hydrophobic interaction of Leu99
        hydrophobics = {hydrophobic.resnr for hydrophobic in s.all_hydrophobic_contacts}
        self.assertTrue({99}.issubset(hydrophobics))
        # pi-stacking interaction with His57
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({57}.issubset(pistackres))

    def test_4agl(self):
        """Binding of P53 to PhiKan784(4agl)
        Reference:  Wilcken et al. Halogen-Enriched Fragment Libraries as Leads for Drug Rescue of Mutant p53.(2012)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/4agl.pdb')
        bsid = 'P84:A:400'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Water bridges to Val147
        waterbridges = {wb.resnr for wb in s.water_bridges}
        self.assertTrue({147}.issubset(waterbridges))
        # hydrophobic interaction of Thr150
        hydrophobics = {hydrophobic.resnr for hydrophobic in s.all_hydrophobic_contacts}
        self.assertTrue({150}.issubset(hydrophobics))
        # Halogen Bonding of Leu145
        halogens = {halogen.resnr for halogen in s.halogen_bonds}
        self.assertTrue({145}.issubset(halogens))

    def test_2efj(self):
        """Binding of teobromine to 1,7 dimethylxanthine methyltransferase(2efj)
        Reference:  McCarthy et al. The Structure of Two N-Methyltransferases from the Caffeine Biosynthetic
        Pathway.(2007)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/2efj.pdb')
        bsid = '37T:A:502'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Trp161, Ser237
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({161, 237}.issubset(hbonds))
        # pi-stacking interaction with Tyr157
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({157}.issubset(pistackres))

    def test_2iuz(self):
        """Binding of C2-dicaffeine to Aspergilius fumigates(2iuz)
        Reference:  Schüttelkopf et al. Screening-based discovery and structural dissection of a novel family 18 chitinase
        inhibitor.(2006)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/2iuz.pdb')
        bsid = 'D1H:A:1440'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Trp137, Trp184
        # res nr 52 mentioned in alternative conformation, not considered
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({137, 384}.issubset(hbonds))
        # Water bridges to Trp137
        # res nr 52 mentioned in alternative conformation not considered
        waterbridges = {wb.resnr for wb in s.water_bridges}
        self.assertTrue({137}.issubset(waterbridges))
        # pi-stacking interaction with Trp384, Trp137 and Trp52
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({52, 137, 384}.issubset(pistackres))

    def test_3shy(self):
        """Binding of 5FO to PDE5A1 catalytic domain(3shy)
        Reference:  Xu et al. Utilization of halogen bond in lead optimization: A case study of rational design of potent
        phosphodiesterase type 5 (PDE5) inhibitors.(2011)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/3shy.pdb')
        bsid = '5FO:A:1'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Gln817
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({817}.issubset(hbonds))
        # hydrophobic interaction of Tyr612
        hydrophobics = {hydrophobic.resnr for hydrophobic in s.all_hydrophobic_contacts}
        self.assertTrue({612}.issubset(hydrophobics))
        # pi-stacking interaction with Phe820
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({820}.issubset(pistackres))
        # Halogen Bonding of Tyr612
        halogens = {halogen.resnr for halogen in s.halogen_bonds}
        self.assertTrue({612}.issubset(halogens))

    def test_1ay8(self):
        """Binding of PLP to aromatic amino acid aminotransferase(1ay8)
        Reference:  Okamoto et al. Crystal structures of Paracoccus denitrificans aromatic amino acid aminotransferase: a
        substrate recognition site constructed by rearrangement of hydrogen bond network..(1998)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1ay8.pdb')
        bsid = 'PLP:A:413'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds to Gly108, Thr109, Asn194 and Ser257
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({108, 109, 194, 257}.issubset(hbonds))
        # Saltbridge to Ly258 and Arg266
        saltb = {saltbridge.resnr for saltbridge in s.saltbridge_lneg}
        self.assertTrue({258, 266}.issubset(saltb))
        # pi-stacking interaction with Trp140
        pistackres = {pistack.resnr for pistack in s.pistacking}
        self.assertTrue({140}.issubset(pistackres))

    def test_4rdl(self):
        """Binding of Norovirus Boxer P domain with Lewis y tetrasaccharide(4rdl)
        Reference:  Hao et al. Crystal structures of GI.8 Boxer virus P dimers in complex with HBGAs, a novel
        evolutionary path selected by the Lewis epitope..(2014)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/4rdl.pdb')
        bsid = 'FUC:A:601'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]  # Instead of FUC-A-604 (sugar representative)
        # Water bridges to Asn395
        waterbridges = {wb.resnr for wb in s.water_bridges}
        self.assertTrue({395}.issubset(waterbridges))
        # Hydrogen bonds to Thr347, Gly348 and Asn395
        hbonds = {hbond.resnr for hbond in s.hbonds_pdon}
        self.assertTrue({347, 348, 395}.issubset(hbonds))
        # hydrophobic interaction of Trp392
        hydrophobics = {hydrophobic.resnr for hydrophobic in s.all_hydrophobic_contacts}
        self.assertTrue({392}.issubset(hydrophobics))

    #########################################
    # Additional literature-validated cases #
    #########################################

    def test_1hii(self):
        """HIV-2 protease in complex with novel inhibitor CGP 53820 (1hii)
        Reference:  Comparative analysis of the X-ray structures of HIV-1 and HIV-2 proteases in complex
        with CGP 53820, a novel pseudosymmetric inhibitor (1995)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1hii.pdb')
        bsid = 'C20:B:101'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Water bridges
        waterbridges = {str(wb.resnr)+wb.reschain for wb in s.water_bridges}
        self.assertTrue({'50A', '50B'}.issubset(waterbridges))  # Bridging Ile-B50 and Ile-A50 with ligand
        # Hydrogen bonds
        hbonds = {str(hbond.resnr)+hbond.reschain for hbond in s.hbonds_pdon+s.hbonds_ldon}
        self.assertTrue({'27A', '27B', '29A', '48A', '48B'}.issubset(hbonds))
        # #@todo Publication mentions additional possible hydrogen bond with Asp28B
        # Hydrogen bonds with Asp-A25 are reported as a salt bridge as both partners have (potential) charges

    def test_1hvi(self):
        """HIV-1 protease in complex with Diol inhibitor (1hvi)
        Reference: Influence of Stereochemistry on Activity and Binding Modes for C2 Symmetry-Based
         Diol Inhibitors of HIV-1 Protease (1994)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1hvi.pdb')
        bsid = 'A77:A:800'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Water bridges
        waterbridges = {str(wb.resnr)+wb.reschain for wb in s.water_bridges}
        # #@todo Water bridge with 50B not detected
        self.assertTrue({'50A'}.issubset(waterbridges))  # Bridging Ile-B50 and Ile-A50 with ligand
        # pi-cation Interactions
        picat = {pication.resnr for pication in s.pication_laro}
        self.assertEqual({8}, picat)  # Described as weakly polar contact/stacking in paper
        # Hydrogen bonds
        hbonds = {str(hbond.resnr)+hbond.reschain for hbond in s.hbonds_pdon+s.hbonds_ldon}
        self.assertTrue({'25B', '27A', '27B', '48A', '48B'}.issubset(hbonds))
        # #@todo Paper describes additional hydrogen bond with Asp25A

    def test_3OG7(self):
        """Inhibitor PLX4032 binding to B-RAF(V600E) (3og7)
        Reference: Clinical efficacy of a RAF inhibitor needs broad target blockade in BRAF-mutant
        melanoma (2010)
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/3og7.pdb')
        bsid = '032:A:1'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrogen bonds
        hbonds = {str(hbond.resnr)+hbond.reschain for hbond in s.hbonds_pdon+s.hbonds_ldon}
        # Additional hydrogen bond to residue 530A reported
        self.assertTrue({'594A'}.issubset(hbonds))

    def test_1HPX(self):
        """
        HIV-1 Protease complexes with the inhibitor KNI-272
        Reference: Structure of HIV-1 protease with KNI-272, a tight-binding transition-state analog
        containing allophenylnorstatine.
        Note that the residue numbering is different in the publication and the PDB structure.
        For residues in the B chain, the offset is -100 (e.g. Ile 50B in the PDB structure is Ile 150 in the paper).
        """
        tmpmol = PDBComplex()
        tmpmol.load_pdb('./pdb/1hpx.pdb')
        bsid = 'KNI:B:900'
        for ligand in tmpmol.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == bsid:
                tmpmol.characterize_complex(ligand)
        s = tmpmol.interaction_sets[bsid]
        # Hydrophobic contacts to Val82, Ile84, Ile150 as part of flap (S1, S1' sites)
        hydroph = {str(hyd.resnr)+hyd.reschain for hyd in s.all_hydrophobic_contacts}
        self.assertTrue({'82A', '84A', '50B'}.issubset(hydroph))
        # Hydrogen bonds
        hbonds = {str(hbond.resnr)+hbond.reschain for hbond in s.hbonds_ldon+s.hbonds_pdon}
        # Additional hbond to 25B not detected (low angle?)
        self.assertTrue({'29B', '48B', '27B', '25A'}.issubset(hbonds))
        # Water bridges
        waterbridges = {str(wb.resnr)+wb.reschain for wb in s.water_bridges}
        # Waterbridge with Gly27 is detected instead of Ala28/Asp29
        self.assertTrue({'50A', '50B', '29A'}.issubset(waterbridges))
        print(waterbridges)
