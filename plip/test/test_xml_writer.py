import unittest

from plip.exchange.report import StructureReport
from plip.structure.preparation import PDBComplex


class XMLWriterTest(unittest.TestCase):
    def test_pi_stacking(self):
        pdb_complex = PDBComplex()
        pdb_complex.load_pdb('./pdb/4dst_protonated.pdb')
        for ligand in pdb_complex.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == 'GCP:A:202':
                pdb_complex.characterize_complex(ligand)
                structure_report = StructureReport(pdb_complex, outputprefix="test_")
                structure_report.write_xml(as_string=True)

    def test_pication(self):
        pdb_complex = PDBComplex()
        pdb_complex.load_pdb('./pdb/6nhb.pdb')
        for ligand in pdb_complex.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == 'H4B:A:802':
                pdb_complex.characterize_complex(ligand)
                structure_report = StructureReport(pdb_complex, outputprefix="test_")
                structure_report.write_xml(as_string=True)