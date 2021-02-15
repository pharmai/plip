import unittest

from plip.basic import config
from plip.exchange.report import StructureReport
from plip.structure.preparation import PDBComplex


class IntraTest(unittest.TestCase):

    def test_4day(self):
        config.PEPTIDES = ['C']
        pdb_complex = PDBComplex()
        pdb_complex.load_pdb('./pdb/4day.pdb')
        for ligand in pdb_complex.ligands:
            pdb_complex.characterize_complex(ligand)
            structure_report = StructureReport(pdb_complex, outputprefix="test_")
            structure_report.write_xml(as_string=True)
        config.PEPTIDES = []
