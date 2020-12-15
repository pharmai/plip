import unittest

from plip.basic import config
from plip.structure.preparation import PDBComplex, PLInteraction


def characterize_complex(pdb_file: str, binding_site_id: str) -> PLInteraction:
    pdb_complex = PDBComplex()
    pdb_complex.load_pdb(pdb_file)
    for ligand in pdb_complex.ligands:
        if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == binding_site_id:
            pdb_complex.characterize_complex(ligand)
    return pdb_complex.interaction_sets[binding_site_id]


class StructureProcessingTestCase(unittest.TestCase):
    def test_nmr(self):
        all_hydrogen_bonds = set()
        for i in range(1, 10):
            config.MODEL = i
            interactions = characterize_complex('./pdb/2ndo.pdb', 'SFQ:A:201')
            all_hbonds = interactions.hbonds_ldon + interactions.hbonds_pdon
            all_hydrogen_bonds.add(len(all_hbonds))
        # models contain from 0-2 hydrogen bonds
        self.assertEqual(all_hydrogen_bonds, {0, 1, 2})

    def test_nmr_invalid_model(self):
        config.MODEL = 11
        interactions = characterize_complex('./pdb/2ndo.pdb', 'SFQ:A:201')
        all_hbonds = interactions.hbonds_ldon + interactions.hbonds_pdon
        self.assertEqual(len(all_hbonds), 1)
