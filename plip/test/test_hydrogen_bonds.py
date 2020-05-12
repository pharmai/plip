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


class HydrogenBondTestCase(unittest.TestCase):

    def test_4dst_nondeterministic_protonation(self):
        config.NOHYDRO = False
        for i in range(0, 10):
            interactions = characterize_complex('./pdb/4dst.pdb', 'GCP:A:202')
            all_hbonds = interactions.hbonds_ldon + interactions.hbonds_pdon
            self.assertTrue(len(all_hbonds) == 16 or len(all_hbonds) == 17)

    def test_4dst_deterministic_protonation(self):
        config.NOHYDRO = True
        for i in range(0, 10):
            interactions = characterize_complex('./pdb/4dst_protonated.pdb', 'GCP:A:202')
            all_hbonds = interactions.hbonds_ldon + interactions.hbonds_pdon
            self.assertTrue(len(all_hbonds) == 16)

    def test_no_protonation(self):
        config.NOHYDRO = True
        interactions1 = characterize_complex('./pdb/1x0n_state_1.pdb', 'DTF:A:174')
        self.assertEqual(len(interactions1.hbonds_ldon), 0)
        config.NOHYDRO = False
        interactions2 = characterize_complex('./pdb/1x0n_state_1.pdb', 'DTF:A:174')
        self.assertEqual(len(interactions2.hbonds_ldon), 1)
