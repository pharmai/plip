import unittest

from plip.modules.preparation import PDBComplex, PLInteraction


def characterize_complex(pdb_file: str, binding_site_id: str) -> PLInteraction:
    pdb_complex = PDBComplex()
    pdb_complex.load_pdb(pdb_file)
    for ligand in pdb_complex.ligands:
        if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == binding_site_id:
            pdb_complex.characterize_complex(ligand)
    return pdb_complex.interaction_sets[binding_site_id]


class HydrogenBondTestCase(unittest.TestCase):

    def test_4dst(self):
        interactions = characterize_complex('./pdb/4dst.pdb', 'GCP:A:202')
        all_hbonds = interactions.hbonds_ldon + interactions.hbonds_pdon
        self.assertTrue(len(all_hbonds) == 16 or len(all_hbonds) == 17)
