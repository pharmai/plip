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


class RingDetectionTest(unittest.TestCase):

    def test_consistent_ring_detection(self):
        config.NOHYDRO = True
        angles = set()
        for i in range(0, 10):
            interactions = characterize_complex('./pdb/4dst_protonated.pdb', 'GCP:A:202')
            angles.add(interactions.pistacking[0].angle)
        self.assertTrue(len(angles) == 1)
        config.NOHYDRO = False
