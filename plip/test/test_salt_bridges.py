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

class SaltBridgeTest(unittest.TestCase):

    def test_4yb0(self):
        '''test salt bridge detection for nucleic acids as part of the receptor'''

        config.DNARECEPTOR = True
        interactions = characterize_complex('./pdb/4yb0.pdb', 'C2E:R:102')
        salt_bridges = interactions.saltbridge_lneg + interactions.saltbridge_pneg
        self.assertTrue(len(salt_bridges) == 1)