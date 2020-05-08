import unittest

from plip.structure.preparation import PDBComplex, PLInteraction


def characterize_complex(pdb_file: str, binding_site_id: str) -> PLInteraction:
    pdb_complex = PDBComplex()
    pdb_complex.load_pdb(pdb_file)
    for ligand in pdb_complex.ligands:
        if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == binding_site_id:
            pdb_complex.characterize_complex(ligand)
    return pdb_complex.interaction_sets[binding_site_id]


class WaterBridgeTest(unittest.TestCase):

    def test_3ems(self):
        interactions = characterize_complex('./pdb/3ems.pdb', 'ARG:A:131')
        water_bridges = interactions.water_bridges
        self.assertEqual(len(water_bridges), 4)
