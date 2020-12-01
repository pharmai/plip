import os
import tempfile
import unittest

from plip.basic import config
from plip.basic.remote import VisualizerData
from plip.structure.preparation import PDBComplex
from plip.visualization.visualize import visualize_in_pymol


class VisualizationTest(unittest.TestCase):

    def setUp(self) -> None:
        self.tmp_dir = tempfile.mkdtemp()

    def test_visualization(self) -> None:

        pdb_file = './pdb/2ndo.pdb'
        binding_site_id = 'SFQ:A:201'
        config.PYMOL = True
        config.MODEL = 2
        config.OUTPATH = str(self.tmp_dir)
        pdb_complex = PDBComplex()
        pdb_complex.load_pdb(pdb_file)
        for ligand in pdb_complex.ligands:
            if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == binding_site_id:
                pdb_complex.characterize_complex(ligand)
        visualizer_complexes = [VisualizerData(pdb_complex, site) for site in sorted(pdb_complex.interaction_sets) if
                                not len(pdb_complex.interaction_sets[site].interacting_res) == 0]
        visualize_in_pymol(visualizer_complexes[0])
        self.assertEqual(1, len(os.listdir(self.tmp_dir)))
