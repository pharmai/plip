"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
test_command_line.py - Unit Tests for special cases.
"""
import os
import subprocess
import sys
import tempfile
import unittest


class CommandLineTest(unittest.TestCase):
    """Checks special and extreme cases for input files."""

    def setUp(self) -> None:
        self.tmp_dir = tempfile.TemporaryDirectory()

    def tearDown(self) -> None:
        self.tmp_dir.cleanup()

    def test_empty_input_file(self):
        """Input file is empty."""
        exitcode = subprocess.call(f'{sys.executable} ../plipcmd.py -f ./special/empty.pdb -o {self.tmp_dir.name}',
                                   shell=True)
        self.assertEqual(exitcode, 1)

    def test_invalid_pdb_id(self):
        """A PDB ID with no valid PDB record is provided."""
        exitcode = subprocess.call(f'{sys.executable} ../plipcmd.py -i xx1x -o {self.tmp_dir.name}', shell=True)
        self.assertEqual(exitcode, 1)

    def test_invalid_input_file(self):
        """A file is provided which is not a PDB file."""
        exitcode = subprocess.call(f'{sys.executable} ../plipcmd.py -f ./special/non-pdb.pdb -o {self.tmp_dir.name}',
                                   shell=True)
        self.assertEqual(exitcode, 1)

    def test_pdb_format_not_available(self):
        """A valid PDB ID is provided, but there is no entry in PDB format from wwPDB"""
        exitcode = subprocess.call(f'{sys.executable} ../plipcmd.py -i 4v59 -o {self.tmp_dir.name}', shell=True)
        self.assertEqual(exitcode, 1)

    def test_valid_pdb(self):
        """A PDB ID with no valid PDB record is provided."""
        exitcode = subprocess.call(f'{sys.executable} ../plipcmd.py -x -f ./pdb/1eve.pdb -o {self.tmp_dir.name}',
                                   shell=True)
        self.assertEqual(len(os.listdir(self.tmp_dir.name)), 2)
        self.assertEqual(exitcode, 0)
