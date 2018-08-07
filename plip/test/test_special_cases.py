# coding=utf-8
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
test_special_cases.py - Unit Tests for special cases.
"""


import unittest
import subprocess


class SpecialCasesTest(unittest.TestCase):
    """Checks special and extreme cases for input files."""

    def test_empty_input_file(self):
        """Input file is empty."""
        exitcode = subprocess.call('python ../plipcmd.py -f ./special/empty.pdb -o /tmp', shell=True)
        self.assertEqual(exitcode, 2)  # Specific exitcode 2

    def test_invalid_pdb_id(self):
        """A PDB ID with no valid PDB record is provided."""
        exitcode = subprocess.call('python ../plipcmd.py -i xx1x -o /tmp', shell=True)
        self.assertEqual(exitcode, 3)  # Specific exitcode 3

    def test_invalid_input_file(self):
        """A file is provided which is not a PDB file."""
        exitcode = subprocess.call('python ../plipcmd.py -f ./special/non-pdb.pdb -o /tmp', shell=True)
        self.assertEqual(exitcode, 4)  # Specific exitcode 4

    def test_pdb_format_not_available(self):
        """A valid PDB ID is provided, but there is no entry in PDB format from wwPDB"""
        exitcode1 = subprocess.call('python ../plipcmd.py -i 4v59 -o /tmp', shell=True)
        self.assertEqual(exitcode1, 5)  # Specific exitcode 5
