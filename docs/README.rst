===========
PLIProfiler
===========

The Protein-Ligand Interaction Profiler (PLIP) is a tool to analyze and visualize protein-ligand interactions in PDB files.

Requirements
============
lxml
pybel, openbabel
numpy
pymol

Features
========
* Direct download of PDB structures from RCSB PDB server if PDB ID is given
* Processing of custom PDB files (e.g. from docking)
* Automatic detection of bound ligands in a PDB file
* No need for special preparation of a PDB file, works out of the box
* Atom-level interaction reports in txt and XML formats for easy parsing
* Generation of [PyMOL] session files (.pse) for each pairing, showing a publication-ready illustration of interactions

Quickstart
==========
Run the PLIP python script with the option `-h` to show a quick start tutorial and list all available command line options.

Publication
===========
A proposal has been submitted to Nucleic Acids Research Web Server Issue 2015.

Documentation
=============


EXIT CODES
----------

1 : Unspecified Error
2 : Empty PDB file as input
3 : Invalid PDB ID
4 : PDB file can't be read by OpenBabel (due to invalid input files)
5 : PDB ID is valid, but wwPDB offers no file in PDB format for download.

Legend for PyMOL visualization
------------------------------

All colors given as RGB values.
<Description> - <RGB> - <PyMOL color> - <Representation>

Structural Elements
"""""""""""""""""""
Protein - [43, 131, 186] - myblue (custom) - sticks
Ligand - [253, 174, 97] - myorange (custom) - sticks
Water - [191, 191, 255] - lightblue - nb_spheres
Charge Center - [255, 255, 0] - yellow - spheres
Aromatic Ring Center - [230, 230, 230] -  grey90 - spheres

Interactions
""""""""""""
Hydrophobic Interaction - [128, 128, 128] - grey50 - dashed Line
Hydrogen Bond - [0, 0, 255] - blue - solid Line
Water Bridges - [191, 191, 255] - lightblue - solid Line
pi-Stacking (parallel) - [0, 255, 0] - green - dashed Line
pi-Stacking (perpendicular) - [140, 179, 102] - smudge - dashed Line
pi-Cation Interaction - [255, 128, 0] - orange - dashed Line
Halogen Bond - [64, 255, 191] - greencyan - solid Line
Salt Bridge - [255, 255, 0] - yellow - dashed Line

Code Contributions
------------------
Sebastian Salentin sebastian.salentin (at) biotec.tu-dresden.de
Joachim Haupt joachim.haupt (at) biotec.tu-dresden.de


## Contact Me

Questions or comments about `PLIProfiler`? Write me an email to sebastian.salentin (at) biotec.tu-dresden.de
