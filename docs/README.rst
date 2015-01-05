==================
PLIProfiler v1.0.0
==================

The Protein-Ligand Interaction Profiler (PLIP) is a tool to analyze and visualize protein-ligand interactions in PDB files.

Requirements
============
lxml
openbabel
pybel
numpy
pymol
imagemagick (optional)

Features
========
* Detection of eight different types of noncovalent interactions
* Automatic detection of relevant ligands in a PDB file
* Direct download of PDB structures from wwPDB server if valid PDB ID is given
* Processing of custom PDB files containing protein-ligand complexes (e.g. from docking)
* No need for special preparation of a PDB file, works out of the box
* Atom-level interaction reports in rST and XML formats for easy parsing
* Generation of PyMOL session files (.pse) for each pairing, enabling easy preparation of images for publications and talks
* Rendering of preview image for each ligand and its interactions with the protein

Quickstart
==========
Run the PLIP python script with the option `-h` to list all available command line options.
To analyze a protein-ligand complex from a Protein Data Bank entry -- e.g. 1vsn --, run `python plip.cmd.py -i 1vsn`.
To analyze a PDB file from your workstation, run `python plip.cmd.py -f path_to_pdbfile.pdb`.

Web Service
===========
A web service for analysis of protein-ligand complexes using PLIP is available at
http://projects.biotec.tu-dresden.de/plip-web
The web site offers advanced functions to search for specific entries from PDB and lists the interaction results in the browser.


Publication
===========
A proposal has been submitted to Nucleic Acids Research Web Server Issue 2015.

Documentation
=============

Algorithm
---------
PLIP uses a rule-based algorithm for the detection of non-covalent interactions. For details on the algorithm, visit
http://projects.biotec.tu-dresden.de/plip-web

Exit codes
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

Output Files
------------
All output files contain information on non-covalent interactions between the protein and all relevant ligands in the PDB structure.

XML Result Files
""""""""""""""""

<pdbid>
Unique identifier for the corresponding entry of the protein structure in Protein Data Bank.

<hetid>
Unique identifier for ligand molecule in Protein Data Bank (PDB).

<chain>
One protein can consist of multiple separate amino acids chains which are named alphabetically

<position>
Position of ligand in PDB numbering. Combination of pdbid, hetid, chain and position gives a unique identifier for
each protein-ligand complex. Same numbering as <resnr>

<interactions>
Contains interaction for protein-ligand complex, organized by interaction type, e.g. hydrophobic interactions

<resnr>
Position of amino acid in protein chain according to PDB numbering

<restype>
Amino acid type in three-letter code

<dist*>
Distance of interacting atoms

<*idx>
Atom ID in original PDB structure

<lig_idx_list>
Atom IDs if several ligand atoms are relevant for a single interaction (e.g. when forming a charge center)

<*angle>
Angle between interacting groups

<protispos>, <protisdon>, <protischarged>
Determines if the protein is positively charged, provides a donor or a charge.
Important for interactions with directionality.

<sidechain>
Is true if a hydrogen bond is formed with the sidechain of the protein and false if it is formed with the backbone.

<ligcoo>, <protcoo>
Coordinates of protein and ligand interacting atoms or interaction centers (e.g. charge centers)


Contact Me
----------

Questions or comments about `PLIProfiler`? Write me an email to sebastian.salentin (at) biotec.tu-dresden.de
