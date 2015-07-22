==================
PLIProfiler v1.2.0
==================

The Protein-Ligand Interaction Profiler (PLIP) is a tool to analyze and visualize protein-ligand interactions in PDB files.


Features
========
* Detection of eight different types of noncovalent interactions, including metal complexes
* Works for complexes of protein with small molecules, ions, polymers, or DNA/RNA (and all combinations)
* Automatic detection and grouping of relevant ligands in a PDB file
* Rich additional information on binding, e.g. unpaired functional groups
* Direct download of PDB structures from wwPDB server if valid PDB ID is given
* Processing of custom PDB files containing protein-ligand complexes (e.g. from docking)
* No need for special preparation of a PDB file, works out of the box
* Atom-level interaction reports in TXT and XML formats for easy parsing
* Generation of PyMOL session files (.pse), enabling easy preparation of images for publications and talks
* Rendering of 3D interaction diagram for each ligand and its interactions with the protein

Quickstart
==========
Run the PLIP python script with the option `-h` to list all available command line options.
To analyze a protein-ligand complex from a Protein Data Bank entry -- e.g. 1vsn --, run
    `plipcmd -i 1vsn`.
To analyze a PDB file from your workstation, run
    `plipcmd -f path_to_pdbfile.pdb`.
The output format(s) can be chosen freely, ranging from ...
* XML report files (`-x`, highest level of detail)
* Text report files (`-t`, medium level of detail)
* PyMOL session files (`-y`)
* Ray-traced images (`-p`)
* Verbose output on command line (`-v`)

Threshold settings
==================
All geometric thresholds used for detection of interactions in PLIP are stored in the `config.py` module and can be
changed there permanently if desired. Another possibility is to change single parameters via command line options when
running plip. The naming of the variables is identical to those in the config file, but all lowercase. Specify the
threshold you want to change and the new value, e.g. to change HYDROPH_DIST_MAX to 5 Angstrom, run PLIP using
    `python plip-cmd.py -i 1vsn --hydroph_dist_max 5`
All distance thresholds can be increased to up to 10 Angstrom. Thresholds for angles can be set between 0 and 180 degree.
If two interdependent thresholds have conflicting values, PLIP will show an error message.

Web Service
===========
A web service for analysis of protein-ligand complexes using PLIP is available at
http://projects.biotec.tu-dresden.de/plip-web
The web site offers advanced functions to search for specific entries from PDB and lists the interaction results in the browser.


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
3 : Invalid PDB ID (wrong format)
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
Ions - [255, 255, 128] - hotpink - spheres

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
Metal Complexation - [140, 64, 153] - violetpurple - dashed Line

Output Files
------------
All output files contain information on non-covalent interactions between the protein and all relevant ligands in the PDB structure.

XML/RST Result Files
""""""""""""""""""""

**report**
Contains all binding site information

**plipversion**
Version of PLIP used for generating the output file

**bindingsite**
Information for one bindingsite. Has a unique ID and attribute `has_interactions`

**identifiers**
Ligand/bindingsite identifiers

**longname**
Long name of ligand, contains all het ids of ligands in one composite cluster

**ligtype**
Classification of ligand, can be SMALLMOLECULE/POLYMER/DNA/RNA/ION or combinations of the first four with ION

**hetid**
PDB hetero ID of the ligand

**chain**
Chain assigned to the ligand in the PDB file

**position**
Position in chain of the ligand in the PDB file

**composite**
Can be True or False depending on whether the ligand consists of several separate subunits or not

**members**
Lists the members of a composite ligand cluster

**smiles**
The SMILES string of the complete (composite) ligand

**lig_properties**
Additional information on the ligand, i.e. number of functional atoms

**num_heavy_atoms**
Number of heavy atoms in the ligand

**num_hbd**
Number of hydrogen bond donors in the ligand

**num_unpaired_hbd**
Number of unpaired hydrogen bond donors in the ligand (not involved as acceptor/donor in hydrogen bonds, salt bridges,
water bridges, metal complexes)

**num_hba**
Number of hydrogen bond acceptors in the ligand

**num_unpaired_hba**
Number of unpaired hydrogen bond acceptors in the ligand (not involved as acceptor/donor in hydrogen bonds, salt bridges,
water bridges, metal complexes)

**num_hal**
Number of halogen bond donors in the ligand

**num_unpaired_hal**
Number of unpaired halogen bond donors in the ligand

**num_aromatic_rings**
Number of aromatic rings in the ligand

**interacting chains**
Lists the chains the ligand interacts with

**bs_residues**
Listing of binding site residues the ligand is near to or interacts with. Contains the type of amino acid, information
on contact, a unique id and the minimal distance to the ligand in Angstrom

**interactions**
Detailed information on all interactions (general attributes documented below)

**resnr**
Residue number of interacting amino acid

**restype**
Residue type of interacting amino acid

**reschain**
Residue chain of interacting amino acid

**dist**
Distance of interacting atoms or groups in Angstrom

**ligcoo**
Coordinates of interacting ligand atom or interaction center in ligand

**protcoo**
Coordinates of interacting protein atom or interaction center in ligand




