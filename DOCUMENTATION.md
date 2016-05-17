# Protein-Ligand Interaction Profiler (PLIP)

## What is PLIP?
The Protein-Ligand Interaction Profiler (PLIP) is a tool to analyze and visualize protein-ligand interactions in PDB files.

## Why should I use PLIP?

#### Comprehensive Detection of Interaction
* Eight types of noncovalent interactions
* Interaction between proteins and
  * small molecules
  * ions
  * polymers
  * DNA / RNA
* Rich additiona information on binding, e.g. unpaired functional groups

#### Everything Is Automated

* Direct download of PDB structures from PDB server
* Automatic detection and grouping of relevant ligands in a PDB file
* No need for special preparation of a PDB file, works out of the box

#### Flexible Usage
* Processing of custom PDB files containing protein-ligand complexes (e.g. from docking)
* Atom-level interaction reports in TXT and XML formats for easy parsing
* Generation of PyMOL session files (.pse), enabling easy preparation of images for publications and talks
* Rendering of 3D interaction diagram for each ligand and its interactions with the protein

## Get PLIP running

#### 1. Install dependencies
The following instructions work on a unix machine.
If you are using Windows or OSX, the process may differ.
Please refer to the webpages of the corresponding tools to get help for installation.

To install all dependencies for PLIP, either just install PyMOL and then use the `pip` installation routine or install all tools and clone from GitHub.
The current version of PLIP depends on
* OpenBabel >=2.3.2
* PyMOL 1.7.x with Python bindings
* Imagemagick >=6.9.x (optional)

and should be executed with Python 2.7.x.

Example command for Ubuntu using apt-get
```bash
sudo apt-get install pymol openbabel python-openbabel imagemagick
```


#### 2. Get PLIP on your machine

There are several options to get PLIP on your machine.
If you have the Python package manager `pip` installed, you can install PLIP in a terminal via the following command:

```bash
pip install plip
```
If you want to clone from Github, open a new system terminal and clone the repository using
```bash
git clone https://github.com/ssalentin/plip.git ~/pliptool
```
You can use any other user directory, but we will use the one given above for the documentation.

#### 3. Simplify access to PLIP

The package manager `pip` will create symlinks for you so you can call PLIP using `plip` in the command line. If you cloned from Github, use the following command to do the same:

```bash
alias plip='~/pliptool/plip/plipcmd'
```

## Analyze a PDB structure with PLIP
Havin PLIP installed, you can run
```bash
plip -h
```
to list all available command line options.
We will go through all important settings in this section.
A typical application of PLIP involves the analysis of protein-ligand interactions in a structure from the Protein Data Bank (PDB).
PLIP can automatically fetch the entry from the PDB server when a valid PDB ID is provided.
```bash
plip -i 1vsn -v
```
The command above will fetch the PDB entry 1vsn from the server, analyze all interactions and print out the results (verbose mode).
No output files are produced at this point.

The same can be done for local PDB files.

```bash
wget http://files.rcsb.org/download/1EVE.pdb
plip -f 1EVE.pdb -v
```

The output formats can be added in any combination, currently including:
* XML report files (`-x`, best for automatic processing)
* Text report files (`-t`, human-readable)
* PyMOL session files (`-y`)
* Ray-traced images (`-p`)

```bash
plip -i 1osn -pyx
```

The command above will fetch the PDB entry 1osn, analye all interactions and produce one XML result file (`-x`) as well as rendered images (`-p`) and PyMOL session files (`-y`) for each binding site.

## Changing detection thresholds
The PLIP algorithm uses a rule-based detection to report non-covalent interaction between proteins and their partners.
The current settings are based on literature values and have been refined based on extensive testing with independent cases from mainly crystallography journals, covering a broad range of structure resolutions.
For most users, it is not recommended to change the standard settings.
However, there may be cases where changing detection thresholds is advisable (e.g. sets with multiple very low resolution structures)
PLIP allows you to change the setting permanently or for one run.

#### Permanent change
PLIP settings are stored in the file `modules/config.py`, which is loaded in each run.

#### Temporary change

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

__Structural Elements__


Protein - [43, 131, 186] - myblue (custom) - sticks

Ligand - [253, 174, 97] - myorange (custom) - sticks

Water - [191, 191, 255] - lightblue - nb_spheres

Charge Center - [255, 255, 0] - yellow - spheres

Aromatic Ring Center - [230, 230, 230] -  grey90 - spheres

Ions - [255, 255, 128] - hotpink - spheres


__Interactions__


Hydrophobic Interaction - [128, 128, 128] - grey50 - dashed Line

Hydrogen Bond - [0, 0, 255] - blue - solid Line

Water Bridges - [191, 191, 255] - lightblue - solid Line

pi-Stacking (parallel) - [0, 255, 0] - green - dashed Line

pi-Stacking (perpendicular) - [140, 179, 102] - smudge - dashed Line

pi-Cation Interaction - [255, 128, 0] - orange - dashed Line

Halogen Bond - [64, 255, 191] - greencyan - solid Line

Salt Bridge - [255, 255, 0] - yellow - dashed Line

Metal Complexation - [140, 64, 153] - violetpurple - dashed Line


XML Report Documentation
------------------------


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

## Testing

TODO

## Troubleshooting

TODO
