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
* Automatic fixiing of errors in PDB files

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
* PyMOL 1.7.x with Python bindings (optional, for visualization only)
* Imagemagick >=6.9.x (optional)

and should be executed with Python 2.7.x.

Example command for Ubuntu using apt-get
```bash
sudo apt-get install pymol openbabel python-openbabel imagemagick swig
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

A third option is to use the [PLIP Conda package](https://anaconda.org/plip/plip).

#### 3. Simplify access to PLIP

The package manager `pip` will create symlinks for you so you can call PLIP using `plip` in the command line. If you cloned from Github, use the following command to do the same:

```bash
alias plip='python ~/pliptool/plip/plipcmd.py'
```

## Analyze a PDB structure with PLIP

### Single Structures
Having PLIP installed, you can run
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

The same can be done for local PDB files (-f <file>) or for reading from stdin (-f -).

```bash
wget http://files.rcsb.org/download/1EVE.pdb
plip -f 1EVE.pdb -v
```

The output formats can be added in any combination, currently including:
* XML report files (`-x`, best for automatic processing)
* Text report files (`-t`, human-readable)
* PyMOL session files (`-y`)
* Ray-traced images (`-p`)
* writing to stdout (`-O`), to be used in combination with XML or text report files

```bash
plip -i 1osn -pyx
```

The command above will fetch the PDB entry 1osn, analyze all interactions and produce one XML result file (`-x`) as well as rendered images (`-p`) and PyMOL session files (`-y`) for each binding site.

### Batch Mode
PLIP can process multiple structures at once, either from PDB or local files.
To activate batch mode, just provide a list of PDB IDs or local file names, e.g.:

```bash
plip -i 1vsn 1osn 2reg -vx
```

PLIP will create subdirectories for each given structure in the output folder.
If in PDB ID mode (`i`), the folder structure will be nested and based on the two middle characters of the PDB ID.
The structure 1vsn in batch processing will have its output files in <outputfolder>/vs/1vsn .

### Detection of protein-peptide interactions
For the detection of ligands, PLIP relies on the separation of ATOM and HETATM entries in the PDB file.
The latter are searched for suitable ligands when running in normal mode.
Peptide ligands, however, are usually deposited as ATOM entries in a separate chain.
PLIP can not detect these entities automatically.
To switch into protein-peptide interaction mode, start PLIP with the option `--peptide`, followed by the peptide chain of interest, e.g.:

```bash
plip -i 5hi4 --peptides I -vx
```

### Detection of intra-chain interactions (__beta__)
Intra-protein interactions are important for the stabilization of a structure and can give valuable insights for protein engineering and drug discovery.
PLIP now supports detection of interactions within one chain.
To switch into intra-chain interaction mode, start PLIP with the option `--intra`, followed by the protein chain of interest, e.g.:

```bash
plip -i 5b2m --intra A -yv
```
Please note that detection within a chain takes much longer than detection of protein-ligand interactions,
especially for large structures.

### Interactions of molecules with DNA/RNA (__beta__)
PLIP can characterize interactions between ligands and DNA/RNA.
A special mode allows to switch from protein to DNA/RNA as the receptor molecule in the structure.
To change the receptor mode, start PLIP with the option `--dnareceptor`.


## Changing detection thresholds
The PLIP algorithm uses a rule-based detection to report non-covalent interaction between proteins and their partners.
The current settings are based on literature values and have been refined based on extensive testing with independent cases from mainly crystallography journals, covering a broad range of structure resolutions.
For most users, it is not recommended to change the standard settings.
However, there may be cases where changing detection thresholds is advisable (e.g. sets with multiple very low resolution structures)
PLIP allows you to change the setting permanently or for one run.

#### Permanent change
PLIP settings are stored in the file `modules/config.py`, which is loaded in each run.
One section in the file includes assignments of standard values to all distance and angle thresholds, together with short descriptions.
Changing the values in this file means they will be changed permanently for PLIP.
Furthermore, filter settings for ligands (metal ions in complexes, unsupported ligands, warnings for possible artifacts) can be changed in this file as well.

#### Temporary change

Thresholds can be changed for each run using command-line parameters.
The naming of the variables is identical to those in the config file, but all lowercase.
Specify the
threshold you want to change and the new value, e.g. to change HYDROPH_DIST_MAX to 5 Angstrom, run PLIP using
```bash
plip -i 1vsn --hydroph_dist_max 5
```
All distance thresholds can be increased to up to 10 Angstrom. Thresholds for angles can be set between 0 and 180 degree.
If two interdependent thresholds have conflicting values, PLIP will show an error message.

## Further Options
PLIP offers further command line options which enables you to switch advanced settings.
* Set number `n` of maximum threads used for parallel processing (`--maxthreads <n>`)
* Do not automatically combine covalently bound ligands (`--breakcomposite`)
* Do not discard alternate locations (`--altlocation`)
* Set debug mode (`--debug`)
* Turn off automatic fixing of errors in PDB files (`--nofix`)
* Keep modified residues as ligands (`--keepmod`)

## Web Service
A web service for analysis of protein-ligand complexes using PLIP is available at
http://plip.biotec.tu-dresden.de/
The web site offers advanced functions to search for specific entries from PDB and lists the interaction results in the browser.
Additionally, the service used the BioLiP database to annotate biologically relevant ligands.
The option to change threshold, ligand filtering, and batch processing is only available in the command line tool and with the Python modules.

## Algorithm
PLIP uses a rule-based system for detection of non-covalent interactions between protein residues and ligands.
Information on chemical groups able to participate in a specific interaction (e.g. requirements for hydrogen bond donors) and interaction geometry (e.g. distance and angle thresholds) from literature are used to detect characteristics of non-covalent interactions between contacting atoms of protein and ligands.
For each binding site, the algorithm searches first for atoms or atom groups in the protein and ligand which could possibly be partner in specific interactions.
In the second step, geometric rules are applied to match groups in protein and ligand forming an interaction.

#### Detection and filtering of ligands
Previous to the detection step for the interactions, PLIP extracts all ligands contained in the structure.
Modified amino acids are identified and excluded using MODRES entries of the PDB files.
Additionally, it uses the BioLiP list of possible artifacts to remove ligands which are in this list and appear 15 times or more in a structure.
Just a few compounds are currently excluded, being listed in the PLIP config file.

#### Preparation of structures
Polar hydrogens are added to the structure and alternative conformations/models/positions removed.
Missing chains are assigned to ligands and non-standard ligand names (with special characters) altered to LIG.

#### Detection of possible interacting groups

##### Binding site atoms
The binding site distance cutoff is determined by adding up BS_DIST_MAX to the maximum extent to the ligand (maximum distance of a ligand atom to ligand centroid).
All protein atoms within this distance cutoff to any binding site atoms are counted as belonging to the binding site.

##### Hydrophobic Atoms
An atom is classified as hydrophobic if it is a carbon and has only carbon or hydrogen atoms as neighbours.

##### Aromatic Rings
OpenBabel is used to identify rings (SSSR perception) and their aromaticity.
In cases where no aromaticity is reported by OpenBabel, the ring is checked for planarity.
To this end, the normals of each atom in the ring to its neighbors is calculated.
The angle between each pair of normals has to be less than AROMATIC_PLANARITY.
If this holds true, the ring is also considered as aromatic.

##### Hydrogen Bond Donors and Acceptors
OpenBabel is used to identify hydrogen bond donor and acceptor atoms.
Halogen atoms are excluded from this group and treated separately (see below).

##### Charged Groups
The detection of charged groups is only exhaustive for the binding site, not the ligands.
For proteins, positive charges are attributed to the side chain nitrogens of Arginine, Histidine and Lysine.
Negative charged are assigned to the carboxyl groups in Aspartic Acid and Glutamic Acid.
In ligands, positive charges are assigned to quaterny ammonium groups, tertiary amines (assuming the nitrogen could pick up a hydrogen and thus get charged), sulfonium and guanidine groups.
Negative charges are reported for phosphate, sulfonate, sulfonic acid and carboxylate.

##### Halogen bonds donors and acceptors
Assuming that halogen atoms are not present in proteins (unless they are artificially modified), halogen bond donors are searched for only in ligands.
All fluorine, chlorine, bromide or iodine atoms connected to a carbon atom qualify as donors.
Halogen bond acceptors in proteins are all carbon, phosphor or sulphur atoms connected to oxygen, phosphor, nitrogen or sulfur.

##### Water
Water atoms are assigned to a ligand-binding site complex if their oxygen atoms are within a certain cutoff to the ligand.
The cutoff is determined by adding up BS_DIST_MAX to the maximum extent to the ligand (maximum distance of a ligand atom to ligand centroid).
This means the farthest distance of a ligand to a water atom is BS_DIST_MAX.

#### Detection of interactions
For an overview on geometric cutoffs used for the prediction of interactions, please refer to the config file.
Note that the threshold can not be changed for jobs running on PLIP Web.
The command line tool (sourcecode available for download), however, allows changing all listed parameters permanently or for single runs.

##### Hydrophobic Interactions
As hydrophobic interactions result from entropic changes rather than attractive forces between atoms, there are no clear geometries of hydrophobic association.
The observed attraction between hydrophobic atoms decays exponentionally with the distance between them.

A generous cutoff was chosen, identifying a prime set of hydrophobic interactions between all pairs of hydrophobic atoms within a distance of HYDROPH_DIST_MAX.
Since the number of hydrophobic interactions with such an one-step approach can easily surpass all other interaction types combined, it may strongly influence subsequent evaluation or applications as interaction fingerprinting.
To overcome this problem, the number of hydrophobic interactions is reduced in several steps.
First, hydrophobic interactions between rings interacting via π-stacking are removed.
This is done because stacking already involves hydrophobic interactions.

Second, two clustering steps are applied.
If a ligand atom interacts with several binding site atoms in the same residue, only the interaction with the closest distance is kept.
Subsequently, the set of hydrophobic interactions is checked from the opposite perspective: if a protein atom interacts with several neighboring ligand atoms, just the interaction with the closest distance is kept.
Together, these reduction steps help to report only the most representative hydrophobic interactions.

##### Hydrogen Bonds
A hydrogen bond between a hydrogen bond donor and acceptor is reported if several geometric requirements are fulfilled.
The distance has to be less than HBOND_DIST_MAX and the angle at the donor group (D-H...A) above HBOND_DON_ANGLE_MIN.

Since salt bridges involve purely electrostatic interactions as well as hydrogen bonds, it is not meaningful to report both interaction types between the same groups.
Thus, hydrogen bonds between atoms are removed if they belong to groups that already form a salt bridge to that atom.

As a general rule, a hydrogen bond donor can take part in only one hydrogen bond, while acceptor atoms can be partners in multiple hydrogen bonds (e.g. bifurcated hydrogen bonds).
For multiple possible hydrogen bonds from one donor, only the contact with the donor angle closer to 180 ° is kept.

##### π-Stacking
π-Stacking for two aromatic rings is reported whenever their centers are within a distance of PISTACK_DIST_MAX, the angle deviates no more than PISTACK_ANG_DEV from the optimal angle of 90 ° for T-stacking or 180 ° for P-stacking.

Additionally, each ring center is projected onto the opposite ring plane.
The distance between the other ring center and the projected point (i.e. the offset) has to be less than PISTACK_OFFSET_MAX.
This value corresponds approximately to the radius of benzene + 0.6 Å.

##### π-Cation Interactions
π-Cation interactions are reported for each pairing of a positive charge and an aromatic ring if the distance between the charge center and the aromatic ring center is less than PICATION_DIST_MAX and the offset between the ring center and the charge is no larger than PISTACK_OFFSET_MAX.
In the case of a putative π-cation interaction with a tertiary amine of the ligand, an additional angle criterion is applied (see documentation in the source code).

##### Salt Bridges
Whenever two centers of opposite charges come within a distance of SALTBRIDGE_DIST_MAX, a salt bridge is reported.
In contrast to hydrogen bonds, there are no additonal geometric restrictions.

##### Water bridges
While residues can be bridged by more than one water molecule, for the prediction in this script the only case considered is one water molecule bridging ligand and protein atoms via hydrogen bonding.

The water molecule has to be positioned between hydrogen bond donor/acceptor pairs of ligand and protein with distances of the water oxygen within WATER_BRIDGE_MINDIST and WATER_BRIDGE_MAXDIST to the corresponding polar atoms of the donor or acceptor groups.
If a constellation with a water atom fulfils these requirements, two angles are checked.
The angle ω between the acceptor atom, the water oxygen and donor hydrogen has to be within WATER_BRIDGE_OMEGA_MIN and WATER_BRIDGE_OMEGA_MAX.
Additionally, the angle θ between the water oxygen, the donor hydrogen and the donor atom has to be larger than WATER_BRIDGE_THETA_MIN.

Similar to standard hydrogen bonds, a water molecule is only allowed to participate as donor in two hydrogen bonds (two hydrogen atoms as donors).
In the case of more than two possible hydrogen bonds for a water molecule as donor, only the two contacts with a water angle closest to 110 ° are kept

##### Halogen Bonds
Halogen bonds are reported for each pairing of halogen bond acceptor and donor group having a distance of less than HALOGEN_DIST_MAX and angles at the donor and acceptor group of HALOGEN_DON_ANGLE and HALOGEN_ACC_ANGLE with a deviation of no more than HALOGEN_ANG_DEV

##### Metal Complexes
For metal complexes, PLIP considers metal ions from a set of more than 50 species (see PLIP config for more details).
Possible interacting groups in the protein are sidechains of cystein (S), histidine (N), asparagine, glutamic acid, serin, threonin, and tyrosin (all O), as well as all main chain oxygens.

In ligands, following groups are considered for metal complexation: alcohols, phenolates, carboxylates, phosphoryls, thiolates, imidazoles, pyrroles, and the iron-sulfur cluster as a special constellation.
For one metal ions, all groups with a maximum distance of METAL_DIST_MAX to the ligand are considered for the complex.

After assigning all target groups to one metal ions, the resulting set of angles of the complex is compared with known sets of angles from common coordination geometries (linear [2], trigonal planar [3], trigonal pyramidal [3], tetrahedral [4], square planar [4], trigonal bipyramidal [5], square pyramidal [5], and octahedral [6]).
The best fit with the least difference in observed targets is chosen as an estimated geometry and targets superfluous to the constellation are removed.

## Additional Information

### Exit error codes

| Exit code  | Description |
| ------------- | ------------- |
| **1**  | Unspecified Error  |
| **2**  | Empty PDB file as input |
| **3**  | Invalid PDB ID (wrong format)  |
| **4**  | PDB file can't be read by OpenBabel  |
| **5**  | Valid PDB ID, but no PDB file on wwwPDB  |

### Legend for PyMOL visualization

#### Structural Elements

| Description  | RGB | PyMOL color | Representation |
| ------------ | --- | ------------| ---------------|
| Protein  | [43, 131, 186] | myblue (custom) | sticks |
| Ligand  | [253, 174, 97] | myorange (custom) | sticks |
| Water  | [191, 191, 255]  | lightblue | nb_spheres |
| Charge Center | [255, 255, 0] | yellow | spheres |
| Aromatic Ring Center  | [230, 230, 230] | grey90 | spheres |
| Ions | [250, 255, 128] | hotpink | spheres |


#### Interactions

| Description  | RGB | PyMOL color | Representation |
| ------------ | --- | ------------| -------------- |
| Hydrophobic Interaction  | [128, 128, 128] | grey50 | dashed |
| Hydrogen Bond  | [0, 0, 255] | blue | solid line |
| Water Bridges  | [191, 191, 255]  | lightblue | solid line |
| pi-Stacking (parallel) | [0, 255, 0] | green | dashed line |
| pi-Stacking (perpendicular)  | [140, 179, 102] | smudge | dashed line |
| pi-Cation Interaction | [255, 128, 0] | orange | dashed line |
| Halogen Bond | [54, 255, 191] | greencyan | solid line |
| Salt Bridge | [255, 255, 0] | yellow | dashed line |
| Metal Complex | [140, 64, 153] | violetpurple | dashed line |


### XML Report Documentation

| Attribute  | Description |
| ------------ | -------------- |
| report | Contains all binding site information |

#### report

| Attribute  | Description |
| ------------ | -------------- |
| plipversion | Version of PLIP used for generating the output file |
| bindingsite | Information for one bindingsite. Has a unique ID and attribute `has_interactions` |
| date_of_creation | Date of the analyis |
| citation_information | How to cite PLIP |
| mode | Documents if PLIP was started in default or any special mode (e.g. intra-protein interactions) |
| pdbid | PDB identifier of the input file |
| pdbfile | Name of the input PDB file |
| pdbfixes | Were any fixes applied automatically to the input file? |
| filename | Filename of the processed PDB file |
| excluded_ligands | List of excluded ligands |

### bindingsite

| Attribute  | Description |
| ------------ | -------------- |
| identifiers | Ligand/bindingsite identifiers |
| lig_properties | Additional information on the ligand, i.e. number of functional atoms |
| interacting chains | Lists the chains the ligand interacts with |
| bs_residues | Listing of binding site residues the ligand is near to or interacts with.  |
| interactions | Detailed information on all interactions |
| mappings | Contains mappings from canonical SMILES to PDB in smiles_to_pdb |

### identifiers

| Attribute  | Description |
| ------------ | -------------- |
| longname | Long name of ligand, contains all het ids of ligands in one composite cluster |
| ligtype | Classification of ligand, can be SMALLMOLECULE/POLYMER/DNA/RNA/ION or combinations of the first four with ION |
| hetid | PDB hetero ID of the ligand |
| chain | Chain assigned to the ligand in the PDB file |
| position | Position in chain of the ligand in the PDB file |
| composite | Can be True or False depending on whether the ligand consists of several separate subunits or not |
| members | Lists the members of a composite ligand cluster |
| smiles | The SMILES string of the complete (composite) ligand |

### lig_properties

| Attribute  | Description |
| ------------ | -------------- |
| num_heavy_atoms | Number of heavy atoms in the ligand |
| num_hbd | Number of hydrogen bond donors in the ligand |
| num_unpaired_hbd | Number of unpaired hydrogen bond donors in the ligand |
| num_hba | Number of hydrogen bond acceptors in the ligand |
| num_unpaired_hba | Number of unpaired hydrogen bond acceptors in the ligand |
| num_hal | Number of halogen bond donors in the ligand |
| num_unpaired_hal | Number of unpaired halogen bond donors in the ligand |
| num_aromatic_rings | Number of aromatic rings in the ligand |
| num_rotatable_bonds | Number of rotatable bonds in the ligand |
| molweight | Molecular weight of the ligand |
| logp | logP value of the ligand |

### bs_residues

Contains a list of bs_residue entries selected with a primary distance cutoff. As IDs, ech bs_residue lists aa (three-letter amino acid code), contact (boolean, indicated if residue interactions with the ligand or not), id (a continuous ID), min_dist (The minimal distance to the ligand).

| Attribute  | Description |
| ------------ | -------------- |
bs_residue | Concatenated chain of position and residue number |

### interactions

Interactions are subdivided by interaction type. Most interaction types share the same information (types and IDs of interacting residue, distances), while some information is specific to certain types (charges, directionality, ring geometries, etc.)

Type | Attribute  | Description |
| ------------| ------------ | -------------- |
All | resnr | Residue number of interacting amino acid |
All | restype | Residue type of interacting amino acid |
All | reschain | Residue chain of interacting amino acid |
All | resnr_lig | Residue number of interacting ligand residue |
All | restype_lig | Residue type of interacting ligand residue |
All | reschain_lig | Residue chain of interacting ligand residue
All | dist | Distance of interacting atoms or groups in Angstrom |
All (except metal_complex) | ligcoo | Coordinates of interacting ligand atom or interaction center in ligand |
All (except metal_complex)| protcoo | Coordinates of interacting protein atom or interaction center in ligand |
hydrogen_bond | sidechain | Is the H-Bond formed with the sidechain of the protein? |
hydrogen_bond | dist_h-a | Distance between H-Bond hydrogen and acceptor atom |
hydrogen_bond | dist_d-a | Distance between H-Bond donor and acceptor atoms |
hydrogen_bond, water_bridge, halogen_bond | don_angle | Angle at the donor |
hydrogen_bond, water_bridge | protisdon | Is protein the donor? |
hydrogen_bond, water_bridge, halogen_bond | donoridx/don_idx | Atom ID of the donor atom |
hydrogen_bond, water_bridge, halogen_bond | donortype | Atom type of the donor atom |
hydrogen_bond, water_bridge, halogen_bond | acceptoridx/acc_idx | Atom ID of the acceptor atom |
hydrogen_bond, water_bridge, halogen_bond | acceptortype | Atom type of the acceptor atom |
water_bridge | dist_a-w | Distance between the acceptor and interacting atom from water |
water_bridge | dist_d-w | Distance between the donor and water interacting atom from water |
water_bridge | water_angle | Angle at the interacting water atoms
water_bridge | water_idx | Atom ID of the water oxygen atom |
halogen_bond | acc_angle | Angle at the aceptor |
salt_bridge | protispos | Does the protein carry the positive charge? |
salt_bridge, pi_cation_interaction | lig_group | Functional group in the ligand |
salt_bridge, pi_stack | lig_idx_list | List of atom IDs from the functional group in the ligand |
pi_stack | cent_dist | Distance between the ring centers |
pi_stack | angle | Angle between the ring planes |
pi_stack, pi_cation_interaction | offset | Offset between the interacting groups |
pi_stack | type | Stacking type (Perpendicular or T-Shaped) |
pi_cation_interaction | protcharged | Does the protein provide the charge? |
metal_complex| metalcoo | Coordinates of interacting metal atom |
metal_complex| targetcoo | Coordinates of interacting protein atom or interaction center in the chelating target group (in protein or ligand) |
metal_complex | metal_idx | Atom ID of the metal ion |
metal_complex | metal_type | Atom type of the metal |
metal_complex | target_idx | Atom ID of the target interacting atom |
metal_complex | target_type | Atom type of the target interacting atom |
metal_complex | coordination | Metal coordination number |
metal_complex | location | Location of the target group |
metal_complex | rms | RMS of the geometry fit |
metal_complex | geometry | Metal coordination type |
metal_complex | complexnum | Continous numbering for the metal complex |

## Testing
To run the tests, run the following set of commands.

```bash
cd ~/pliptool/plip/test
python -m unittest test*
```
