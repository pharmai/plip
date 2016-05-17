# PLIP

**Protein-Ligand Interaction Profiler (PLIP)**


 Analyze non-covalent protein-ligand interactions in 3D structures

## How to use PLIP

#### 1. Clone the repository

Open a new system terminal and clone this repository using
```bash
git clone https://github.com/ssalentin/plip.git ~/pliptool
```
#### 2. Run PLIP

##### As a command line tool

Run the `plipcmd` script inside the PLIP folder to detect, report, and visualize interactions. The following example creates a PYMOL visualization for the interactions
between the inhibitor NFT and its target protein in the PDB structure 1VSN.

```bash
alias plip='~/pliptool/plip/plipcmd'
mkdir /tmp/1vsn && cd /tmp/1vsn
plip -i 1vsn -yv
pymol 1VSN_NFT_A_283.pse
```

##### As a python library

In your terminal, add the PLIP repository to your PYTHONPATH variable.
For our example, we also download a PDB file for testing.
```bash
export PYTHONPATH=~/pliptool/plip:${PYTHONPATH}
cd /tmp && wget http://files.rcsb.org/download/1EVE.pdb
python
```
In python, import the PLIP modules, load a PDB structure and run the analysis.
This small example shows how to print all numbers of residues involved in pi-stacking:

```python
from plip.modules.preparation import PDBComplex

my_mol = PDBComplex()
my_mol.load_pdb('/tmp/1EVE.pdb') # Load the PDB file into PLIP class
print my_mol # Shows name of structure and ligand binding sites
my_bsid = 'E20:A:2001' # Unique binding site identifier (HetID:Chain:Position)
my_mol.analyze()
my_interactions = my_mol.interaction_sets[my_bsid] # Contains all interaction data

# Now print numbers of all residues taking part in pi-stacking
print [pistack.resnr for pistack in s.pistacking] # Prints [84, 129]
```

#### 3. View and process the results

Interpretation of output files and running options, read README

## Versions and Branches
The latest commits may contain newer, but untested features. The last version marked as release is always stable.

## Installation
Previous to the installation of PLIP, make sure you have the following tools and libraries installed:
* Python 2.x
* OpenBabel >=2.3.2
* PyMOL 1.7.x with Python bindings
* pip (Python package managaer)
* Imagemagick (optional, needed for automatic cropping of images)

To install PLIP, simply run the following command using a terminal
> python setup.py install

## Contributions
Sebastian Salentin sebastian.salentin (at) biotec.tu-dresden.de

Joachim Haupt joachim.haupt (at) biotec.tu-dresden.de | @vjhaupt

Melissa F. Adasme Mora melissa.adasme (at) biotec.tu-dresden.de

## PLIP Web Server
Visit our PLIP Web Server on http://projects.biotec.tu-dresden.de/plip-web

## Contact Me
Feature requests, bugs, want to use in project ...
Questions or comments about `PLIP`? Write me an email to sebastian.salentin (at) biotec.tu-dresden.de

## License Information

## Citation Information
If you are using PLIP in your work, please cite
> Salentin,S. et al. PLIP: fully automated protein-ligand interaction profiler.
> Nucl. Acids Res. (1 July 2015) 43 (W1): W443-W447. doi: 10.1093/nar/gkv315
