# PLIP

**Protein-Ligand Interaction Profiler (PLIP)**

 Analyze non-covalent protein-ligand interactions in 3D structures

 <img src="pliplogo.png"  alt="PLIP Logo" height="100">

## How to use PLIP

This README provides instructions for setup and using basic functions of PLIP.
For more details, see the [Documentation](DOCUMENTATION.md).

#### 1. Clone the repository

Open a new system terminal and clone this repository using
```bash
git clone https://github.com/ssalentin/plip.git ~/pliptool
```
#### 2. Run PLIP

##### As a command line tool

Run the `plipcmd.py` script inside the PLIP folder to detect, report, and visualize interactions. The following example creates a PYMOL visualization for the interactions
between the inhibitor NFT and its target protein in the PDB structure 1VSN.

```bash
alias plip='python ~/pliptool/plip/plipcmd.py'
mkdir /tmp/1vsn && cd /tmp/1vsn
plip -i 1vsn -yv
pymol 1VSN_NFT_A_283.pse
```

##### As a python library

In your terminal, add the PLIP repository to your PYTHONPATH variable.
For our example, we also download a PDB file for testing.
```bash
export PYTHONPATH=~/pliptool:${PYTHONPATH}
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
print [pistack.resnr for pistack in my_interactions.pistacking] # Prints [84, 129]
```

#### 3. View and process the results
PLIP offers various output formats, ranging from renderes images and PyMOL session files to human-readable text files and XML files.
By default, all files are deposited in the working directory unless and output path is provided.
For a full documentation of running options and output formats, please refear to the documentation.

## Versions and Branches
For production environments, you should use the latest versioned commit from the stable branch.
Newer commits from the stable and development branch may contain new but untested and not documented features.

## Requirements
Previous to using PLIP, make sure you have the following tools and libraries installed:
* Python 2.7.x
* OpenBabel >=2.3.2
* PyMOL 1.7.x with Python bindings (optional, for visualization only)
* Imagemagick >=6.9.x (optional)

## Contributions
Sebastian Salentin

Joachim Haupt joachim.haupt (at) biotec.tu-dresden.de | https://github.com/vjhaupt

Melissa F. Adasme Mora melissa.adasme (at) biotec.tu-dresden.de

## PLIP Web Server
Visit our PLIP Web Server on https://plip.biotec.tu-dresden.de/plip-web

## Contact Me
Do you have feature requests, found a bug or want to use `PLIP` in your project?
Write me an email to joachim.haupt (at) biotec.tu-dresden.de

## License Information
PLIP is published under the Apache License. For more information, please read the LICENSE.txt file.
Using PLIP in your commercial or non-commercial project is generally possible when giving a proper reference to this project and the publication in NAR.
If you are unsure about usage options, don't hesitate to contact me.

## Citation Information
If you are using PLIP in your work, please cite
> Salentin,S. et al. PLIP: fully automated protein-ligand interaction profiler.
> Nucl. Acids Res. (1 July 2015) 43 (W1): W443-W447. doi: 10.1093/nar/gkv315

## FAQ
> I try to run PLIP, but I'm getting an error message saying:
> ValueError: [...] is not a recognised Open Babel descriptor type

Make sure Open Babel is correctly installed. This error can occur if the installed Python bindings don't match the OpenBabel version on your machine.
We don't offer technical support for installation of third-party packages.
For an instruction how to install Open Babel, please refer to their [website](https://openbabel.org/docs/dev/Installation/install.html).
