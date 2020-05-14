# Protein-Ligand Interaction Profiler (PLIP)

![PLIP Build](https://github.com/pharmai/plip/workflows/PLIP%20Build/badge.svg)
![GitHub](https://img.shields.io/github/license/pharmai/plip?style=social)
![GitHub All Releases](https://img.shields.io/github/downloads/pharmai/plip/total?style=social)
![Docker Image Size (tag)](https://img.shields.io/docker/image-size/pharmai/plip/latest?style=social)

Analyze noncovalent protein-ligand interactions in 3D structures with ease.

<img src="pliplogo.png"  alt="PLIP Logo" height="100">

## Quickstart
If you have Docker installed, you can run a PLIP analysis for the structure `1vsn` with the following shell command:

On Linux / MacOS:
```bash
$ docker run --rm \
    -v ${PWD}:/results \
    -w /results \
    -u $(id -u ${USER}):$(id -g ${USER}) \
    pharmai/plip:latest -i 1s3v -yv
```

On Windows:
```bash
$ docker run --rm \
    -v ${PWD}:/results \
    -w /results \
    -u $(id -u ${USER}):$(id -g ${USER}) \
    pharmai/plip:latest -i 1s3v -yv
```

The equivalent command for our pre-built [Singularity](https://singularity.lbl.gov/) image for Linux (available under [Releases](https://github.com/pharmai/plip/releases)) is as follows:

```
./plip.simg -i 1vsn -yv
```

Singularity allows to use PLIP with ease in HPC environments.

## How to use PLIP

This README provides instructions for setup and using basic functions of PLIP.
For more details, see the [Documentation](DOCUMENTATION.md).

### 1. (optional) Clone the repository

Open a new system terminal and clone this repository using
```bash
git clone https://github.com/ssalentin/plip.git ~/pliptool
```

### 2. Install PLIP

#### Containerized Image (no installation) :exclamation:
We ship PLIP as a pre-built Docker container, available on the Docker Hub ([https://hub.docker.com/](https://hub.docker.com/)) or as pre-built Singularity image.

#### Python Module
If you cannot use the Docker bundle or want to use PLIP sources, make sure you have the following requirements installed:
- Python >= 3.6.9
- OpenBabel >= 3.0.0
- PyMOL >= 2.3.0 with Python bindings (optional, for visualization only)
- ImageMagick >= 6.9 (optional)

Set your `PYTHONPATH` environment variable to the root directory of this repository.

### 3. Run PLIP

#### As a command line tool

Run the `plipcmd.py` script inside the PLIP folder to detect, report, and visualize interactions. The following example creates a PYMOL visualization for the interactions
between the inhibitor NFT and its target protein in the PDB structure 1VSN.

```bash
alias plip='python ~/pliptool/plip/plipcmd.py'
mkdir /tmp/1vsn && cd /tmp/1vsn
plip -i 1vsn -yv
pymol 1VSN_NFT_A_283.pse
```

#### As a python library
In your terminal, add the PLIP repository to your PYTHONPATH variable. For our example, we also download a PDB file for testing.
```bash
export PYTHONPATH=~/plip:${PYTHONPATH}
cd /tmp && wget http://files.rcsb.org/download/1EVE.pdb
python
```
In python, import the PLIP modules, load a PDB structure and run the analysis.
This small example shows how to print all numbers of residues involved in pi-stacking:

```python
from plip.structure.preparation import PDBComplex

my_mol = PDBComplex()
my_mol.load_pdb('/tmp/1EVE.pdb') # Load the PDB file into PLIP class
print(my_mol) # Shows name of structure and ligand binding sites
my_bsid = 'E20:A:2001' # Unique binding site identifier (HetID:Chain:Position)
my_mol.analyze()
my_interactions = my_mol.interaction_sets[my_bsid] # Contains all interaction data

# Now print numbers of all residues taking part in pi-stacking
print([pistack.resnr for pistack in my_interactions.pistacking]) # Prints [84, 129]
```

#### 3. View and process the results
PLIP offers various output formats, ranging from renderes images and PyMOL session files to human-readable text files and XML files.
By default, all files are deposited in the working directory unless and output path is provided.
For a full documentation of running options and output formats, please refear to the documentation.

## Versions and Branches
For production environments, you should use the latest tagged commit from the `master` branch or refer to the  [Releases](https://github.com/pharmai/plip/releases)) page. Newer commits from the `master` and `development` branch may contain new but untested and not documented features.

## Contributors
- Sebastian Salentin (original author) | [github.com/ssalentin](https://github.com/ssalentin)
- Joachim Haupt | [github.com/vjhaupt](https://github.com/vjhaupt)
- Melissa F. Adasme Mora |  [github.com/madasme](https://github.com/madasme)
- Alexandre Mestiashvili | [github.com/mestia](https://github.com/mestia)
- Christoph Leberecht  | [github.com/cleberecht](https://github.com/cleberecht)
- Florian Kaiser  | [github.com/fkaiserbio](https://github.com/fkaiserbio)

## PLIP Web Server
Visit our PLIP Web Server on https://plip.biotec.tu-dresden.de/plip-web

Do you have feature requests, found a bug or want to use  PLIP in your project?

Write an email to `contact@pharm.ai`.

## License Information
PLIP is published under the GNU GPLv2. For more information, please read the LICENSE.txt file.
Using PLIP in your commercial or non-commercial project is generally possible when giving a proper reference to this project and the publication in NAR.
If you are unsure about usage options, don't hesitate to contact me.

## Citation Information
If you are using PLIP in your work, please cite
> Salentin,S. et al. PLIP: fully automated protein-ligand interaction profiler.
> Nucl. Acids Res. (1 July 2015) 43 (W1): W443-W447. doi: 10.1093/nar/gkv315

## FAQ
> I try to run PLIP, but I'm getting an error message saying:
> ValueError: [...] is not a recognised Open Babel descriptor type
>
Make sure OpenBabel is correctly installed. This error can occur if the installed Python bindings don't match the OpenBabel version on your machine.
We don't offer technical support for installation of third-party packages.
For an instruction how to install Open Babel, please refer to their [website](https://openbabel.org/docs/dev/Installation/install.html).

> I'm unsure on how to run PLIP and don't have much Linux experience.
>
You should consider running PLIP as Docker image, as we describe above.

> PLIP is reporting different interactions on several runs!
>
Due to the non-deterministic nature on how hydrogen atoms can be added to the input structure, it cannot be guaranteed that each run returns exactly the same set of interactions. If you want to make sure to achieve consistent results, you can:

- protonate the input structure once with PLIP or your tool of preference
- run PLIP with `--nohydro`

## Contact / Maintainer
As of April 2020 PLIP is now officially maintained by [PharmAI GmbH](https://pharm.ai). Do you have feature requests, found a bug or want to use  PLIP in your project? Commercial support is available upon request.

 ![](https://www.pharm.ai/wp-content/uploads/2020/04/PharmAI_logo_color_no_slogan_500px.png)
 
 Please get in touch: `contact@pharm.ai`
