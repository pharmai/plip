PLIP
====

Protein-Ligand Interaction Profiler (PLIP) - Analyze non-covalent protein-ligand interactions in 3D structures

Changelog
---------
#### 1.0.2
* __Automatic grouping of composite ligands (e.g. polysaccharides)__
* __Proper handling of alternative conformations in PDB structures__
* __Exclusion of modified residues as ligands__
* Adds atom type description in the output files
* Basic support for usage on Windows (without multithreading)
* Option to turn multithreading off by setting maxthreads to 0
* Improved detection of hydrogen bond donors in ligands
* Adaption of standard parameters
* Fixes a bug in PyMOL visualization script leading to missing or wrong interactions with pseudoatoms
* Fixes a bug leading to duplicate or triplicate detection of identical pi-cation interactions with guanidine
* Adds now unit tests
* Small changes to existing unit tests for new features

#### 1.0.1
* __Option to change detection thresholds permanently or for single runs__
* Option to (de)activate output for images, PyMOL session files and XML files
* Changed standard behaviour to output of RST report only
* Information on sidechain/backbone hydrogen bond type
* Sorted output
* Detection of more flavors of halogen bonds
* Fixed bug leading to duplicate interactions with quartamine groups

#### 1.0.0
* __Initial Release__
