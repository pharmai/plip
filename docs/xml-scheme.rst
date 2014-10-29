The XML files contain information on non-covalent interactions between proteins and ligands.
One protein can have many binding sites for ligands, i.e. one protein structure (PDB file) can contain multiple ligands.

TAGS
====

<pdbid>
-------
Unique identifier for protein structure. Can be linked to corresponding PDB entry
e.g. http://www.rcsb.org/pdb/explore/explore.do?structureId=1vsn

<hetid>
-------
Unique identifier for ligand molecule. Can be linked to corresponding PDB entry
e.g. http://www.rcsb.org/pdb/ligand/ligandsummary.do?hetId=NFT&sid=1VSN

<chain>
-------
One protein can consist of multiple separate amino acids chains which are counted using A,B,C,...

<position>
----------
Position of ligand in PDB numbering. Combination of pdbid, hetid, chain and position gives a unique identifier for
each protein-ligand complex. Same numbering as <resnr>

<interactions>
--------------
Contains interaction for protein-ligand complex, organized by interaction type, e.g. hydrophobic interactions

<resnr>
-------
Position of amino acid in protein chain according to PDB numbering

<restype>
---------
Amino acid type in three-letter code

<dist*>
-------
Distance of interacting atoms

<*idx>
------
Atom ID in PDB structure (all atoms in one PDB file are enumerated)

<lig_idx_list>
--------------
Atom IDs if several ligand atoms are relevant for a single interaction (e.g. when forming a charge center)

<*angle>
--------
Angle between interacting groups

<protispos>, <protisdon>, <protischarged>
-----------------------------------------
Determines if the protein is positively charged, provides a donor or a charge.
Important for interactions with directionality.

<ligcoo>, <protcoo>, <itype>
----------------------------
Irrelevant for output, should not be displayed