PLIP Ergebnistabellen
=====================
<Attribut> - <Name> - <Beschreibung>

Hydrophobic Interactions (hydrophobic_interactions)
---------------------------------------------------
id - # - Index
resnr -  Residue - Number of residue in PDB file
restype - AA - Amino acid type
dist -  Distance - Distance between interactions carbon atoms
ligcarbonidx - Ligand Atom - ID of ligand carbon atom
protcarbonidx - Protein Atom - ID of protein carbon atom

Hydrogen Bonds (hydrogen_bonds)
-------------------------------
id - # - Index
resnr -  Residue - Number of residue in PDB file
restype - AA - Amino acid type
dist_h-a -  Distance H-A - Distance between hydrogen and acceptor atoms
dist_d-a - Distance D-A - Distance between donor and acceptor atoms
don_angle - Donor Angle - Angle between donor, acceptor and hydrogen atoms
protisdon - Is protein donor? - Does the protein provide the donor group?
donoridx - Donor Atom - ID of donor atom
acceptoridx - Acceptor Atom - ID of acceptor atom

Water Bridges (water_bridges)
-----------------------------
id - # - Index
resnr -  Residue - Number of residue in PDB file
restype - AA - Amino acid type
dist_a-w - Distance A-W - Distance between acceptor atom and water oxygen
dist_d-w - Distance D-W - Distance between donor atom and water oxygen
don_angle - Donor Angle - Angle between donor, hydrogen and water oxygen atoms
water_angle - Water Angle - Angle between acceptor, water oxygen and donor hydrogen atoms
protisdon - Is protein donor? Does the protein provide the donor group?
donor_idx - Donor Atom - ID of donor atom
acceptor_idx - Acceptor Atom - ID of acceptor atom
water_idx - Water Atom - ID of water oxygen atom

π-Stacking (pi_stacks)
----------------------
id - # - Index
resnr -  Residue - Number of residue in PDB file
restype - AA - Amino acid type
centdist - Distance - Distance between ring centers
angle - Angle - Angle between ring planes
offset - Offset - Offset between ring centers
type - Type - Type of stacking: Parallel (P) or T-Shaped (T)
lig_idx_list - Ligand Atoms - List of ligand atom IDs constituting the aromatic ring

Salt Bridges (salt_bridges)
---------------------------
id - # - Index
resnr -  Residue - Number of residue in PDB file
restype - AA - Amino acid type
dist - Distance -  Distance between centers of charge
protispos - Is protein positive? - Does the protein provide the positive charge?
lig_group - Ligand Group - Functional group in the ligand providing the charge
lig_idx_list - Ligand Atoms - List of ligand atom IDs constituting the aromatic ring

π-Cation Interactions (pi_cation_interactions)
----------------------------------------------
id - # - Index
resnr -  Residue - Number of residue in PDB file
restype - AA - Amino acid type
dist - Distance - Distance between center of charge and aromatic ring center
offset - Offset - Offset between center of charge and aromatic ring center
protcharged - Is protein charged? - Does the protein provide the charge?
lig_group - Ligand Group - Functional group in the ligand providing the charge
lig_idx_list - Ligand Atoms - List of ligand atom IDs constituting the aromatic ring or providing the charge

Halogen Bonds (halogen_bonds)
-----------------------------
id - # - Index
resnr -  Residue - Number of residue in PDB file
restype - AA - Amino acid type
dist - Distance - Distance between acceptor oxygen and donor halogen atom
don_angle - Donor Angle - Angle between donor halogen, donor carbon and acceptor oxygen atom
acc_angle - Acceptor Angle - Angle between acceptor oxygen, acceptor base and donor halogen atom
donortype - Halogen - Type of halogen participating in the interaction
don_idx - Donor Atom - ID of donor atom
acc_idx - Acceptor Atom - ID of acceptor atom