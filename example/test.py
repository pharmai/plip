from plip.structure.preparation import PDBComplex

my_mol = PDBComplex()
my_mol.load_pdb('./tmp/1EVE.pdb') # Load the PDB file into PLIP class
print(my_mol) # Shows name of structure and ligand binding sites
my_bsid = 'E20:A:2001' # Unique binding site identifier (HetID:Chain:Position)
my_mol.analyze()
my_interactions = my_mol.interaction_sets[my_bsid] # Contains all interaction data

# Now print numbers of all residues taking part in pi-stacking
print([pistack.resnr for pistack in my_interactions.pistacking]) # Prints [84, 129]
