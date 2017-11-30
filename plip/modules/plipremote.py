"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
plipremote.py - Modules involved in multiprocessing and remote computation.
"""

# Python Standard Library
from collections import namedtuple

hbonds_info = namedtuple('hbonds_info', 'ldon_id lig_don_id prot_acc_id pdon_id prot_don_id lig_acc_id')
hydrophobic_info = namedtuple('hydrophobic_info', 'bs_ids lig_ids pairs_ids')
halogen_info = namedtuple('halogen_info', 'don_id acc_id')
pistack_info = namedtuple('pistack_info', 'proteinring_atoms, proteinring_center ligandring_atoms '
                                          'ligandring_center type')
pication_info = namedtuple('pication_info', 'ring_center charge_center ring_atoms charge_atoms, protcharged')
sbridge_info = namedtuple('sbridge_info', 'positive_atoms negative_atoms positive_center negative_center protispos')
wbridge_info = namedtuple('wbridge_info', 'don_id acc_id water_id protisdon')
metal_info = namedtuple('metal_info', 'metal_id, target_id location')

class VisualizerData:
    """Contains all information on a complex relevant for visualization. Can be pickled"""
    def __init__(self, mol, site):
        pcomp = mol
        pli = mol.interaction_sets[site]
        ligand = pli.ligand



        # General Information
        self.lig_members = sorted(pli.ligand.members)
        self.sourcefile = pcomp.sourcefiles['pdbcomplex']
        self.corrected_pdb = pcomp.corrected_pdb
        self.pdbid = mol.pymol_name
        self.hetid = ligand.hetid
        self.ligandtype = ligand.type
        self.chain = ligand.chain if not ligand.chain == "0" else ""  # #@todo Fix this
        self.position = str(ligand.position)
        self.uid = ":".join([self.hetid, self.chain, self.position])
        self.outpath = mol.output_path
        self.metal_ids = [x.m_orig_idx for x in pli.ligand.metals]
        self.unpaired_hba_idx = pli.unpaired_hba_orig_idx
        self.unpaired_hbd_idx = pli.unpaired_hbd_orig_idx
        self.unpaired_hal_idx = pli.unpaired_hal_orig_idx

        # Information on Interactions

        # Hydrophobic Contacts
        # Contains IDs of contributing binding site, ligand atoms and the pairings
        hydroph_pairs_id = [(h.bsatom_orig_idx, h.ligatom_orig_idx) for h in pli.hydrophobic_contacts]
        self.hydrophobic_contacts = hydrophobic_info(bs_ids=[hp[0] for hp in hydroph_pairs_id],
                                                     lig_ids=[hp[1] for hp in hydroph_pairs_id],
                                                     pairs_ids=hydroph_pairs_id)

        # Hydrogen Bonds
        # #@todo Don't use indices, simplify this code here
        hbonds_ldon, hbonds_pdon = pli.hbonds_ldon, pli.hbonds_pdon
        hbonds_ldon_id = [(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_ldon]
        hbonds_pdon_id = [(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_pdon]
        self.hbonds = hbonds_info(ldon_id=[(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_ldon],
                                  lig_don_id=[hb[1] for hb in hbonds_ldon_id],
                                  prot_acc_id=[hb[0] for hb in hbonds_ldon_id],
                                  pdon_id=[(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_pdon],
                                  prot_don_id=[hb[1] for hb in hbonds_pdon_id],
                                  lig_acc_id=[hb[0] for hb in hbonds_pdon_id])

        # Halogen Bonds
        self.halogen_bonds = [halogen_info(don_id=h.don_orig_idx, acc_id=h.acc_orig_idx)
                              for h in pli.halogen_bonds]

        # Pistacking
        self.pistacking = [pistack_info(proteinring_atoms=pistack.proteinring.atoms_orig_idx,
                                        proteinring_center=pistack.proteinring.center,
                                        ligandring_atoms=pistack.ligandring.atoms_orig_idx,
                                        ligandring_center=pistack.ligandring.center,
                                        type=pistack.type) for pistack in pli.pistacking]

        # Pi-cation interactions
        self.pication = [pication_info(ring_center=picat.ring.center,
                                       charge_center=picat.charge.center,
                                       ring_atoms=picat.ring.atoms_orig_idx,
                                       charge_atoms=picat.charge.atoms_orig_idx,
                                       protcharged=picat.protcharged)
                         for picat in pli.pication_paro+pli.pication_laro]

        # Salt Bridges
        self.saltbridges = [sbridge_info(positive_atoms=sbridge.positive.atoms_orig_idx,
                                         negative_atoms=sbridge.negative.atoms_orig_idx,
                                         positive_center=sbridge.positive.center,
                                         negative_center=sbridge.negative.center,
                                         protispos=sbridge.protispos)
                            for sbridge in pli.saltbridge_lneg+pli.saltbridge_pneg]

        # Water Bridgese('wbridge_info', 'don_id acc_id water_id protisdon')
        self.waterbridges = [wbridge_info(don_id=wbridge.d_orig_idx,
                                          acc_id=wbridge.a_orig_idx,
                                          water_id=wbridge.water_orig_idx,
                                          protisdon=wbridge.protisdon) for wbridge in pli.water_bridges]

        # Metal Complexes
        self.metal_complexes = [metal_info(metal_id=metalc.metal_orig_idx,
                                           target_id=metalc.target_orig_idx,
                                           location=metalc.location) for metalc in pli.metal_complexes]
