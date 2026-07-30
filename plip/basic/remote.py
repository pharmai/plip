from typing import NamedTuple

from plip.structure.preparation import PDBComplex
from plip.structure.records import Coordinate


class HydrogenBondVisualization(NamedTuple):
    ldon_id: list[tuple[int, int]]
    lig_don_id: list[int]
    prot_acc_id: list[int]
    pdon_id: list[tuple[int, int]]
    prot_don_id: list[int]
    lig_acc_id: list[int]


class HydrophobicVisualization(NamedTuple):
    bs_ids: list[int]
    lig_ids: list[int]
    pairs_ids: list[tuple[int, int]]


class HalogenVisualization(NamedTuple):
    don_id: int
    acc_id: int


class PiStackVisualization(NamedTuple):
    proteinring_atoms: list[int]
    proteinring_center: Coordinate
    ligandring_atoms: list[int]
    ligandring_center: Coordinate
    type: str


class PiCationVisualization(NamedTuple):
    ring_center: Coordinate
    charge_center: Coordinate
    ring_atoms: list[int]
    charge_atoms: list[int]
    protcharged: bool


class SaltBridgeVisualization(NamedTuple):
    positive_atoms: list[int]
    negative_atoms: list[int]
    positive_center: Coordinate
    negative_center: Coordinate
    protispos: bool


class WaterBridgeVisualization(NamedTuple):
    don_id: int
    acc_id: int
    water_id: int
    protisdon: bool


class MetalVisualization(NamedTuple):
    metal_id: int
    target_id: int
    location: str


# Backwards-compatible aliases for the historical public record names.
hbonds_info = HydrogenBondVisualization
hydrophobic_info = HydrophobicVisualization
halogen_info = HalogenVisualization
pistack_info = PiStackVisualization
pication_info = PiCationVisualization
sbridge_info = SaltBridgeVisualization
wbridge_info = WaterBridgeVisualization
metal_info = MetalVisualization


class VisualizerData:
    """Contains all information on a complex relevant for visualization. Can be pickled"""

    def __init__(self, mol: PDBComplex, site: str) -> None:
        pcomp = mol
        pli = mol.interaction_sets[site]
        ligand = pli.ligand

        # General Information
        self.lig_members = sorted(pli.ligand.members)
        self.source_pdb_file_content = pcomp.sourcefiles['pdbstring'] # store pdb file content as string
        self.corrected_pdb = pcomp.corrected_pdb
        self.pdbid = mol.pymol_name
        self.hetid = ligand.hetid
        self.ligandtype = ligand.type
        self.regions = ligand.regions
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
        self.hydrophobic_contacts = HydrophobicVisualization(bs_ids=[hp[0] for hp in hydroph_pairs_id],
                                                             lig_ids=[hp[1] for hp in hydroph_pairs_id],
                                                             pairs_ids=hydroph_pairs_id)

        # Hydrogen Bonds
        # #@todo Don't use indices, simplify this code here
        hbonds_ldon, hbonds_pdon = pli.hbonds_ldon, pli.hbonds_pdon
        hbonds_ldon_id = [(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_ldon]
        hbonds_pdon_id = [(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_pdon]
        self.hbonds = HydrogenBondVisualization(
            ldon_id=[(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_ldon],
            lig_don_id=[hb[1] for hb in hbonds_ldon_id],
            prot_acc_id=[hb[0] for hb in hbonds_ldon_id],
            pdon_id=[(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_pdon],
            prot_don_id=[hb[1] for hb in hbonds_pdon_id],
            lig_acc_id=[hb[0] for hb in hbonds_pdon_id],
        )

        # Halogen Bonds
        self.halogen_bonds = [HalogenVisualization(don_id=h.don_orig_idx, acc_id=h.acc_orig_idx)
                              for h in pli.halogen_bonds]

        # Pistacking
        self.pistacking = [PiStackVisualization(proteinring_atoms=pistack.proteinring.atoms_orig_idx,
                                                proteinring_center=pistack.proteinring.center,
                                                ligandring_atoms=pistack.ligandring.atoms_orig_idx,
                                                ligandring_center=pistack.ligandring.center,
                                                type=pistack.type) for pistack in pli.pistacking]

        # Pi-cation interactions
        self.pication = [PiCationVisualization(ring_center=picat.ring.center,
                                               charge_center=picat.charge.center,
                                               ring_atoms=picat.ring.atoms_orig_idx,
                                               charge_atoms=picat.charge.atoms_orig_idx,
                                               protcharged=picat.protcharged)
                         for picat in pli.pication_paro + pli.pication_laro]

        # Salt Bridges
        self.saltbridges = [SaltBridgeVisualization(positive_atoms=sbridge.positive.atoms_orig_idx,
                                                    negative_atoms=sbridge.negative.atoms_orig_idx,
                                                    positive_center=sbridge.positive.center,
                                                    negative_center=sbridge.negative.center,
                                                    protispos=sbridge.protispos)
                            for sbridge in pli.saltbridge_lneg + pli.saltbridge_pneg]

        # Water Bridgese('wbridge_info', 'don_id acc_id water_id protisdon')
        self.waterbridges = [WaterBridgeVisualization(don_id=wbridge.d_orig_idx,
                                                      acc_id=wbridge.a_orig_idx,
                                                      water_id=wbridge.water_orig_idx,
                                                      protisdon=wbridge.protisdon)
                             for wbridge in pli.water_bridges]

        # Metal Complexes
        self.metal_complexes = [MetalVisualization(metal_id=metalc.metal_orig_idx,
                                                   target_id=metalc.target_orig_idx,
                                                   location=metalc.location)
                                for metalc in pli.metal_complexes]
