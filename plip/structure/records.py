"""Typed immutable records shared by structure preparation and detection."""

from typing import NamedTuple, TypeAlias

import numpy as np
from numpy.typing import NDArray
from openbabel import pybel
from openbabel.openbabel import OBResidue, OBRing


Coordinate: TypeAlias = list[float] | tuple[float, float, float] | NDArray[np.floating]
LigandMember: TypeAlias = tuple[str, str, int]
Region: TypeAlias = dict[str, list[int]]
RegionPair: TypeAlias = tuple[Region, Region | None]


class CovalentLink(NamedTuple):
    id1: str
    chain1: str
    pos1: int
    conf1: str
    id2: str
    chain2: str
    pos2: int
    conf2: str


class LigandRecord(NamedTuple):
    mol: pybel.Molecule
    hetid: str
    chain: str
    position: int
    water: list[OBResidue]
    members: list[LigandMember]
    longname: str
    type: str
    atomorder: list[int] | None
    can_to_pdb: dict[int, int]
    regions: RegionPair | None


class HydrophobicAtom(NamedTuple):
    atom: pybel.Atom
    orig_atom: pybel.Atom
    orig_idx: int


class HydrogenBondAcceptor(NamedTuple):
    a: pybel.Atom
    a_orig_atom: pybel.Atom
    a_orig_idx: int
    type: str


class HydrogenBondDonor(NamedTuple):
    d: pybel.Atom | HydrophobicAtom
    d_orig_atom: pybel.Atom
    d_orig_idx: int
    h: pybel.Atom
    type: str


class AromaticRing(NamedTuple):
    atoms: list[pybel.Atom]
    orig_atoms: list[pybel.Atom]
    atoms_orig_idx: list[int]
    normal: NDArray[np.floating]
    obj: OBRing
    center: Coordinate
    type: str


class HalogenBondAcceptor(NamedTuple):
    o: pybel.Atom
    o_orig_idx: int
    y: pybel.Atom
    y_orig_idx: int


class ProteinCharge(NamedTuple):
    atoms: list[pybel.Atom]
    atoms_orig_idx: list[int]
    type: str
    center: Coordinate
    restype: str
    resnr: int
    reschain: str


class ProteinMetalBinding(NamedTuple):
    atom: pybel.Atom
    atom_orig_idx: int
    type: str
    restype: str
    resnr: int
    reschain: str
    location: str


class WaterMolecule(NamedTuple):
    oxy: pybel.Atom
    oxy_orig_idx: int


class MetalAtom(NamedTuple):
    m: pybel.Atom
    orig_m: pybel.Atom
    m_orig_idx: int


class HalogenBondDonor(NamedTuple):
    x: pybel.Atom
    orig_x: pybel.Atom
    x_orig_idx: int
    c: pybel.Atom
    c_orig_idx: list[int]


class LigandCharge(NamedTuple):
    atoms: list[pybel.Atom]
    orig_atoms: list[pybel.Atom]
    atoms_orig_idx: list[int]
    type: str
    center: Coordinate
    fgroup: str


class LigandMetalBinding(NamedTuple):
    atom: pybel.Atom
    orig_atom: pybel.Atom
    atom_orig_idx: int
    type: str
    fgroup: str
    restype: str
    resnr: int
    reschain: str
    location: str


MetalBinding: TypeAlias = ProteinMetalBinding | LigandMetalBinding


class HydrophobicInteraction(NamedTuple):
    bsatom: pybel.Atom
    bsatom_orig_idx: int
    ligatom: pybel.Atom
    ligatom_orig_idx: int
    distance: float
    restype: str
    resnr: int
    reschain: str
    restype_l: str
    resnr_l: int
    reschain_l: str


class HydrogenBond(NamedTuple):
    a: pybel.Atom
    a_orig_idx: int
    d: pybel.Atom
    d_orig_idx: int
    h: pybel.Atom
    distance_ah: float
    distance_ad: float
    angle: float
    type: str
    protisdon: bool
    resnr: int
    restype: str
    reschain: str
    resnr_l: int
    restype_l: str
    reschain_l: str
    sidechain: bool
    atype: str
    dtype: str


class PiStack(NamedTuple):
    proteinring: AromaticRing
    ligandring: AromaticRing
    distance: float
    angle: float
    offset: float
    type: str
    restype: str
    resnr: int
    reschain: str
    restype_l: str
    resnr_l: int
    reschain_l: str


class PiCationInteraction(NamedTuple):
    ring: AromaticRing
    charge: ProteinCharge | LigandCharge
    distance: float
    offset: float
    type: str
    restype: str
    resnr: int
    reschain: str
    restype_l: str
    resnr_l: int
    reschain_l: str
    protcharged: bool


class SaltBridge(NamedTuple):
    positive: ProteinCharge | LigandCharge
    negative: ProteinCharge | LigandCharge
    distance: float
    protispos: bool
    resnr: int
    restype: str
    reschain: str
    resnr_l: int
    restype_l: str
    reschain_l: str


class HalogenBond(NamedTuple):
    acc: HalogenBondAcceptor
    acc_orig_idx: int
    don: HalogenBondDonor
    don_orig_idx: int
    distance: float
    don_angle: float
    acc_angle: float
    restype: str
    resnr: int
    reschain: str
    restype_l: str
    resnr_l: int
    reschain_l: str
    donortype: str
    acctype: str
    sidechain: bool


class WaterBridge(NamedTuple):
    a: pybel.Atom
    a_orig_idx: int
    atype: str
    d: pybel.Atom
    d_orig_idx: int
    dtype: str
    h: pybel.Atom
    water: pybel.Atom
    water_orig_idx: int
    distance_aw: float
    distance_dw: float
    d_angle: float
    w_angle: float
    type: str
    resnr: int
    restype: str
    reschain: str
    resnr_l: int
    restype_l: str
    reschain_l: str
    protisdon: bool


class MetalComplex(NamedTuple):
    metal: pybel.Atom
    metal_orig_idx: int
    metal_type: str
    target: MetalBinding
    target_orig_idx: int
    target_type: str
    coordination_num: int
    distance: float
    resnr: int
    restype: str
    reschain: str
    restype_l: str
    reschain_l: str
    resnr_l: int
    location: str
    rms: float
    geometry: str
    num_partners: int
    complexnum: int


class GeometryFit(NamedTuple):
    geometry: str
    rms: float
    coordination: int
    excluded: list[int]
    diff_targets: int
