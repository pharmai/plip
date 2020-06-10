from typing import Optional, Any, List, Dict, TypeVar, Callable, Type, cast
from datetime import datetime
import dateutil.parser
import json
import hashlib
import plip.basic.config
from plip.basic import config

from plip.structure.preparation import PDBComplex

T = TypeVar("T")


def from_int(x: Any) -> int:
    assert isinstance(x, int) and not isinstance(x, bool)
    return x


def from_none(x: Any) -> Any:
    assert x is None
    return x


def from_union(fs, x):
    for f in fs:
        try:
            return f(x)
        except:
            pass
    assert False


def from_str(x: Any) -> str:
    assert isinstance(x, str)
    return x


def from_bool(x: Any) -> bool:
    assert isinstance(x, bool)
    return x


def from_float(x: Any) -> float:
    assert isinstance(x, (float, int)) and not isinstance(x, bool)
    return float(x)


def to_float(x: Any) -> float:
    assert isinstance(x, float)
    return x


def from_list(f: Callable[[Any], T], x: Any) -> List[T]:
    assert isinstance(x, list)
    return [f(y) for y in x]


def to_class(c: Type[T], x: Any) -> dict:
    assert isinstance(x, c)
    return cast(Any, x).to_dict()


def from_dict(f: Callable[[Any], T], x: Any) -> Dict[str, T]:
    assert isinstance(x, dict)
    return {k: f(v) for (k, v) in x.items()}


def from_datetime(x: Any) -> datetime:
    return dateutil.parser.parse(x)


class BindingSiteMember:
    binding_site_member_id: Optional[int]
    amino_acid: Optional[str]
    ligand_contact: Optional[bool]
    minimal_distance: Optional[float]
    chain_id: Optional[str]
    residue_number: Optional[int]
    insertion_code: Optional[str]

    def __init__(self, binding_site_member_id: Optional[int], amino_acid: Optional[str], ligand_contact: Optional[bool],
                 minimal_distance: Optional[float], chain_id: Optional[str], residue_number: Optional[int],
                 insertion_code: Optional[str]) -> None:
        self.binding_site_member_id = binding_site_member_id
        self.amino_acid = amino_acid
        self.ligand_contact = ligand_contact
        self.minimal_distance = minimal_distance
        self.chain_id = chain_id
        self.residue_number = residue_number
        self.insertion_code = insertion_code

    @staticmethod
    def from_dict(obj: Any) -> 'BindingSiteMember':
        assert isinstance(obj, dict)
        binding_site_member_id = from_union([from_int, from_none], obj.get("binding_site_member_id"))
        amino_acid = from_union([from_str, from_none], obj.get("amino_acid"))
        ligand_contact = from_union([from_bool, from_none], obj.get("ligand_contact"))
        minimal_distance = from_union([from_float, from_none], obj.get("minimal_distance"))
        chain_id = from_union([from_str, from_none], obj.get("chain_id"))
        residue_number = from_union([from_int, from_none], obj.get("residue_number"))
        insertion_code = from_union([from_str, from_none], obj.get("insertion_code"))
        return BindingSiteMember(binding_site_member_id, amino_acid, ligand_contact, minimal_distance, chain_id,
                                 residue_number, insertion_code)

    def to_dict(self) -> dict:
        result: dict = {}
        result["binding_site_member_id"] = from_union([from_int, from_none], self.binding_site_member_id)
        result["amino_acid"] = from_union([from_str, from_none], self.amino_acid)
        result["ligand_contact"] = from_union([from_bool, from_none], self.ligand_contact)
        result["minimal_distance"] = from_union([to_float, from_none], self.minimal_distance)
        result["chain_id"] = from_union([from_str, from_none], self.chain_id)
        result["residue_number"] = from_union([from_int, from_none], self.residue_number)
        result["insertion_code"] = from_union([from_str, from_none], self.insertion_code)
        return result


class HalogenBond:
    binding_site_member_id: Optional[int]
    ligand_member_id: Optional[int]
    interacting_with_sidechain: Optional[bool]
    distance: Optional[float]
    donor_angle: Optional[float]
    acceptor_angle: Optional[float]
    ligand_atom_id: Optional[int]
    ligand_group_type: Optional[str]
    protein_atom_id: Optional[int]
    protein_group_type: Optional[str]
    ligand_representative_position: Optional[List[float]]
    protein_representative_position: Optional[List[float]]

    def __init__(self, binding_site_member_id: Optional[int], ligand_member_id: Optional[int],
                 interacting_with_sidechain: Optional[bool], distance: Optional[float], donor_angle: Optional[float],
                 acceptor_angle: Optional[float], ligand_atom_id: Optional[int], ligand_group_type: Optional[str],
                 protein_atom_id: Optional[int], protein_group_type: Optional[str],
                 ligand_representative_position: Optional[List[float]],
                 protein_representative_position: Optional[List[float]]) -> None:
        self.binding_site_member_id = binding_site_member_id
        self.ligand_member_id = ligand_member_id
        self.interacting_with_sidechain = interacting_with_sidechain
        self.distance = distance
        self.donor_angle = donor_angle
        self.acceptor_angle = acceptor_angle
        self.ligand_atom_id = ligand_atom_id
        self.ligand_group_type = ligand_group_type
        self.protein_atom_id = protein_atom_id
        self.protein_group_type = protein_group_type
        self.ligand_representative_position = ligand_representative_position
        self.protein_representative_position = protein_representative_position

    @staticmethod
    def from_dict(obj: Any) -> 'HalogenBond':
        assert isinstance(obj, dict)
        binding_site_member_id = from_union([from_int, from_none], obj.get("binding_site_member_id"))
        ligand_member_id = from_union([from_int, from_none], obj.get("ligand_member_id"))
        interacting_with_sidechain = from_union([from_bool, from_none], obj.get("interacting_with_sidechain"))
        distance = from_union([from_float, from_none], obj.get("distance"))
        donor_angle = from_union([from_float, from_none], obj.get("donor_angle"))
        acceptor_angle = from_union([from_float, from_none], obj.get("acceptor_angle"))
        ligand_atom_id = from_union([from_int, from_none], obj.get("ligand_atom_id"))
        ligand_group_type = from_union([from_str, from_none], obj.get("ligand_group_type"))
        protein_atom_id = from_union([from_int, from_none], obj.get("protein_atom_id"))
        protein_group_type = from_union([from_str, from_none], obj.get("protein_group_type"))
        ligand_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                    obj.get("ligand_representative_position"))
        protein_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                     obj.get("protein_representative_position"))
        return HalogenBond(binding_site_member_id, ligand_member_id, interacting_with_sidechain, distance, donor_angle,
                           acceptor_angle, ligand_atom_id, ligand_group_type, protein_atom_id, protein_group_type,
                           ligand_representative_position, protein_representative_position)

    def to_dict(self) -> dict:
        result: dict = {}
        result["binding_site_member_id"] = from_union([from_int, from_none], self.binding_site_member_id)
        result["ligand_member_id"] = from_union([from_int, from_none], self.ligand_member_id)
        result["interacting_with_sidechain"] = from_union([from_bool, from_none], self.interacting_with_sidechain)
        result["distance"] = from_union([to_float, from_none], self.distance)
        result["donor_angle"] = from_union([to_float, from_none], self.donor_angle)
        result["acceptor_angle"] = from_union([to_float, from_none], self.acceptor_angle)
        result["ligand_atom_id"] = from_union([from_int, from_none], self.ligand_atom_id)
        result["ligand_group_type"] = from_union([from_str, from_none], self.ligand_group_type)
        result["protein_atom_id"] = from_union([from_int, from_none], self.protein_atom_id)
        result["protein_group_type"] = from_union([from_str, from_none], self.protein_group_type)
        result["ligand_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                              self.ligand_representative_position)
        result["protein_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                               self.protein_representative_position)
        return result


class HydrogenBond:
    binding_site_member_id: Optional[int]
    ligand_member_id: Optional[int]
    distance: Optional[float]
    interacting_with_sidechain: Optional[bool]
    distance_hydrogen_acceptor: Optional[float]
    distance_donor_acceptor: Optional[float]
    donor_angle: Optional[float]
    protein_is_donor: Optional[bool]
    donor_atom_id: Optional[int]
    donor_group_type: Optional[str]
    acceptor_atom_id: Optional[int]
    acceptor_group_type: Optional[str]
    ligand_representative_position: Optional[List[float]]
    protein_representative_position: Optional[List[float]]

    def __init__(self, binding_site_member_id: Optional[int], ligand_member_id: Optional[int],
                 distance: Optional[float], interacting_with_sidechain: Optional[bool],
                 distance_hydrogen_acceptor: Optional[float], distance_donor_acceptor: Optional[float],
                 donor_angle: Optional[float], protein_is_donor: Optional[bool], donor_atom_id: Optional[int],
                 donor_group_type: Optional[str], acceptor_atom_id: Optional[int], acceptor_group_type: Optional[str],
                 ligand_representative_position: Optional[List[float]],
                 protein_representative_position: Optional[List[float]]) -> None:
        self.binding_site_member_id = binding_site_member_id
        self.ligand_member_id = ligand_member_id
        self.distance = distance
        self.interacting_with_sidechain = interacting_with_sidechain
        self.distance_hydrogen_acceptor = distance_hydrogen_acceptor
        self.distance_donor_acceptor = distance_donor_acceptor
        self.donor_angle = donor_angle
        self.protein_is_donor = protein_is_donor
        self.donor_atom_id = donor_atom_id
        self.donor_group_type = donor_group_type
        self.acceptor_atom_id = acceptor_atom_id
        self.acceptor_group_type = acceptor_group_type
        self.ligand_representative_position = ligand_representative_position
        self.protein_representative_position = protein_representative_position

    @staticmethod
    def from_dict(obj: Any) -> 'HydrogenBond':
        assert isinstance(obj, dict)
        binding_site_member_id = from_union([from_int, from_none], obj.get("binding_site_member_id"))
        ligand_member_id = from_union([from_int, from_none], obj.get("ligand_member_id"))
        distance = from_union([from_float, from_none], obj.get("distance"))
        interacting_with_sidechain = from_union([from_bool, from_none], obj.get("interacting_with_sidechain"))
        distance_hydrogen_acceptor = from_union([from_float, from_none], obj.get("distance_hydrogen_acceptor"))
        distance_donor_acceptor = from_union([from_float, from_none], obj.get("distance_donor_acceptor"))
        donor_angle = from_union([from_float, from_none], obj.get("donor_angle"))
        protein_is_donor = from_union([from_bool, from_none], obj.get("protein_is_donor"))
        donor_atom_id = from_union([from_int, from_none], obj.get("donor_atom_id"))
        donor_group_type = from_union([from_str, from_none], obj.get("donor_group_type"))
        acceptor_atom_id = from_union([from_int, from_none], obj.get("acceptor_atom_id"))
        acceptor_group_type = from_union([from_str, from_none], obj.get("acceptor_group_type"))
        ligand_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                    obj.get("ligand_representative_position"))
        protein_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                     obj.get("protein_representative_position"))
        return HydrogenBond(binding_site_member_id, ligand_member_id, distance, interacting_with_sidechain,
                            distance_hydrogen_acceptor, distance_donor_acceptor, donor_angle, protein_is_donor,
                            donor_atom_id, donor_group_type, acceptor_atom_id, acceptor_group_type,
                            ligand_representative_position, protein_representative_position)

    def to_dict(self) -> dict:
        result: dict = {}
        result["binding_site_member_id"] = from_union([from_int, from_none], self.binding_site_member_id)
        result["ligand_member_id"] = from_union([from_int, from_none], self.ligand_member_id)
        result["distance"] = from_union([to_float, from_none], self.distance)
        result["interacting_with_sidechain"] = from_union([from_bool, from_none], self.interacting_with_sidechain)
        result["distance_hydrogen_acceptor"] = from_union([to_float, from_none], self.distance_hydrogen_acceptor)
        result["distance_donor_acceptor"] = from_union([to_float, from_none], self.distance_donor_acceptor)
        result["donor_angle"] = from_union([to_float, from_none], self.donor_angle)
        result["protein_is_donor"] = from_union([from_bool, from_none], self.protein_is_donor)
        result["donor_atom_id"] = from_union([from_int, from_none], self.donor_atom_id)
        result["donor_group_type"] = from_union([from_str, from_none], self.donor_group_type)
        result["acceptor_atom_id"] = from_union([from_int, from_none], self.acceptor_atom_id)
        result["acceptor_group_type"] = from_union([from_str, from_none], self.acceptor_group_type)
        result["ligand_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                              self.ligand_representative_position)
        result["protein_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                               self.protein_representative_position)
        return result


class HydrophobicContact:
    binding_site_member_id: Optional[int]
    ligand_member_id: Optional[int]
    distance: Optional[float]
    ligand_atom_id: Optional[int]
    ligand_representative_position: Optional[List[float]]
    protein_atom_id: Optional[int]
    protein_representative_position: Optional[List[float]]

    def __init__(self, binding_site_member_id: Optional[int], ligand_member_id: Optional[int],
                 distance: Optional[float], ligand_atom_id: Optional[int],
                 ligand_representative_position: Optional[List[float]], protein_atom_id: Optional[int],
                 protein_representative_position: Optional[List[float]]) -> None:
        self.binding_site_member_id = binding_site_member_id
        self.ligand_member_id = ligand_member_id
        self.distance = distance
        self.ligand_atom_id = ligand_atom_id
        self.ligand_representative_position = ligand_representative_position
        self.protein_atom_id = protein_atom_id
        self.protein_representative_position = protein_representative_position

    @staticmethod
    def from_dict(obj: Any) -> 'HydrophobicContact':
        assert isinstance(obj, dict)
        binding_site_member_id = from_union([from_int, from_none], obj.get("binding_site_member_id"))
        ligand_member_id = from_union([from_int, from_none], obj.get("ligand_member_id"))
        distance = from_union([from_float, from_none], obj.get("distance"))
        ligand_atom_id = from_union([from_int, from_none], obj.get("ligand_atom_id"))
        ligand_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                    obj.get("ligand_representative_position"))
        protein_atom_id = from_union([from_int, from_none], obj.get("protein_atom_id"))
        protein_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                     obj.get("protein_representative_position"))
        return HydrophobicContact(binding_site_member_id, ligand_member_id, distance, ligand_atom_id,
                                  ligand_representative_position, protein_atom_id, protein_representative_position)

    def to_dict(self) -> dict:
        result: dict = {}
        result["binding_site_member_id"] = from_union([from_int, from_none], self.binding_site_member_id)
        result["ligand_member_id"] = from_union([from_int, from_none], self.ligand_member_id)
        result["distance"] = from_union([to_float, from_none], self.distance)
        result["ligand_atom_id"] = from_union([from_int, from_none], self.ligand_atom_id)
        result["ligand_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                              self.ligand_representative_position)
        result["protein_atom_id"] = from_union([from_int, from_none], self.protein_atom_id)
        result["protein_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                               self.protein_representative_position)
        return result


class MetalComplex:
    binding_site_member_ids: Optional[List[int]]
    ligand_member_id: Optional[int]
    coordination_order: Optional[int]
    geometry: Optional[str]
    geometry_fit_rmsd: Optional[float]
    ligand_atom_id: Optional[int]
    protein_atom_ids: Optional[List[int]]
    ligand_representative_position: Optional[List[float]]
    protein_representative_positions: Optional[List[List[float]]]

    def __init__(self, binding_site_member_ids: Optional[List[int]], ligand_member_id: Optional[int],
                 coordination_order: Optional[int], geometry: Optional[str], geometry_fit_rmsd: Optional[float],
                 ligand_atom_id: Optional[int], protein_atom_ids: Optional[List[int]],
                 ligand_representative_position: Optional[List[float]],
                 protein_representative_positions: Optional[List[List[float]]]) -> None:
        self.binding_site_member_ids = binding_site_member_ids
        self.ligand_member_id = ligand_member_id
        self.coordination_order = coordination_order
        self.geometry = geometry
        self.geometry_fit_rmsd = geometry_fit_rmsd
        self.ligand_atom_id = ligand_atom_id
        self.protein_atom_ids = protein_atom_ids
        self.ligand_representative_position = ligand_representative_position
        self.protein_representative_positions = protein_representative_positions

    @staticmethod
    def from_dict(obj: Any) -> 'MetalComplex':
        assert isinstance(obj, dict)
        binding_site_member_ids = from_union([lambda x: from_list(from_int, x), from_none],
                                             obj.get("binding_site_member_ids"))
        ligand_member_id = from_union([from_int, from_none], obj.get("ligand_member_id"))
        coordination_order = from_union([from_int, from_none], obj.get("coordination_order"))
        geometry = from_union([from_str, from_none], obj.get("geometry"))
        geometry_fit_rmsd = from_union([from_float, from_none], obj.get("geometry_fit_rmsd"))
        ligand_atom_id = from_union([from_int, from_none], obj.get("ligand_atom_id"))
        protein_atom_ids = from_union([lambda x: from_list(from_int, x), from_none], obj.get("protein_atom_ids"))
        ligand_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                    obj.get("ligand_representative_position"))
        protein_representative_positions = from_union(
            [lambda x: from_list(lambda x: from_list(from_float, x), x), from_none],
            obj.get("protein_representative_positions"))
        return MetalComplex(binding_site_member_ids, ligand_member_id, coordination_order, geometry, geometry_fit_rmsd,
                            ligand_atom_id, protein_atom_ids, ligand_representative_position,
                            protein_representative_positions)

    def to_dict(self) -> dict:
        result: dict = {}
        result["binding_site_member_ids"] = from_union([lambda x: from_list(from_int, x), from_none],
                                                       self.binding_site_member_ids)
        result["ligand_member_id"] = from_union([from_int, from_none], self.ligand_member_id)
        result["coordination_order"] = from_union([from_int, from_none], self.coordination_order)
        result["geometry"] = from_union([from_str, from_none], self.geometry)
        result["geometry_fit_rmsd"] = from_union([to_float, from_none], self.geometry_fit_rmsd)
        result["ligand_atom_id"] = from_union([from_int, from_none], self.ligand_atom_id)
        result["protein_atom_ids"] = from_union([lambda x: from_list(from_int, x), from_none], self.protein_atom_ids)
        result["ligand_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                              self.ligand_representative_position)
        result["protein_representative_positions"] = from_union(
            [lambda x: from_list(lambda x: from_list(to_float, x), x), from_none],
            self.protein_representative_positions)
        return result


class Pi:
    binding_site_member_id: Optional[int]
    ligand_member_id: Optional[int]
    ring_center_distance: Optional[float]
    ring_offset: Optional[float]
    protein_is_positive_charge: Optional[bool]
    ligand_group_type: Optional[str]
    ligand_atom_ids: Optional[List[int]]
    protein_atom_ids: Optional[List[int]]
    ligand_representative_position: Optional[List[float]]
    protein_representative_position: Optional[List[float]]
    ring_plane_angle: Optional[float]
    stacking_type: Optional[str]

    def __init__(self, binding_site_member_id: Optional[int], ligand_member_id: Optional[int],
                 ring_center_distance: Optional[float], ring_offset: Optional[float],
                 protein_is_positive_charge: Optional[bool], ligand_group_type: Optional[str],
                 ligand_atom_ids: Optional[List[int]], protein_atom_ids: Optional[List[int]],
                 ligand_representative_position: Optional[List[float]],
                 protein_representative_position: Optional[List[float]], ring_plane_angle: Optional[float],
                 stacking_type: Optional[str]) -> None:
        self.binding_site_member_id = binding_site_member_id
        self.ligand_member_id = ligand_member_id
        self.ring_center_distance = ring_center_distance
        self.ring_offset = ring_offset
        self.protein_is_positive_charge = protein_is_positive_charge
        self.ligand_group_type = ligand_group_type
        self.ligand_atom_ids = ligand_atom_ids
        self.protein_atom_ids = protein_atom_ids
        self.ligand_representative_position = ligand_representative_position
        self.protein_representative_position = protein_representative_position
        self.ring_plane_angle = ring_plane_angle
        self.stacking_type = stacking_type

    @staticmethod
    def from_dict(obj: Any) -> 'Pi':
        assert isinstance(obj, dict)
        binding_site_member_id = from_union([from_int, from_none], obj.get("binding_site_member_id"))
        ligand_member_id = from_union([from_int, from_none], obj.get("ligand_member_id"))
        ring_center_distance = from_union([from_float, from_none], obj.get("ring_center_distance"))
        ring_offset = from_union([from_float, from_none], obj.get("ring_offset"))
        protein_is_positive_charge = from_union([from_bool, from_none], obj.get("protein_is_positive_charge"))
        ligand_group_type = from_union([from_str, from_none], obj.get("ligand_group_type"))
        ligand_atom_ids = from_union([lambda x: from_list(from_int, x), from_none], obj.get("ligand_atom_ids"))
        protein_atom_ids = from_union([lambda x: from_list(from_int, x), from_none], obj.get("protein_atom_ids"))
        ligand_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                    obj.get("ligand_representative_position"))
        protein_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                     obj.get("protein_representative_position"))
        ring_plane_angle = from_union([from_float, from_none], obj.get("ring_plane_angle"))
        stacking_type = from_union([from_str, from_none], obj.get("stacking_type"))
        return Pi(binding_site_member_id, ligand_member_id, ring_center_distance, ring_offset,
                  protein_is_positive_charge, ligand_group_type, ligand_atom_ids, protein_atom_ids,
                  ligand_representative_position, protein_representative_position, ring_plane_angle, stacking_type)

    def to_dict(self) -> dict:
        result: dict = {}
        result["binding_site_member_id"] = from_union([from_int, from_none], self.binding_site_member_id)
        result["ligand_member_id"] = from_union([from_int, from_none], self.ligand_member_id)
        result["ring_center_distance"] = from_union([to_float, from_none], self.ring_center_distance)
        result["ring_offset"] = from_union([to_float, from_none], self.ring_offset)
        result["protein_is_positive_charge"] = from_union([from_bool, from_none], self.protein_is_positive_charge)
        result["ligand_group_type"] = from_union([from_str, from_none], self.ligand_group_type)
        result["ligand_atom_ids"] = from_union([lambda x: from_list(from_int, x), from_none], self.ligand_atom_ids)
        result["protein_atom_ids"] = from_union([lambda x: from_list(from_int, x), from_none], self.protein_atom_ids)
        result["ligand_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                              self.ligand_representative_position)
        result["protein_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                               self.protein_representative_position)
        result["ring_plane_angle"] = from_union([to_float, from_none], self.ring_plane_angle)
        result["stacking_type"] = from_union([from_str, from_none], self.stacking_type)
        return result


class SaltBridge:
    binding_site_member_id: Optional[int]
    ligand_member_id: Optional[int]
    distance: Optional[float]
    protein_is_positive_charge: Optional[bool]
    ligand_group_type: Optional[str]
    ligand_atom_ids: Optional[List[int]]
    ligand_representative_position: Optional[List[float]]
    protein_representative_position: Optional[List[float]]

    def __init__(self, binding_site_member_id: Optional[int], ligand_member_id: Optional[int],
                 distance: Optional[float], protein_is_positive_charge: Optional[bool],
                 ligand_group_type: Optional[str], ligand_atom_ids: Optional[List[int]],
                 ligand_representative_position: Optional[List[float]],
                 protein_representative_position: Optional[List[float]]) -> None:
        self.binding_site_member_id = binding_site_member_id
        self.ligand_member_id = ligand_member_id
        self.distance = distance
        self.protein_is_positive_charge = protein_is_positive_charge
        self.ligand_group_type = ligand_group_type
        self.ligand_atom_ids = ligand_atom_ids
        self.ligand_representative_position = ligand_representative_position
        self.protein_representative_position = protein_representative_position

    @staticmethod
    def from_dict(obj: Any) -> 'SaltBridge':
        assert isinstance(obj, dict)
        binding_site_member_id = from_union([from_int, from_none], obj.get("binding_site_member_id"))
        ligand_member_id = from_union([from_int, from_none], obj.get("ligand_member_id"))
        distance = from_union([from_float, from_none], obj.get("distance"))
        protein_is_positive_charge = from_union([from_bool, from_none], obj.get("protein_is_positive_charge"))
        ligand_group_type = from_union([from_str, from_none], obj.get("ligand_group_type"))
        ligand_atom_ids = from_union([lambda x: from_list(from_int, x), from_none], obj.get("ligand_atom_ids"))
        ligand_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                    obj.get("ligand_representative_position"))
        protein_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                     obj.get("protein_representative_position"))
        return SaltBridge(binding_site_member_id, ligand_member_id, distance, protein_is_positive_charge,
                          ligand_group_type, ligand_atom_ids, ligand_representative_position,
                          protein_representative_position)

    def to_dict(self) -> dict:
        result: dict = {}
        result["binding_site_member_id"] = from_union([from_int, from_none], self.binding_site_member_id)
        result["ligand_member_id"] = from_union([from_int, from_none], self.ligand_member_id)
        result["distance"] = from_union([to_float, from_none], self.distance)
        result["protein_is_positive_charge"] = from_union([from_bool, from_none], self.protein_is_positive_charge)
        result["ligand_group_type"] = from_union([from_str, from_none], self.ligand_group_type)
        result["ligand_atom_ids"] = from_union([lambda x: from_list(from_int, x), from_none], self.ligand_atom_ids)
        result["ligand_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                              self.ligand_representative_position)
        result["protein_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                               self.protein_representative_position)
        return result


class WaterBridge:
    binding_site_member_id: Optional[int]
    ligand_member_id: Optional[int]
    distance_acceptor_water: Optional[float]
    distance_donor_water: Optional[float]
    donor_angle: Optional[float]
    water_angle: Optional[float]
    protein_is_donor: Optional[bool]
    donor_atom_id: Optional[int]
    donor_group_type: Optional[str]
    acceptor_atom_id: Optional[int]
    acceptor_group_type: Optional[str]
    water_atom_id: Optional[str]
    water_atom_position: Optional[List[float]]
    ligand_representative_position: Optional[List[float]]
    protein_representative_position: Optional[List[float]]

    def __init__(self, binding_site_member_id: Optional[int], ligand_member_id: Optional[int],
                 distance_acceptor_water: Optional[float], distance_donor_water: Optional[float],
                 donor_angle: Optional[float], water_angle: Optional[float], protein_is_donor: Optional[bool],
                 donor_atom_id: Optional[int], donor_group_type: Optional[str], acceptor_atom_id: Optional[int],
                 acceptor_group_type: Optional[str], water_atom_id: Optional[str],
                 water_atom_position: Optional[List[float]], ligand_representative_position: Optional[List[float]],
                 protein_representative_position: Optional[List[float]]) -> None:
        self.binding_site_member_id = binding_site_member_id
        self.ligand_member_id = ligand_member_id
        self.distance_acceptor_water = distance_acceptor_water
        self.distance_donor_water = distance_donor_water
        self.donor_angle = donor_angle
        self.water_angle = water_angle
        self.protein_is_donor = protein_is_donor
        self.donor_atom_id = donor_atom_id
        self.donor_group_type = donor_group_type
        self.acceptor_atom_id = acceptor_atom_id
        self.acceptor_group_type = acceptor_group_type
        self.water_atom_id = water_atom_id
        self.water_atom_position = water_atom_position
        self.ligand_representative_position = ligand_representative_position
        self.protein_representative_position = protein_representative_position

    @staticmethod
    def from_dict(obj: Any) -> 'WaterBridge':
        assert isinstance(obj, dict)
        binding_site_member_id = from_union([from_int, from_none], obj.get("binding_site_member_id"))
        ligand_member_id = from_union([from_int, from_none], obj.get("ligand_member_id"))
        distance_acceptor_water = from_union([from_float, from_none], obj.get("distance_acceptor_water"))
        distance_donor_water = from_union([from_float, from_none], obj.get("distance_donor_water"))
        donor_angle = from_union([from_float, from_none], obj.get("donor_angle"))
        water_angle = from_union([from_float, from_none], obj.get("water_angle"))
        protein_is_donor = from_union([from_bool, from_none], obj.get("protein_is_donor"))
        donor_atom_id = from_union([from_int, from_none], obj.get("donor_atom_id"))
        donor_group_type = from_union([from_str, from_none], obj.get("donor_group_type"))
        acceptor_atom_id = from_union([from_int, from_none], obj.get("acceptor_atom_id"))
        acceptor_group_type = from_union([from_str, from_none], obj.get("acceptor_group_type"))
        water_atom_id = from_union([from_str, from_none], obj.get("water_atom_id"))
        water_atom_position = from_union([lambda x: from_list(from_float, x), from_none],
                                         obj.get("water_atom_position"))
        ligand_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                    obj.get("ligand_representative_position"))
        protein_representative_position = from_union([lambda x: from_list(from_float, x), from_none],
                                                     obj.get("protein_representative_position"))
        return WaterBridge(binding_site_member_id, ligand_member_id, distance_acceptor_water, distance_donor_water,
                           donor_angle, water_angle, protein_is_donor, donor_atom_id, donor_group_type,
                           acceptor_atom_id, acceptor_group_type, water_atom_id, water_atom_position,
                           ligand_representative_position, protein_representative_position)

    def to_dict(self) -> dict:
        result: dict = {}
        result["binding_site_member_id"] = from_union([from_int, from_none], self.binding_site_member_id)
        result["ligand_member_id"] = from_union([from_int, from_none], self.ligand_member_id)
        result["distance_acceptor_water"] = from_union([to_float, from_none], self.distance_acceptor_water)
        result["distance_donor_water"] = from_union([to_float, from_none], self.distance_donor_water)
        result["donor_angle"] = from_union([to_float, from_none], self.donor_angle)
        result["water_angle"] = from_union([to_float, from_none], self.water_angle)
        result["protein_is_donor"] = from_union([from_bool, from_none], self.protein_is_donor)
        result["donor_atom_id"] = from_union([from_int, from_none], self.donor_atom_id)
        result["donor_group_type"] = from_union([from_str, from_none], self.donor_group_type)
        result["acceptor_atom_id"] = from_union([from_int, from_none], self.acceptor_atom_id)
        result["acceptor_group_type"] = from_union([from_str, from_none], self.acceptor_group_type)
        result["water_atom_id"] = from_union([from_str, from_none], self.water_atom_id)
        result["water_atom_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                   self.water_atom_position)
        result["ligand_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                              self.ligand_representative_position)
        result["protein_representative_position"] = from_union([lambda x: from_list(to_float, x), from_none],
                                                               self.protein_representative_position)
        return result


class Interactions:
    hydrogen_bonds: Optional[List[HydrogenBond]]
    hydrophobic_contacts: Optional[List[HydrophobicContact]]
    water_bridges: Optional[List[WaterBridge]]
    salt_bridges: Optional[List[SaltBridge]]
    pi_stackings: Optional[List[Pi]]
    pi_cations: Optional[List[Pi]]
    halogen_bonds: Optional[List[HalogenBond]]
    metal_complexes: Optional[List[MetalComplex]]

    def __init__(self, hydrogen_bonds: Optional[List[HydrogenBond]],
                 hydrophobic_contacts: Optional[List[HydrophobicContact]], water_bridges: Optional[List[WaterBridge]],
                 salt_bridges: Optional[List[SaltBridge]], pi_stackings: Optional[List[Pi]],
                 pi_cations: Optional[List[Pi]], halogen_bonds: Optional[List[HalogenBond]],
                 metal_complexes: Optional[List[MetalComplex]]) -> None:
        self.hydrogen_bonds = hydrogen_bonds
        self.hydrophobic_contacts = hydrophobic_contacts
        self.water_bridges = water_bridges
        self.salt_bridges = salt_bridges
        self.pi_stackings = pi_stackings
        self.pi_cations = pi_cations
        self.halogen_bonds = halogen_bonds
        self.metal_complexes = metal_complexes

    @staticmethod
    def from_dict(obj: Any) -> 'Interactions':
        assert isinstance(obj, dict)
        hydrogen_bonds = from_union([lambda x: from_list(HydrogenBond.from_dict, x), from_none],
                                    obj.get("hydrogen_bonds"))
        hydrophobic_contacts = from_union([lambda x: from_list(HydrophobicContact.from_dict, x), from_none],
                                          obj.get("hydrophobic_contacts"))
        water_bridges = from_union([lambda x: from_list(WaterBridge.from_dict, x), from_none], obj.get("water_bridges"))
        salt_bridges = from_union([lambda x: from_list(SaltBridge.from_dict, x), from_none], obj.get("salt_bridges"))
        pi_stackings = from_union([lambda x: from_list(Pi.from_dict, x), from_none], obj.get("pi_stackings"))
        pi_cations = from_union([lambda x: from_list(Pi.from_dict, x), from_none], obj.get("pi_cations"))
        halogen_bonds = from_union([lambda x: from_list(HalogenBond.from_dict, x), from_none], obj.get("halogen_bonds"))
        metal_complexes = from_union([lambda x: from_list(MetalComplex.from_dict, x), from_none],
                                     obj.get("metal_complexes"))
        return Interactions(hydrogen_bonds, hydrophobic_contacts, water_bridges, salt_bridges, pi_stackings, pi_cations,
                            halogen_bonds, metal_complexes)

    def to_dict(self) -> dict:
        result: dict = {}
        result["hydrogen_bonds"] = from_union([lambda x: from_list(lambda x: to_class(HydrogenBond, x), x), from_none],
                                              self.hydrogen_bonds)
        result["hydrophobic_contacts"] = from_union(
            [lambda x: from_list(lambda x: to_class(HydrophobicContact, x), x), from_none], self.hydrophobic_contacts)
        result["water_bridges"] = from_union([lambda x: from_list(lambda x: to_class(WaterBridge, x), x), from_none],
                                             self.water_bridges)
        result["salt_bridges"] = from_union([lambda x: from_list(lambda x: to_class(SaltBridge, x), x), from_none],
                                            self.salt_bridges)
        result["pi_stackings"] = from_union([lambda x: from_list(lambda x: to_class(Pi, x), x), from_none],
                                            self.pi_stackings)
        result["pi_cations"] = from_union([lambda x: from_list(lambda x: to_class(Pi, x), x), from_none],
                                          self.pi_cations)
        result["halogen_bonds"] = from_union([lambda x: from_list(lambda x: to_class(HalogenBond, x), x), from_none],
                                             self.halogen_bonds)
        result["metal_complexes"] = from_union([lambda x: from_list(lambda x: to_class(MetalComplex, x), x), from_none],
                                               self.metal_complexes)
        return result


class LigandMember:
    ligand_member_id: Optional[int]
    het_id: Optional[str]
    chain_id: Optional[str]
    residue_number: Optional[int]
    insertion_code: Optional[str]

    def __init__(self, ligand_member_id: Optional[int], het_id: Optional[str], chain_id: Optional[str],
                 residue_number: Optional[int], insertion_code: Optional[str]) -> None:
        self.ligand_member_id = ligand_member_id
        self.het_id = het_id
        self.chain_id = chain_id
        self.residue_number = residue_number
        self.insertion_code = insertion_code

    @staticmethod
    def from_dict(obj: Any) -> 'LigandMember':
        assert isinstance(obj, dict)
        ligand_member_id = from_union([from_int, from_none], obj.get("ligand_member_id"))
        het_id = from_union([from_str, from_none], obj.get("het_id"))
        chain_id = from_union([from_str, from_none], obj.get("chain_id"))
        residue_number = from_union([from_int, from_none], obj.get("residue_number"))
        insertion_code = from_union([from_str, from_none], obj.get("insertion_code"))
        return LigandMember(ligand_member_id, het_id, chain_id, residue_number, insertion_code)

    def to_dict(self) -> dict:
        result: dict = {}
        result["ligand_member_id"] = from_union([from_int, from_none], self.ligand_member_id)
        result["het_id"] = from_union([from_str, from_none], self.het_id)
        result["chain_id"] = from_union([from_str, from_none], self.chain_id)
        result["residue_number"] = from_union([from_int, from_none], self.residue_number)
        result["insertion_code"] = from_union([from_str, from_none], self.insertion_code)
        return result


class BindingSite:
    smiles: Optional[str]
    inchi: Optional[str]
    inchi_key: Optional[str]
    ligand_properties: Optional[Dict[str, float]]
    ligand_members: Optional[List[LigandMember]]
    binding_site_members: Optional[List[BindingSiteMember]]
    interactions: Optional[Interactions]

    def __init__(self, smiles: Optional[str], inchi: Optional[str], inchi_key: Optional[str],
                 ligand_properties: Optional[Dict[str, float]], ligand_members: Optional[List[LigandMember]],
                 binding_site_members: Optional[List[BindingSiteMember]], interactions: Optional[Interactions]) -> None:
        self.smiles = smiles
        self.inchi = inchi
        self.inchi_key = inchi_key
        self.ligand_properties = ligand_properties
        self.ligand_members = ligand_members
        self.binding_site_members = binding_site_members
        self.interactions = interactions

    @staticmethod
    def from_dict(obj: Any) -> 'BindingSite':
        assert isinstance(obj, dict)
        smiles = from_union([from_str, from_none], obj.get("smiles"))
        inchi = from_union([from_str, from_none], obj.get("inchi"))
        inchi_key = from_union([from_str, from_none], obj.get("inchi_key"))
        ligand_properties = from_union([lambda x: from_dict(from_float, x), from_none], obj.get("ligand_properties"))
        ligand_members = from_union([lambda x: from_list(LigandMember.from_dict, x), from_none],
                                    obj.get("ligand_members"))
        binding_site_members = from_union([lambda x: from_list(BindingSiteMember.from_dict, x), from_none],
                                          obj.get("binding_site_members"))
        interactions = from_union([Interactions.from_dict, from_none], obj.get("interactions"))
        return BindingSite(smiles, inchi, inchi_key, ligand_properties, ligand_members, binding_site_members,
                           interactions)

    def to_dict(self) -> dict:
        result: dict = {}
        result["smiles"] = from_union([from_str, from_none], self.smiles)
        result["inchi"] = from_union([from_str, from_none], self.inchi)
        result["inchi_key"] = from_union([from_str, from_none], self.inchi_key)
        result["ligand_properties"] = from_union([lambda x: from_dict(to_float, x), from_none], self.ligand_properties)
        result["ligand_members"] = from_union([lambda x: from_list(lambda x: to_class(LigandMember, x), x), from_none],
                                              self.ligand_members)
        result["binding_site_members"] = from_union(
            [lambda x: from_list(lambda x: to_class(BindingSiteMember, x), x), from_none], self.binding_site_members)
        result["interactions"] = from_union([lambda x: to_class(Interactions, x), from_none], self.interactions)
        return result


class PLIPAnalysis:
    plip_version: Optional[str]
    time_stamp: Optional[datetime]
    input_file: Optional[str]
    input_file_md5: Optional[str]
    excluded_ligands: Optional[List[str]]
    binding_sites: Optional[List[BindingSite]]

    def __init__(self, plip_version: Optional[str], time_stamp: Optional[datetime], input_file: Optional[str],
                 input_file_md5: Optional[str], excluded_ligands: Optional[List[str]],
                 binding_sites: Optional[List[BindingSite]]) -> None:
        self.plip_version = plip_version
        self.time_stamp = time_stamp
        self.input_file = input_file
        self.input_file_md5 = input_file_md5
        self.excluded_ligands = excluded_ligands
        self.binding_sites = binding_sites

    @staticmethod
    def from_dict(obj: Any) -> 'PLIPAnalysis':
        assert isinstance(obj, dict)
        plip_version = from_union([from_str, from_none], obj.get("plip_version"))
        time_stamp = from_union([from_datetime, from_none], obj.get("time_stamp"))
        input_file = from_union([from_str, from_none], obj.get("input_file"))
        input_file_md5 = from_union([from_str, from_none], obj.get("input_file_md5"))
        excluded_ligands = from_union([lambda x: from_list(from_str, x), from_none], obj.get("excluded_ligands"))
        binding_sites = from_union([lambda x: from_list(BindingSite.from_dict, x), from_none], obj.get("binding_sites"))
        return PLIPAnalysis(plip_version, time_stamp, input_file, input_file_md5, excluded_ligands, binding_sites)

    def to_dict(self) -> dict:
        result: dict = {}
        result["plip_version"] = from_union([from_str, from_none], self.plip_version)
        result["time_stamp"] = from_union([lambda x: x.isoformat(), from_none], self.time_stamp)
        result["input_file"] = from_union([from_str, from_none], self.input_file)
        result["input_file_md5"] = from_union([from_str, from_none], self.input_file_md5)
        result["excluded_ligands"] = from_union([lambda x: from_list(from_str, x), from_none], self.excluded_ligands)
        result["binding_sites"] = from_union([lambda x: from_list(lambda x: to_class(BindingSite, x), x), from_none],
                                             self.binding_sites)
        return result


def plip_analysis_from_dict(s: Any) -> PLIPAnalysis:
    return PLIPAnalysis.from_dict(s)


def plip_analysis_to_dict(x: PLIPAnalysis) -> Any:
    return to_class(PLIPAnalysis, x)


def plip_analysis_from_pdb_complex(pdb_complex: PDBComplex) -> PLIPAnalysis:
    # compute hash of input file
    pdb_input = pdb_complex.sourcefiles['pdbcomplex.original']
    hash_md5 = hashlib.md5()
    with open(pdb_input, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)

    plip_analysis = PLIPAnalysis(config.__version__, datetime.now(), pdb_input, hash_md5.hexdigest(), None, None)

    # create ligand members
    for ligand in pdb_complex.ligands:
        ligand_members = []
        for i, member in enumerate(ligand.members):
            ligand_member = LigandMember(i, member[0], member[1], member[2], '')  # todo insertion code support
            ligand_members.append(ligand_member)
