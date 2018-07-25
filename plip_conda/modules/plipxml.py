"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
plipxml.py - Read in PLIP XML files for further analysis.
"""

# Python standard library
from collections import defaultdict
from lxml import etree
import itertools
from itertools import groupby
import sys

try: # Python 3
    from urllib.request import urlopen
except ImportError: # Fallback Python 2.x
    from urllib2 import urlopen


class XMLStorage:
    """Generic class for storing XML data from PLIP XML files."""

    def getdata(self, tree, location, force_string=False):
        """Gets XML data from a specific element and handles types."""
        found = tree.xpath('%s/text()' % location)
        if not found:
            return None
        else:
            data = found[0]
        if force_string:
            return data
        if data == 'True':
            return True
        elif data == 'False':
            return False
        else:
            try:
                return int(data)
            except ValueError:
                try:
                    return float(data)
                except ValueError:
                    # It's a string
                    return data

    def getcoordinates(self, tree, location):
        """Gets coordinates from a specific element in PLIP XML"""
        return tuple(float(x) for x in tree.xpath('.//%s/*/text()' % location))


class Interaction(XMLStorage):
    """Stores information on a specific interaction type"""

    def __init__(self, interaction_part):
        self.id = interaction_part.get('id')
        self.resnr = self.getdata(interaction_part, 'resnr')
        self.restype = self.getdata(interaction_part, 'restype', force_string=True)
        self.reschain = self.getdata(interaction_part, 'reschain', force_string=True)
        self.resnr_lig = self.getdata(interaction_part, 'resnr_lig')
        self.restype_lig = self.getdata(interaction_part, 'restype_lig', force_string=True)
        self.reschain_lig = self.getdata(interaction_part, 'reschain_lig', force_string=True)
        self.ligcoo = self.getcoordinates(interaction_part, 'ligcoo')
        self.protcoo = self.getcoordinates(interaction_part, 'protcoo')



class HydrophobicInteraction(Interaction):
    """Stores information on a hydrophobic interaction"""

    def __init__(self, hydrophobic_part):
        Interaction.__init__(self, hydrophobic_part)
        self.dist = self.getdata(hydrophobic_part, 'dist')
        self.ligcarbonidx = self.getdata(hydrophobic_part, 'ligcarbonidx')
        self.protcarbonidx = self.getdata(hydrophobic_part, 'protcarbonidx')


class HydrogenBond(Interaction):
    """Stores information on a hydrogen bond interaction"""

    def __init__(self, hbond_part):
        Interaction.__init__(self, hbond_part)
        self.sidechain = self.getdata(hbond_part, 'sidechain')
        self.dist_h_a = self.getdata(hbond_part, 'dist_h-a')
        self.dist_d_a = self.getdata(hbond_part, 'dist_d-a')
        self.dist = self.dist_d_a

        self.don_angle = self.getdata(hbond_part, 'don_angle')
        self.protisdon = self.getdata(hbond_part, 'protisdon')
        self.donoridx = self.getdata(hbond_part, 'donoridx')
        self.acceptoridx = self.getdata(hbond_part, 'acceptoridx')
        self.donortype = self.getdata(hbond_part, 'donortype', force_string=True)
        self.acceptortype = self.getdata(hbond_part, 'acceptortype', force_string=True)


class WaterBridge(Interaction):
    """Stores information on a water bridge interaction"""

    def __init__(self, wbridge_part):
        Interaction.__init__(self, wbridge_part)
        self.dist_a_w = self.getdata(wbridge_part, 'dist_a-w')
        self.dist_d_w = self.getdata(wbridge_part, 'dist_d-w')
        self.don_angle = self.getdata(wbridge_part, 'don_angle')
        self.water_angle = self.getdata(wbridge_part, 'water_angle')
        self.protisdon = self.getdata(wbridge_part, 'protisdon')
        self.dist = self.dist_a_w if self.protisdon else self.dist_d_w

        self.donor_idx = self.getdata(wbridge_part, 'donor_idx')
        self.acceptor_idx = self.getdata(wbridge_part, 'acceptor_idx')
        self.donortype = self.getdata(wbridge_part, 'donortype', force_string=True)
        self.acceptortype = self.getdata(wbridge_part, 'acceptortype', force_string=True)
        self.water_idx = self.getdata(wbridge_part, 'water_idx')
        self.watercoo = self.getcoordinates(wbridge_part, 'watercoo')


class SaltBridge(Interaction):
    """Stores information on a salt bridge interaction"""

    def __init__(self, sbridge_part):
        Interaction.__init__(self, sbridge_part)
        self.dist =  self.getdata(sbridge_part, 'dist')
        self.protispos = self.getdata(sbridge_part, 'protispos')
        self.lig_group = self.getdata(sbridge_part, 'lig_group', force_string=True)
        self.lig_idx_list = [int(tagpart.text) for tagpart in
                             sbridge_part.xpath('lig_idx_list/idx')]


class PiStacking(Interaction):
    """Stores information on a pi stacking interaction"""

    def __init__(self, pistack_part):
        Interaction.__init__(self, pistack_part)
        self.centdist = self.getdata(pistack_part, 'centdist')
        self.dist = self.centdist
        self.angle = self.getdata(pistack_part, 'angle')
        self.offset = self.getdata(pistack_part, 'offset')
        self.type = self.getdata(pistack_part, 'type')
        self.lig_idx_list = [int(tagpart.text) for tagpart in
                             pistack_part.xpath('lig_idx_list/idx')]


class PiCation(Interaction):
    """Stores information on a pi cation interaction"""

    def __init__(self, pication_part):
        Interaction.__init__(self, pication_part)
        self.dist = self.getdata(pication_part, 'dist')
        self.offset = self.getdata(pication_part, 'offset')
        self.protcharged = self.getdata(pication_part, 'protcharged')
        self.lig_group = self.getdata(pication_part, 'lig_group')
        self.lig_idx_list = [int(tag.text) for tag in pication_part.xpath('.//lig_idx_list/idx')]


class HalogenBond(Interaction):
    """Stores information on a halogen bond interaction"""

    def __init__(self, halogen_part):
        Interaction.__init__(self, halogen_part)
        self.dist = self.getdata(halogen_part, 'dist')
        self.don_angle = self.getdata(halogen_part, 'don_angle')
        self.acc_angle = self.getdata(halogen_part, 'acc_angle')
        self.donortype = self.getdata(halogen_part, 'donortype', force_string=True)
        self.acceptortype = self.getdata(halogen_part, 'acceptortype', force_string=True)
        self.don_idx = self.getdata(halogen_part, 'don_idx')
        self.acc_idx = self.getdata(halogen_part, 'acc_idx')
        self.sidechain = self.getdata(halogen_part, 'sidechain')


class MetalComplex(Interaction):
    """Stores information on a metal complexe interaction"""

    def __init__(self, metalcomplex_part):
        Interaction.__init__(self, metalcomplex_part)
        self.metal_idx = self.getdata(metalcomplex_part, 'metal_idx')
        self.metal_type = self.getdata(metalcomplex_part, 'metal_type', force_string=True)
        self.target_idx = self.getdata(metalcomplex_part, 'target_idx')
        self.target_type = self.getdata(metalcomplex_part, 'target_type', force_string=True)
        self.coordination = self.getdata(metalcomplex_part, 'coordination')
        self.dist = self.getdata(metalcomplex_part, 'dist')
        self.location = self.getdata(metalcomplex_part, 'location', force_string=True)
        self.rms = self.getdata(metalcomplex_part, 'rms')
        self.geometry = self.getdata(metalcomplex_part, 'geometry', force_string=True)
        self.complexnum = self.getdata(metalcomplex_part, 'complexnum')
        self.targetcoo = self.getcoordinates(metalcomplex_part, 'targetcoo')
        self.metalcoo = self.getcoordinates(metalcomplex_part, 'metalcoo')

class BSite(XMLStorage):
    """Stores all information about an specific binding site."""

    def __init__(self, bindingsite, pdbid):
        self.bindingsite = bindingsite
        self.pdbid = pdbid
        self.bsid = ":".join(bindingsite.xpath('identifiers/*/text()')[2:5])
        self.uniqueid = ":".join([self.pdbid, self.bsid])
        self.hetid = self.getdata(bindingsite, 'identifiers/hetid', force_string=True)
        self.longname = self.getdata(bindingsite, 'identifiers/longname', force_string=True)
        self.ligtype = self.getdata(bindingsite, 'identifiers/ligtype', force_string=True)
        self.smiles = self.getdata(bindingsite, 'identifiers/smiles', force_string=True)
        self.inchikey = self.getdata(bindingsite, 'identifiers/inchikey', force_string=True)
        self.position = self.getdata(bindingsite, 'identifiers/position')
        self.chain = self.getdata(bindingsite, 'identifiers/chain', force_string=True)

        # Information on binding site members
        self.members = []
        for member in bindingsite.xpath('identifiers/members/member'):
            self.members += member.xpath('text()')

        self.composite = self.getdata(bindingsite, 'identifiers/composite')

        # Ligand Properties
        self.heavy_atoms = self.getdata(bindingsite, 'lig_properties/num_heavy_atoms')
        self.hbd = self.getdata(bindingsite, 'lig_properties/num_hbd')
        self.unpaired_hbd = self.getdata(bindingsite, 'lig_properties/num_unpaired_hbd')
        self.hba = self.getdata(bindingsite, 'lig_properties/num_hba')
        self.unpaired_hba = self.getdata(bindingsite, 'lig_properties/num_unpaired_hba')
        self.hal = self.getdata(bindingsite, 'lig_properties/num_hal')
        self.unpaired_hal = self.getdata(bindingsite, 'lig_properties/num_unpaired_hal')
        self.molweight = self.getdata(bindingsite, 'lig_properties/molweight')
        self.logp = self.getdata(bindingsite, 'lig_properties/logp')
        self.rotatable_bonds = self.getdata(bindingsite, 'lig_properties/num_rotatable_bonds')
        self.rings = self.getdata(bindingsite, 'lig_properties/num_aromatic_rings')


        # Binding Site residues
        self.bs_res = []
        for tagpart in bindingsite.xpath('bs_residues/bs_residue'):
            resnumber, reschain = tagpart.text[:-1], tagpart.text[-1]
            aa, contact, min_dist = tagpart.get('aa'), tagpart.get('contact'), tagpart.get('min_dist')
            new_bs_res = {'resnr': int(resnumber), 'reschain': reschain, 'aa': aa, 'contact': True if contact == 'True' else False, 'min_dist': float(min_dist)}
            self.bs_res.append(new_bs_res)

        # Interacting chains
        self.interacting_chains = []
        for chain in bindingsite.xpath('interacting_chains/interacting_chain'):
            self.interacting_chains += chain.xpath('text()')

        # Interactions
        interactions = bindingsite.xpath('interactions')[0]
        self.hydrophobics = [HydrophobicInteraction(x) for x in
                             interactions.xpath('hydrophobic_interactions/hydrophobic_interaction')]
        self.hbonds = [HydrogenBond(x) for x in interactions.xpath('hydrogen_bonds/hydrogen_bond')]
        self.wbridges = [WaterBridge(x) for x in interactions.xpath('water_bridges/water_bridge')]
        self.sbridges = [SaltBridge(x) for x in interactions.xpath('salt_bridges/salt_bridge')]
        self.pi_stacks = [PiStacking(x) for x in interactions.xpath('pi_stacks/pi_stack')]
        self.pi_cations = [PiCation(x) for x in interactions.xpath('pi_cation_interactions/pi_cation_interaction')]
        self.halogens = [HalogenBond(x) for x in interactions.xpath('halogen_bonds/halogen_bond')]
        self.metal_complexes = [MetalComplex(x) for x in interactions.xpath('metal_complexes/metal_complex')]
        self.num_contacts = len(self.hydrophobics) + len(self.hbonds) + len(self.wbridges) + len(self.sbridges) + len(self.pi_stacks) + len(self.pi_cations) + len(self.halogens) + len(self.metal_complexes)
        self.has_interactions = self.num_contacts > 0

        self.get_atom_mapping()
        self.counts = self.get_counts()

    def get_atom_mapping(self):
        """Parses the ligand atom mapping."""
        # Atom mappings
        smiles_to_pdb_mapping = self.bindingsite.xpath('mappings/smiles_to_pdb/text()')
        if smiles_to_pdb_mapping == []:
            self.mappings = {'smiles_to_pdb': None, 'pdb_to_smiles': None}
        else:
            smiles_to_pdb_mapping = {int(y[0]): int(y[1]) for y in [x.split(':') for x in smiles_to_pdb_mapping[0].split(',')]}
            self.mappings = {'smiles_to_pdb': smiles_to_pdb_mapping}
            self.mappings['pdb_to_smiles'] = {v: k for k, v in self.mappings['smiles_to_pdb'].items()}

    def get_counts(self):
        """counts the interaction types and backbone hydrogen bonding in a binding site"""

        hbondsback = len([hb for hb in self.hbonds if not hb.sidechain])
        counts = {'hydrophobics': len(self.hydrophobics), 'hbonds': len(self.hbonds),
                      'wbridges': len(self.wbridges), 'sbridges': len(self.sbridges), 'pistacks': len(self.pi_stacks),
                      'pications': len(self.pi_cations), 'halogens': len(self.halogens), 'metal': len(self.metal_complexes),
                      'hbond_back': hbondsback, 'hbond_nonback': (len(self.hbonds) - hbondsback)}
        counts['total'] = counts['hydrophobics'] + counts['hbonds'] + counts['wbridges'] + counts['sbridges'] + counts['pistacks'] + counts['pications'] + counts['halogens'] + counts['metal']
        return counts

class PLIPXML(XMLStorage):
    """Parses and stores all information from a PLIP XML file."""
    def __init__(self, xmlfile):
        self.load_data(xmlfile)

        # Parse general information
        self.version = self.getdata(self.doc, '/report/plipversion/')
        self.pdbid = self.getdata(self.doc, '/report/pdbid', force_string=True)
        self.filetype = self.getdata(self.doc, '/report/filetype')
        self.fixed = self.getdata(self.doc, '/report/pdbfixes/')
        self.filename = self.getdata(self.doc, '/report/filename')
        self.excluded = self.doc.xpath('/report/excluded_ligands/excluded_ligand/text()')

        # Parse binding site information
        self.bsites = {BSite(bs, self.pdbid).bsid: BSite(bs, self.pdbid) for bs in self.doc.xpath('//bindingsite')}
        self.num_bsites = len(self.bsites)


    def load_data(self, xmlfile):
        """Loads/parses an XML file and saves it as a tree if successful."""
        self.doc = etree.parse(xmlfile)

class PLIPXMLREST(PLIPXML):
    """Parses and stores all from a PLIP XML file from the PLIP REST service"""
    def __init__(self, pdbid):
        PLIPXML.__init__(self, pdbid)

    def load_data(self, pdbid):
        """Loads and parses an XML resource and saves it as a tree if successful"""
        #TODO Implement error handling
        f = urlopen("http://projects.biotec.tu-dresden.de/plip-rest/pdb/%s?format=xml" % pdbid.lower())
        self.doc = etree.parse(f)
