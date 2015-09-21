"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
report.py - Write PLIP results to output files.
Copyright 2014-2015 Sebastian Salentin

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""


# Python Standard Library
import time
from operator import itemgetter

# External libraries
import lxml.etree as et


class TextOutput:
    """Gather report data and generate reports for one binding site in different formats"""
    def __init__(self, plcomplex):

        ################
        # GENERAL DATA #
        ################

        self.complex = plcomplex
        self.ligand = self.complex.ligand
        self.bindingsite = self.complex.bindingsite
        self.output_path = self.complex.output_path
        self.bsid = ':'.join([self.ligand.hetid, self.ligand.chain, str(self.ligand.position)])
        self.longname = self.ligand.longname
        self.ligtype = self.ligand.type
        self.bs_res = self.bindingsite.bs_res
        self.min_dist = self.bindingsite.min_dist
        self.bs_res_interacting = self.complex.interacting_res
        self.pdbid = self.complex.pdbid.upper()
        self.lig_members = self.complex.lig_members
        self.interacting_chains = self.complex.interacting_chains

        ############################
        # HYDROPHOBIC INTERACTIONS #
        ############################

        self.hydrophobic_features = ('RESNR', 'RESTYPE', 'RESCHAIN', 'DIST', 'LIGCARBONIDX', 'PROTCARBONIDX', 'LIGCOO',
                                     'PROTCOO')
        self.hydrophobic_info = []
        for hydroph in self.complex.hydrophobic_contacts:
            self.hydrophobic_info.append((hydroph.resnr, hydroph.restype, hydroph.reschain, '%.2f' % hydroph.distance,
                                          hydroph.ligatom_orig_idx, hydroph.bsatom_orig_idx, hydroph.ligatom.coords,
                                          hydroph.bsatom.coords))

        ##################
        # HYDROGEN BONDS #
        ##################

        self.hbond_features = ('RESNR', 'RESTYPE', 'RESCHAIN', 'SIDECHAIN', 'DIST_H-A', 'DIST_D-A', 'DON_ANGLE',
                               'PROTISDON', 'DONORIDX', 'DONORTYPE', 'ACCEPTORIDX', 'ACCEPTORTYPE', 'LIGCOO', 'PROTCOO')
        self.hbond_info = []
        for hbond in self.complex.hbonds_pdon + self.complex.hbonds_ldon:
            ligatom, protatom = (hbond.a, hbond.d) if hbond.protisdon else (hbond.d, hbond.a)
            self.hbond_info.append((hbond.resnr, hbond.restype, hbond.reschain, hbond.sidechain,
                                    '%.2f' % hbond.distance_ah, '%.2f' % hbond.distance_ad, '%.2f' % hbond.angle,
                                    hbond.protisdon, hbond.d_orig_idx, hbond.dtype, hbond.a_orig_idx, hbond.atype,
                                    ligatom.coords, protatom.coords))

        #################
        # WATER-BRIDGES #
        #################

        self.waterbridge_features = ('RESNR', 'RESTYPE', 'RESCHAIN', 'DIST_A-W', 'DIST_D-W', 'DON_ANGLE', 'WATER_ANGLE',
                                     'PROTISDON', 'DONOR_IDX', 'DONORTYPE', 'ACCEPTOR_IDX', 'ACCEPTORTYPE', 'WATER_IDX',
                                     'LIGCOO', 'PROTCOO', 'WATERCOO')
        # The coordinate format is an exception here, since the interaction is not only between ligand and protein
        self.waterbridge_info = []
        for wbridge in self.complex.water_bridges:
            lig, prot = (wbridge.a, wbridge.d) if wbridge.protisdon else (wbridge.d, wbridge.a)
            self.waterbridge_info.append((wbridge.resnr, wbridge.restype, wbridge.reschain,
                                          '%.2f' % wbridge.distance_aw, '%.2f' % wbridge.distance_dw,
                                          '%.2f' % wbridge.d_angle, '%.2f' % wbridge.w_angle, wbridge.protisdon,
                                          wbridge.d_orig_idx, wbridge.dtype, wbridge.a_orig_idx, wbridge.atype,
                                          wbridge.water_orig_idx, lig.coords, prot.coords, wbridge.water.coords))

        ################
        # SALT BRIDGES #
        ################

        self.saltbridge_features = ('RESNR', 'RESTYPE', 'RESCHAIN', 'DIST', 'PROTISPOS', 'LIG_GROUP', 'LIG_IDX_LIST',
                                    'LIGCOO', 'PROTCOO')
        self.saltbridge_info = []
        for sb in self.complex.saltbridge_lneg + self.complex.saltbridge_pneg:
            if sb.protispos:
                group, ids = sb.negative.fgroup, [str(x) for x in sb.negative.atoms_orig_idx]
                self.saltbridge_info.append((sb.resnr, sb.restype, sb.reschain, '%.2f' % sb.distance, sb.protispos,
                                             group.capitalize(), ",".join(ids),
                                             tuple(sb.negative.center), tuple(sb.positive.center)))
            else:
                group, ids = sb.positive.fgroup, [str(x) for x in sb.positive.atoms_orig_idx]
                self.saltbridge_info.append((sb.resnr, sb.restype, sb.reschain, '%.2f' % sb.distance, sb.protispos,
                                             group.capitalize(), ",".join(ids),
                                             tuple(sb.positive.center), tuple(sb.negative.center)))

        ###############
        # PI-STACKING #
        ###############

        self.pistacking_features = ('RESNR', 'RESTYPE', 'RESCHAIN', 'CENTDIST', 'ANGLE', 'OFFSET', 'TYPE',
                                    'LIG_IDX_LIST', 'LIGCOO', 'PROTCOO')
        self.pistacking_info = []
        for stack in self.complex.pistacking:
            ids = [str(x) for x in stack.ligandring.atoms_orig_idx]
            self.pistacking_info.append((stack.resnr, stack.restype, stack.reschain, '%.2f' % stack.distance,
                                         '%.2f' % stack.angle, '%.2f' % stack.offset, stack.type, ",".join(ids),
                                         tuple(stack.ligandring.center), tuple(stack.proteinring.center)))

        ##########################
        # PI-CATION INTERACTIONS #
        ##########################

        self.pication_features = ('RESNR', 'RESTYPE', 'RESCHAIN', 'DIST', 'OFFSET', 'PROTCHARGED', 'LIG_GROUP',
                                  'LIG_IDX_LIST', 'LIGCOO', 'PROTCOO')
        self.pication_info = []
        for picat in self.complex.pication_laro + self.complex.pication_paro:
            if picat.protcharged:
                ids = [str(x) for x in picat.ring.atoms_orig_idx]
                group = 'Aromatic'
                self.pication_info.append((picat.resnr, picat.restype, picat.reschain, '%.2f' % picat.distance,
                                           '%.2f' % picat.offset, picat.protcharged, group, ",".join(ids),
                                           tuple(picat.ring.center), tuple(picat.charge.center)))
            else:
                ids = [str(x) for x in picat.charge.atoms_orig_idx]
                group = picat.charge.fgroup
                self.pication_info.append((picat.resnr, picat.restype, picat.reschain, '%.2f' % picat.distance,
                                           '%.2f' % picat.offset, picat.protcharged, group, ",".join(ids),
                                           tuple(picat.charge.center), tuple(picat.ring.center)))

        #################
        # HALOGEN BONDS #
        #################

        self.halogen_features = ('RESNR', 'RESTYPE', 'RESCHAIN', 'SIDECHAIN', 'DIST', 'DON_ANGLE', 'ACC_ANGLE',
                                 'DON_IDX', 'DONORTYPE', 'ACC_IDX', 'ACCEPTORTYPE', 'LIGCOO', 'PROTCOO')
        self.halogen_info = []
        for halogen in self.complex.halogen_bonds:
            self.halogen_info.append((halogen.resnr, halogen.restype, halogen.reschain, halogen.sidechain,
                                      '%.2f' % halogen.distance, '%.2f' % halogen.don_angle, '%.2f' % halogen.acc_angle,
                                      halogen.don_orig_idx, halogen.donortype,
                                      halogen.acc_orig_idx, halogen.acctype,
                                      halogen.acc.o.coords, halogen.don.x.coords))

        ###################
        # METAL COMPLEXES #
        ###################

        self.metal_features = ('RESNR', 'RESTYPE', 'RESCHAIN', 'METAL_IDX', 'METAL_TYPE', 'TARGET_IDX', 'TARGET_TYPE',
                               'COORDINATION', 'DIST', 'LOCATION', 'RMS', 'GEOMETRY', 'COMPLEXNUM', 'METALCOO',
                               'TARGETCOO')
        self.metal_info = []
        # Coordinate format here is non-standard since the interaction partner can be either ligand or protein
        for m in self.complex.metal_complexes:
            self.metal_info.append((m.resnr, m.restype, m.reschain, m.metal_orig_idx, m.metal_type,
                                    m.target_orig_idx, m.target_type, m.coordination_num, '%.2f' % m.distance,
                                    m.location, '%.2f' % m.rms, m.geometry, str(m.complexnum), m.metal.coords,
                                    m.target.atom.coords))

    def write_section(self, name, features, info, f):
        """Provides formatting for one section (e.g. hydrogen bonds)"""
        if not len(info) == 0:
            f.write('\n\n### %s ###\n' % name)
            f.write('%s\n' % '\t'.join(features))
            for line in info:
                f.write('%s\n' % '\t'.join(map(str, line)))

    def rst_table(self, array):
        """Given an array, the function formats and returns and table in rST format."""
        # Determine cell width for each column
        cell_dict = {}
        for i, row in enumerate(array):
            for j, val in enumerate(row):
                if j not in cell_dict:
                    cell_dict[j] = []
                cell_dict[j].append(val)
        for item in cell_dict:
            cell_dict[item] = max([len(x) for x in cell_dict[item]]) + 1  # Contains adapted width for each column

        # Format top line
        num_cols = len(array[0])
        form = '+'
        for col in range(num_cols):
            form += (cell_dict[col] + 1) * '-'
            form += '+'
        form += '\n'

        # Format values
        for i, row in enumerate(array):
            form += '| '
            for j, val in enumerate(row):
                cell_width = cell_dict[j]
                form += str(val) + (cell_width - len(val)) * ' ' + '| '
            form.rstrip()
            form += '\n'

            # Seperation lines
            form += '+'
            if i == 0:
                sign = '='
            else:
                sign = '-'
            for col in range(num_cols):
                form += (cell_dict[col] + 1) * sign
                form += '+'
            form += '\n'
        return form

    def generate_txt(self):
        """Generates an flat text report for a single binding site"""

        txt = []
        txt.append('%s (%s) - %s' % (self.bsid, self.longname, self.ligtype))
        for i, member in enumerate(self.lig_members[1:]):
            txt.append('  + %s' % ":".join(str(element) for element in member))
        txt.append("-" * len(self.bsid))
        txt.append("Interacting chain(s): %s\n" % ','.join([chain for chain in self.interacting_chains]))
        for section in [['Hydrophobic Interactions', self.hydrophobic_features, self.hydrophobic_info],
                        ['Hydrogen Bonds', self.hbond_features, self.hbond_info],
                        ['Water Bridges', self.waterbridge_features, self.waterbridge_info],
                        ['Salt Bridges', self.saltbridge_features, self.saltbridge_info],
                        ['pi-Stacking', self.pistacking_features, self.pistacking_info],
                        ['pi-Cation Interactions', self.pication_features, self.pication_info],
                        ['Halogen Bonds', self.halogen_features, self.halogen_info],
                        ['Metal Complexes', self.metal_features, self.metal_info]]:
            iname, features, interaction_information = section
            # Sort results first by res number, then by distance and finally ligand coordinates to get a unique order
            interaction_information = sorted(interaction_information, key=itemgetter(0, 2, -2))
            if not len(interaction_information) == 0:

                txt.append('\n**%s**' % iname)
                table = [features, ]
                for single_contact in interaction_information:
                    values = []
                    for x in single_contact:
                        if type(x) == str:
                            values.append(x)
                        elif type(x) == tuple and len(x) == 3:  # Coordinates
                            values.append("%.3f, %.3f, %.3f" % x)
                        else:
                            values.append(str(x))
                    table.append(values)
                txt.append(self.rst_table(table))
        txt.append('\n')
        return txt

    def generate_xml(self):
        """Generates an XML-formatted report for a single binding site"""
        report = et.Element('bindingsite')
        identifiers = et.SubElement(report, 'identifiers')
        longname = et.SubElement(identifiers, 'longname')
        ligtype = et.SubElement(identifiers, 'ligtype')
        hetid = et.SubElement(identifiers, 'hetid')
        chain = et.SubElement(identifiers, 'chain')
        position = et.SubElement(identifiers, 'position')
        composite = et.SubElement(identifiers, 'composite')
        members = et.SubElement(identifiers, 'members')
        smiles = et.SubElement(identifiers, 'smiles')

        # Ligand properties. Number of (unpaired) functional atoms and rings.
        lig_properties = et.SubElement(report, 'lig_properties')
        num_heavy_atoms = et.SubElement(lig_properties, 'num_heavy_atoms')
        num_hbd = et.SubElement(lig_properties, 'num_hbd')
        num_hbd.text = str(self.ligand.num_hbd)
        num_unpaired_hbd = et.SubElement(lig_properties, 'num_unpaired_hbd')
        num_unpaired_hbd.text = str(self.complex.num_unpaired_hbd)
        num_hba = et.SubElement(lig_properties, 'num_hba')
        num_hba.text = str(self.ligand.num_hba)
        num_unpaired_hba = et.SubElement(lig_properties, 'num_unpaired_hba')
        num_unpaired_hba.text = str(self.complex.num_unpaired_hba)
        num_hal = et.SubElement(lig_properties, 'num_hal')
        num_hal.text = str(self.ligand.num_hal)
        num_unpaired_hal = et.SubElement(lig_properties, 'num_unpaired_hal')
        num_unpaired_hal.text = str(self.complex.num_unpaired_hal)
        num_aromatic_rings = et.SubElement(lig_properties, 'num_aromatic_rings')
        num_aromatic_rings.text = str(self.ligand.num_rings)
        num_rot_bonds = et.SubElement(lig_properties, 'num_rotatable_bonds')
        num_rot_bonds.text = str(self.ligand.num_rot_bonds)
        molweight = et.SubElement(lig_properties, 'molweight')
        molweight.text = str(self.ligand.molweight)
        logp = et.SubElement(lig_properties, 'logp')
        logp.text = str(self.ligand.logp)

        ichains = et.SubElement(report, 'interacting_chains')
        bsresidues = et.SubElement(report, 'bs_residues')
        for i, ichain in enumerate(self.interacting_chains):
            c = et.SubElement(ichains, 'interacting_chain', id=str(i + 1))
            c.text = ichain
        for i, bsres in enumerate(self.bs_res):
            contact = 'True' if bsres in self.bs_res_interacting else 'False'
            distance = '%.1f' % self.min_dist[bsres][0]
            aatype = self.min_dist[bsres][1]
            c = et.SubElement(bsresidues, 'bs_residue', id=str(i + 1), contact=contact, min_dist=distance, aa=aatype)
            c.text = bsres
        hetid.text, chain.text, position.text = self.ligand.hetid, self.ligand.chain, str(self.ligand.position)
        composite.text = 'True' if len(self.lig_members) > 1 else 'False'
        longname.text = self.longname
        ligtype.text = self.ligtype
        smiles.text = self.ligand.smiles
        num_heavy_atoms.text = str(self.ligand.heavy_atoms)  # Number of heavy atoms in ligand
        for i, member in enumerate(self.lig_members):
            bsid = ":".join(str(element) for element in member)
            m = et.SubElement(members, 'member', id=str(i + 1))
            m.text = bsid
        interactions = et.SubElement(report, 'interactions')

        def format_interactions(element_name, features, interaction_information):
            """Returns a formatted element with interaction information."""
            interaction = et.Element(element_name)
            # Sort results first by res number, then by distance and finally ligand coordinates to get a unique order
            interaction_information = sorted(interaction_information, key=itemgetter(0, 2, -2))
            for j, single_contact in enumerate(interaction_information):
                if not element_name == 'metal_complexes':
                    new_contact = et.SubElement(interaction, element_name[:-1], id=str(j + 1))
                else:  # Metal Complex[es]
                    new_contact = et.SubElement(interaction, element_name[:-2], id=str(j + 1))
                for i, feature in enumerate(single_contact):
                    # Just assign the value unless it's an atom list, use subelements in this case
                    if features[i] == 'LIG_IDX_LIST':
                        feat = et.SubElement(new_contact, features[i].lower())
                        for k, atm_idx in enumerate(feature.split(',')):
                            idx = et.SubElement(feat, 'idx', id=str(k + 1))
                            idx.text = str(atm_idx)
                    elif features[i].endswith('COO'):
                        feat = et.SubElement(new_contact, features[i].lower())
                        xc, yc, zc = feature
                        xcoo = et.SubElement(feat, 'x')
                        xcoo.text = '%.3f' % xc
                        ycoo = et.SubElement(feat, 'y')
                        ycoo.text = '%.3f' % yc
                        zcoo = et.SubElement(feat, 'z')
                        zcoo.text = '%.3f' % zc
                    else:
                        feat = et.SubElement(new_contact, features[i].lower())
                        feat.text = str(feature)
            return interaction

        interactions.append(format_interactions('hydrophobic_interactions', self.hydrophobic_features,
                                                self.hydrophobic_info))
        interactions.append(format_interactions('hydrogen_bonds', self.hbond_features, self.hbond_info))
        interactions.append(format_interactions('water_bridges', self.waterbridge_features, self.waterbridge_info))
        interactions.append(format_interactions('salt_bridges', self.saltbridge_features, self.saltbridge_info))
        interactions.append(format_interactions('pi_stacks', self.pistacking_features, self.pistacking_info))
        interactions.append(format_interactions('pi_cation_interactions', self.pication_features, self.pication_info))
        interactions.append(format_interactions('halogen_bonds', self.halogen_features, self.halogen_info))
        interactions.append(format_interactions('metal_complexes', self.metal_features, self.metal_info))
        return report
