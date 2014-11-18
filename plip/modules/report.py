"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
Module for generation of textual reports.
Copyright (C) 2014  Sebastian Salentin

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
"""


# Python Standard Library
import time

# Own modules
from supplemental import *

# External libraries
import lxml.etree as et


class TextOutput():
    """Create and save a report for the pliprofiler of a binding site."""
    def __init__(self, pli_class):
        ### GENERAL
        self.output_path = pli_class.output_path
        self.name = pli_class.name
        self.pdbid = pli_class.pdbid.upper()
        mapping = pli_class.idx_to_pdb
        lig_to_pdb = {key: mapping[pli_class.lig_to_pdb[key]] for key in pli_class.lig_to_pdb}  # Atom mapping ligand
        self.header = ['#PREDICTION OF NONCOVALENT INTERACTIONS FOR %s:%s' % (self.pdbid, self.name),
                       '#Created on %s' % time.strftime("%Y/%m/%d")]

        ### HYDROPHOBIC
        self.hydrophobic_features = ('RESNR', 'RESTYPE', 'DIST', 'LIGCARBONIDX', 'PROTCARBONIDX', 'LIGCOO',
                                     'PROTCOO')
        self.hydrophobic_info = []
        for hydroph in pli_class.hydrophobic_contacts:
            self.hydrophobic_info.append((hydroph.resnr, hydroph.restype, '%.2f' % hydroph.distance,
                                          lig_to_pdb[hydroph.ligatom.idx], mapping[hydroph.bsatom.idx]
                                          , hydroph.ligatom.coords, hydroph.bsatom.coords))

        ### HBONDS
        self.hbond_features = ('RESNR', 'RESTYPE', 'DIST_H-A', 'DIST_D-A', 'DON_ANGLE', 'PROTISDON', 'DONORIDX',
                               'ACCEPTORIDX', 'LIGCOO', 'PROTCOO')
        self.hbond_info = []
        for hbond in pli_class.hbonds_pdon+pli_class.hbonds_ldon:
            if hbond.protisdon:
                donidx, accidx = mapping[hbond.d.idx], lig_to_pdb[hbond.a.idx]
                self.hbond_info.append((hbond.resnr, hbond.restype, '%.2f' % hbond.distance_ah,
                                        '%.2f' % hbond.distance_ad, '%.2f' % hbond.angle, hbond.protisdon, donidx,
                                        accidx, hbond.a.coords, hbond.d.coords))
            else:
                donidx, accidx = lig_to_pdb[hbond.d.idx], mapping[hbond.a.idx]
                self.hbond_info.append((hbond.resnr, hbond.restype, '%.2f' % hbond.distance_ah,
                                        '%.2f' % hbond.distance_ad, '%.2f' % hbond.angle, hbond.protisdon, donidx,
                                        accidx, hbond.d.coords, hbond.a.coords))

        ### WATER-BRIDGED HYDROGEN BONDS
        self.waterbridge_features = ('RESNR', 'RESTYPE', 'DIST_A-W', 'DIST_D-W', 'DON_ANGLE', 'WATER_ANGLE',
                                     'PROTISDON', 'DONOR_IDX', 'ACCEPTOR_IDX', 'WATER_IDX')
        self.waterbridge_info = []
        for wbridge in pli_class.water_bridges:
            if wbridge.protisdon:
                donidx, accidx = mapping[wbridge.d.idx], lig_to_pdb[wbridge.a.idx]
            else:
                donidx, accidx = lig_to_pdb[wbridge.d.idx], mapping[wbridge.a.idx]
            self.waterbridge_info.append((wbridge.resnr, wbridge.restype, '%.2f' % wbridge.distance_aw,
                                          '%.2f' % wbridge.distance_dw, '%.2f' % wbridge.d_angle,
                                          '%.2f' % wbridge.w_angle, wbridge.protisdon, donidx, accidx,
                                          mapping[wbridge.water.idx]))

        ### SALTBRIDGES
        self.saltbridge_features = ('RESNR', 'RESTYPE', 'DIST', 'PROTISPOS', 'LIG_GROUP', 'LIG_IDX_LIST',
                                    'LIGCOO', 'PROTCOO')
        self.saltbridge_info = []
        for sb in pli_class.saltbridge_lneg+pli_class.saltbridge_pneg:
            if sb.protispos:
                group, ids = sb.negative.fgroup, [str(lig_to_pdb[x.idx]) for x in sb.negative.atoms]
                self.saltbridge_info.append((sb.resnr, sb.restype, '%.2f' % sb.distance, sb.protispos,
                                             group.capitalize(), ",".join(ids),
                                             tuple(sb.negative.center), tuple(sb.positive.center)))
            else:
                group, ids = sb.positive.fgroup, [str(lig_to_pdb[x.idx]) for x in sb.positive.atoms]
                self.saltbridge_info.append((sb.resnr, sb.restype, '%.2f' % sb.distance, sb.protispos,
                                             group.capitalize(), ",".join(ids),
                                             tuple(sb.positive.center), tuple(sb.negative.center)))

        ### PISTACKING
        self.pistacking_features = ('RESNR', 'RESTYPE', 'CENTDIST', 'ANGLE', 'OFFSET', 'TYPE', 'LIG_IDX_LIST',
                                    'LIGCOO', 'PROTCOO')
        self.pistacking_info = []
        for stack in pli_class.pistacking:
            ids = [str(lig_to_pdb[x.idx]) for x in stack.ligandring.atoms]
            self.pistacking_info.append((stack.resnr, stack.restype, '%.2f' % stack.distance,
                                         '%.2f' % stack.angle, '%.2f' % stack.offset, stack.type, ",".join(ids)
                                         , tuple(stack.ligandring.center), tuple(stack.proteinring.center)))

        ### PI-CATION
        self.pication_features = ('RESNR', 'RESTYPE', 'DIST', 'OFFSET', 'PROTCHARGED', 'LIG_GROUP', 'LIG_IDX_LIST',
                                  'LIGCOO', 'PROTCOO')
        self.pication_info = []
        for picat in pli_class.pication_laro+pli_class.pication_paro:
            if picat.protcharged:
                ids = [str(lig_to_pdb[x.idx]) for x in picat.ring.atoms]
                group = 'Aromatic'
                self.pication_info.append((picat.resnr, picat.restype, '%.2f' % picat.distance,
                                           '%.2f' % picat.offset, picat.protcharged, group, ",".join(ids),
                                           tuple(picat.ring.center), tuple(picat.charge.center)))
            else:
                ids = [str(lig_to_pdb[x.idx]) for x in picat.charge.atoms]
                group = picat.charge.fgroup
                self.pication_info.append((picat.resnr, picat.restype, '%.2f' % picat.distance,
                                           '%.2f' % picat.offset, picat.protcharged, group, ",".join(ids),
                                           tuple(picat.charge.center), tuple(picat.ring.center)))

        ### HALOGEN BONDS
        self.halogen_features = ('RESNR', 'RESTYPE', 'DIST', 'DON_ANGLE', 'ACC_ANGLE', 'DONORTYPE', 'DON_IDX',
                                 'ACC_IDX', 'LIGCOO', 'PROTCOO')
        self.halogen_info = []
        for halogen in pli_class.halogen_bonds:
            self.halogen_info.append((halogen.resnr, halogen.restype, '%.2f' % halogen.distance,
                                      '%.2f' % halogen.don_angle, '%.2f' % halogen.acc_angle, halogen.donortype,
                                      lig_to_pdb[halogen.don.x.idx],
                                      mapping[halogen.acc.o.idx], halogen.acc.o.coords, halogen.don.x.coords))

    def write_section(self, name, features, info, f):
        """Provides formatting for one section (e.g. hydrogen bonds)"""
        if not len(info) == 0:
            f.write('\n\n### %s ###\n' % name)
            f.write('%s\n' % '\t'.join(features))
            for line in info:
                f.write('%s\n' % '\t'.join(map(str, line)))

    #@todo Write own code for tables
    def rst_table(self, array):
        #@todo Make cell_width dependent on maximum value in table column, not row
        cell_width = 2 + max(reduce(lambda x, y: x+y, [[len(item) for item in row] for row in array], []))
        num_cols = len(array[0])
        rst = num_cols*('+' + cell_width*'-') + '+\n'
        for i, row in enumerate(array):
            rst = rst + '| ' + '| '.join([(x + ((cell_width-1 - len(x)) * ' ')) for x in row]) + '|\n'
            if i == 0:
                rst += num_cols*('+' + cell_width*'=') + '+\n'
            else:
                rst += num_cols*('+' + cell_width*'-') + '+\n'
        return rst

    def generate_rst(self):
        """Generates an flat text report for a single binding site"""
        txt = []
        txt.append('%s' % self.name)
        txt.append("-"*len(self.name))
        for section in [['Hydrophobic Interactions', self.hydrophobic_features, self.hydrophobic_info],
                        ['Hydrogen Bonds', self.hbond_features, self.hbond_info],
                        ['Water-bridged Hydrogend Bonds', self.waterbridge_features, self.waterbridge_info],
                        ['Salt Bridges', self.saltbridge_features, self.saltbridge_info],
                        ['pi-Stacking', self.pistacking_features, self.pistacking_info],
                        ['pi-Cation Interaction', self.pication_features, self.pication_info],
                        ['Halogen Bonds', self.halogen_features, self.halogen_info]]:
            iname, features, info = section
            if not len(info) == 0:
                txt.append('\n**%s**' % iname)
                table = [features, ]
                for interaction in info:
                    values = []
                    for x in interaction:
                        if type(x) == str:
                            values.append(x)
                        elif type(x) == tuple and len(x) == 3:  # Coordinates
                            values.append("%.2f,%.2f,%.2f" % x)
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
        hetid = et.SubElement(identifiers, 'hetid')
        chain = et.SubElement(identifiers, 'chain')
        position = et.SubElement(identifiers, 'position')
        hetid.text, chain.text, position.text = self.name.split('-')
        interactions = et.SubElement(report, 'interactions')

        def format_interactions(element_name, features, interaction_information):
            """Returns a formatted element with interaction information."""
            interaction = et.Element(element_name)
            for j, single_contact in enumerate(interaction_information):
                new_contact = et.SubElement(interaction, element_name[:-1], id=str(j+1))
                for i, feature in enumerate(single_contact):
                    # Just assign the value unless it's an atom list, use subelements in this case
                    if features[i] == 'LIG_IDX_LIST':
                        feat = et.SubElement(new_contact, features[i].lower())
                        for k, atm_idx in enumerate(feature.split(',')):
                            idx = et.SubElement(feat, 'idx', id=str(k+1))
                            idx.text = str(atm_idx)
                    elif features[i] in ['LIGCOO', 'PROTCOO']:
                        feat = et.SubElement(new_contact, features[i].lower())
                        xc, yc, zc = feature
                        xcoo = et.SubElement(feat, 'x')
                        xcoo.text = str(xc)
                        ycoo = et.SubElement(feat, 'y')
                        ycoo.text = str(yc)
                        zcoo = et.SubElement(feat, 'z')
                        zcoo.text = str(zc)
                    else:
                        feat = et.SubElement(new_contact, features[i].lower())
                        feat.text = str(feature)
            return interaction

        interactions.append(format_interactions('hydrophobic_interactions', self.hydrophobic_features,
                                                self.hydrophobic_info))
        interactions.append(format_interactions('hydrogen_bonds', self.hbond_features, self.hbond_info))
        interactions.append(format_interactions('water_bridges', self.waterbridge_features, self.waterbridge_info))
        interactions.append(format_interactions('salt_bridges', self.saltbridge_features, self.saltbridge_info))
        interactions.append(format_interactions('pi_stacks', self.pication_features, self.pication_info))
        interactions.append(format_interactions('pi_cation_interactions', self.pication_features, self.pication_info))
        interactions.append(format_interactions('halogen_bonds', self.halogen_features, self.halogen_info))
        return report