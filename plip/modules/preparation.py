"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
preparation.py - Prepare PDB input files for processing.
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
from operator import itemgetter

# Own modules
from detection import *
from supplemental import *
import config

################
# MAIN CLASSES #
################


class Mol:

    def __init__(self, altconf):
        self.rings = None
        self.hydroph_atoms = None
        self.charged = None
        self.hbond_don_atom_pairs = None
        self.hbond_acc_atoms = None
        self.altconf = altconf

    def hydrophobic_atoms(self, all_atoms):
        """Select all carbon atoms which have only carbons and/or hydrogens as direct neighbors."""
        data = namedtuple('hydrophobic', 'atoms')
        atm = [a for a in all_atoms if a.atomicnum == 6 and set([natom.GetAtomicNum() for natom
                                                                in pybel.ob.OBAtomAtomIter(a.OBAtom)]).issubset({1, 6})]
        atm = [a for a in atm if a.idx not in self.altconf]
        return data(atoms=atm)

    def find_hba(self, all_atoms):
        """Find all possible hydrogen bond acceptors"""
        data = namedtuple('hbondacceptor', 'a type')
        a_set = []
        for atom in itertools.ifilter(lambda at: at.OBAtom.IsHbondAcceptor(), all_atoms):
            if atom.atomicnum not in [9, 17, 35, 53] and atom.idx not in self.altconf:  # Exclude halogen atoms
                a_set.append(data(a=atom, type='regular'))
        return a_set

    def find_hbd(self, all_atoms, hydroph_atoms):
        """Find all possible strong and weak hydrogen bonds donors (all hydrophobic C-H pairings)"""
        donor_pairs = []
        data = namedtuple('hbonddonor', 'd h type')
        for donor in [a for a in all_atoms if a.OBAtom.IsHbondDonor() and a.idx not in self.altconf]:
            in_ring = False
            if not in_ring:
                for adj_atom in [a for a in pybel.ob.OBAtomAtomIter(donor.OBAtom) if a.IsHbondDonorH()]:
                    donor_pairs.append(data(d=donor, h=pybel.Atom(adj_atom), type='regular'))
        for carbon in hydroph_atoms.atoms:
            for adj_atom in [a for a in pybel.ob.OBAtomAtomIter(carbon.OBAtom) if a.GetAtomicNum() == 1]:
                donor_pairs.append(data(d=carbon, h=pybel.Atom(adj_atom), type='weak'))
        return donor_pairs

    def find_rings(self, mol, all_atoms):
        """Find rings and return only aromatic."""
        data = namedtuple('aromatic_ring', 'atoms normal obj center type')
        rings, arings = [], []
        # Check here first for ligand rings not being detected as aromatic by Babel and check for planarity
        if len(mol.title) > 0:  # it's the ligand
            for ring in [r for r in mol.OBMol.GetSSSR()]:
                r_atoms = [a for a in all_atoms if ring.IsMember(a.OBAtom)]
                normals = []
                aromatic = True
                for a in r_atoms:
                    adj = pybel.ob.OBAtomAtomIter(a.OBAtom)
                    # Check for neighboring atoms in the ring
                    n_coords = [pybel.Atom(neigh).coords for neigh in adj if ring.IsMember(neigh)]
                    vec1, vec2 = vector(a.coords, n_coords[0]), vector(a.coords, n_coords[1])
                    normals.append(np.cross(vec1, vec2))
                # Given all normals of ring atoms and their neighbors, the angle between any has to be 5.0 deg or less
                for n1, n2 in itertools.product(normals, repeat=2):
                    arom_angle = vecangle(n1, n2)
                    if all([arom_angle > config.AROMATIC_PLANARITY, arom_angle < 180.0 - config.AROMATIC_PLANARITY]):
                        aromatic = False
                        break
                # Ring is aromatic either by OpenBabel's criteria or if sufficiently planar
                if aromatic or ring.IsAromatic():
                    arings.append(ring)
        else:
            # Detection for proteins should work more reliably just with OpenBabel
            arings = [r for r in mol.OBMol.GetSSSR() if r.IsAromatic()]

        # Store all rings which are detected as aromatic by Babel or are sufficiently planar
        for r in arings:
            r_atoms = [a for a in all_atoms if r.IsMember(a.OBAtom)]
            # Only consider rings with a minimum size of 5 atoms and restrict selection to avoid problems in
            # covalently bound ligands
            if 4 < len(r_atoms) <= 6:
                typ = r.GetType() if not r.GetType() == '' else 'unknown'
                ring_atms = [r_atoms[a].coords for a in [0, 2, 4]]  # Probe atoms for normals, assuming planarity
                ringv1 = vector(ring_atms[0], ring_atms[1])
                ringv2 = vector(ring_atms[2], ring_atms[0])
                rings.append(data(atoms=r_atoms,
                                  normal=normalize_vector(np.cross(ringv1, ringv2)),
                                  obj=r,
                                  center=centroid([ra.coords for ra in r_atoms]),
                                  type=typ))
        return rings

    def get_hydrophobic_atoms(self):
        return self.hydroph_atoms

    def get_hba(self):
        return self.hbond_acc_atoms

    def get_hbd(self):
        return [don_pair for don_pair in self.hbond_don_atom_pairs if don_pair.type == 'regular']

    def get_weak_hbd(self):
        return [don_pair for don_pair in self.hbond_don_atom_pairs if don_pair.type == 'weak']

    def get_pos_charged(self):
        return [charge for charge in self.charged if charge.type == 'positive']

    def get_neg_charged(self):
        return [charge for charge in self.charged if charge.type == 'negative']


class PLInteraction:
    """Class to store a ligand, a protein and their interactions."""
    def __init__(self, lig_obj, bs_obj, protcomplex):
        """Detect all interactions when initializing"""
        self.ligand = lig_obj
        self.name = lig_obj.name
        self.lig_members = lig_obj.members
        self.pdbid = protcomplex.pymol_name
        self.bindingsite = bs_obj
        self.idx_to_pdb = protcomplex.idx_to_pdb_mapping
        self.lig_to_pdb = lig_obj.pymol_data.maptopdb
        self.output_path = protcomplex.output_path
        self.altconf = protcomplex.altconf
        # #@todo Refactor code to combine different directionality

        self.saltbridge_lneg = saltbridge(self.bindingsite.get_pos_charged(), self.ligand.get_neg_charged(), True)
        self.saltbridge_pneg = saltbridge(self.ligand.get_pos_charged(), self.bindingsite.get_neg_charged(), False)

        self.all_hbonds_ldon = hbonds(self.bindingsite.get_hba(),
                                      self.ligand.get_hbd(), False, 'strong')
        self.all_hbonds_pdon = hbonds(self.ligand.get_hba(),
                                      self.bindingsite.get_hbd(), True, 'strong')

        self.hbonds_ldon = self.refine_hbonds_ldon(self.all_hbonds_ldon, self.saltbridge_lneg,
                                                   self.saltbridge_pneg)
        self.hbonds_pdon = self.refine_hbonds_pdon(self.all_hbonds_pdon, self.saltbridge_lneg,
                                                   self.saltbridge_pneg)

        self.pistacking = pistacking(self.bindingsite.rings, self.ligand.rings)

        self.all_pi_cation_laro = pication(self.ligand.rings, self.bindingsite.get_pos_charged(), True)
        self.pication_paro = pication(self.bindingsite.rings, self.ligand.get_pos_charged(), False)

        self.pication_laro = self.refine_pi_cation_laro(self.all_pi_cation_laro, self.pistacking)

        self.all_hydrophobic_contacts = hydrophobic_interactions(self.bindingsite.get_hydrophobic_atoms(),
                                                                 self.ligand.get_hydrophobic_atoms())
        self.hydrophobic_contacts = self.refine_hydrophobic(self.all_hydrophobic_contacts, self.pistacking)
        self.halogen_bonds = halogen(self.bindingsite.halogenbond_acc, self.ligand.halogenbond_don)
        self.water_bridges = water_bridges(self.bindingsite.get_hba(), self.ligand.get_hba(),
                                           self.bindingsite.get_hbd(), self.ligand.get_hbd(),
                                           self.ligand.water)

        self.water_bridges = self.refine_water_bridges(self.water_bridges, self.hbonds_ldon, self.hbonds_pdon)

        self.metal_complexes = metal_complexation(self.ligand.metals, self.ligand.metal_binding,
                                                  self.bindingsite.metal_binding)

        self.all_itypes = self.saltbridge_lneg + self.saltbridge_pneg + self.hbonds_pdon + self.hbonds_ldon \
                          + self.pistacking + self.pication_laro + self.pication_paro + self.hydrophobic_contacts \
                          + self.halogen_bonds + self.water_bridges + self.metal_complexes

        self.no_interactions = all(len(i) == 0 for i in self.all_itypes)

        # Exclude empty chains (coming from ligand as a target, from metal complexes)
        self.interacting_chains = sorted(list(set([i.reschain for i in self.all_itypes if not i.reschain == ' '])))

        self.interacting_res = list(set([''.join([str(i.resnr), str(i.reschain)]) for i in self.all_itypes]))
        if len(self.interacting_res) != 0:
            message('Ligand interacts with %i binding site residue(s) in chain(s) %s.\n'
                    % (len(self.interacting_res), '/'.join(self.interacting_chains)), indent=True)
            interactions_list = []
            num_saltbridges = len(self.saltbridge_lneg + self.saltbridge_pneg)
            num_hbonds = len(self.hbonds_ldon + self.hbonds_pdon)
            num_pication = len(self.pication_laro + self.pication_paro)
            num_pistack = len(self.pistacking)
            num_halogen = len(self.halogen_bonds)
            num_waterbridges = len(self.water_bridges)
            if num_saltbridges != 0:
                interactions_list.append('%i salt bridge(s)' % num_saltbridges)
            if num_hbonds != 0:
                interactions_list.append('%i hydrogen bond(s)' % num_hbonds)
            if num_pication != 0:
                interactions_list.append('%i pi-cation interaction(s)' % num_pication)
            if num_pistack != 0:
                interactions_list.append('%i pi-stacking(s)' % num_pistack)
            if num_halogen != 0:
                interactions_list.append('%i halogen bond(s)' % num_halogen)
            if num_waterbridges != 0:
                interactions_list.append('%i water bridge(s)' % num_waterbridges)
            if not len(interactions_list) == 0:
                message('Complex uses %s.\n' % ', '.join(interactions_list), indent=True)
        else:
            message('No interactions for this ligand.\n', indent=True)

    def refine_hydrophobic(self, all_h, pistacks):
        """Apply several rules to reduce the number of hydrophobic interactions."""
        sel = {}
        #  1. Rings interacting via stacking can't have additional hydrophobic pliprofiler between each other.
        for pistack, h in itertools.product(pistacks, all_h):
            h1, h2 = h.bsatom.idx, h.ligatom.idx
            brs, lrs = [p1.idx for p1 in pistack.proteinring.atoms], [p2.idx for p2 in pistack.ligandring.atoms]
            if h1 in brs and h2 in lrs:
                sel[(h1, h2)] = "EXCLUDE"
        hydroph = [h for h in all_h if not (h.bsatom.idx, h.ligatom.idx) in sel]
        sel2 = {}
        #  2. If a ligand atom interacts with several binding site atoms in the same residue,
        #  keep only the one with the closest distance
        for h in hydroph:
            if not (h.ligatom.idx, h.resnr) in sel2:
                sel2[(h.ligatom.idx, h.resnr)] = h
            else:
                if sel2[(h.ligatom.idx, h.resnr)].distance > h.distance:
                    sel2[(h.ligatom.idx, h.resnr)] = h
        hydroph = [h for h in sel2.values()]
        hydroph_final = []
        bsclust = {}
        #  3. If a protein atom interacts with several neighboring ligand atoms, just keep the one with the closest dist
        for h in hydroph:
            if h.bsatom.idx not in bsclust:
                bsclust[h.bsatom.idx] = [h, ]
            else:
                bsclust[h.bsatom.idx].append(h)

        idx_to_h = {}
        for bs in [a for a in bsclust if len(bsclust[a]) == 1]:
            hydroph_final.append(bsclust[bs][0])

        # A list of tuples with the idx of an atom and one of its neighbours is created
        for bs in [a for a in bsclust if not len(bsclust[a]) == 1]:
            tuples = []
            all_idx = [i.ligatom.idx for i in bsclust[bs]]
            for b in bsclust[bs]:
                idx = b.ligatom.idx
                neigh = [na for na in pybel.ob.OBAtomAtomIter(b.ligatom.OBAtom)]
                for n in neigh:
                    n_idx = n.GetIdx()
                    if n_idx in all_idx:
                        if n_idx < idx:
                            tuples.append((n_idx, idx))
                        else:
                            tuples.append((idx, n_idx))
                        idx_to_h[idx] = b

            tuples = list(set(tuples))
            tuples = sorted(tuples, key=itemgetter(1))
            clusters = cluster_doubles(tuples)  # Cluster connected atoms (i.e. find hydrophobic patches)

            for cluster in clusters:
                min_dist = float('inf')
                min_h = None
                for atm_idx in cluster:
                    h = idx_to_h[atm_idx]
                    if h.distance < min_dist:
                        min_dist = h.distance
                        min_h = h
                hydroph_final.append(min_h)
        before, reduced = len(all_h), len(hydroph_final)
        if not before == 0 and not before == reduced:
            message('Reduced number of hydrophobic contacts from %i to %i.\n' % (before, reduced), indent=True)
        return hydroph_final

    def refine_hbonds_ldon(self, all_hbonds, salt_lneg, salt_pneg):
        """Refine selection of hydrogen bonds. Do not allow groups which already form salt bridges to form H-Bonds."""
        i_set = {}
        for hbond in all_hbonds:
            i_set[hbond] = False
            for salt in salt_pneg:
                protidx, ligidx = [at.idx for at in salt.negative.atoms], [at.idx for at in salt.positive.atoms]
                if hbond.d.idx in ligidx and hbond.a.idx in protidx:
                    i_set[hbond] = True
            for salt in salt_lneg:
                protidx, ligidx = [at.idx for at in salt.positive.atoms], [at.idx for at in salt.negative.atoms]
                if hbond.d.idx in ligidx and hbond.a.idx in protidx:
                    i_set[hbond] = True

        # Allow only one hydrogen bond per donor, select interaction with larger donor angle
        second_set = {}
        hbls = [k for k in i_set.keys() if not i_set[k]]
        for hbl in hbls:
            if hbl.d.idx not in second_set:
                second_set[hbl.d.idx] = (hbl.angle, hbl)
            else:
                if second_set[hbl.d.idx][0] < hbl.angle:
                    second_set[hbl.d.idx] = (hbl.angle, hbl)
        return [hb[1] for hb in second_set.values()]

    def refine_hbonds_pdon(self, all_hbonds, salt_lneg, salt_pneg):
        """Refine selection of hydrogen bonds. Do not allow groups which already form salt bridges to form H-Bonds with
        atoms of the same group.
        """
        i_set = {}
        for hbond in all_hbonds:
            i_set[hbond] = False
            for salt in salt_lneg:
                protidx, ligidx = [at.idx for at in salt.positive.atoms], [at.idx for at in salt.negative.atoms]
                if hbond.a.idx in ligidx and hbond.d.idx in protidx:
                    i_set[hbond] = True
            for salt in salt_pneg:
                protidx, ligidx = [at.idx for at in salt.negative.atoms], [at.idx for at in salt.positive.atoms]
                if hbond.a.idx in ligidx and hbond.d.idx in protidx:
                    i_set[hbond] = True

        # Allow only one hydrogen bond per donor, select interaction with larger donor angle
        second_set = {}
        hbps = [k for k in i_set.keys() if not i_set[k]]
        for hbp in hbps:
            if hbp.d.idx not in second_set:
                second_set[hbp.d.idx] = (hbp.angle, hbp)
            else:
                if second_set[hbp.d.idx][0] < hbp.angle:
                    second_set[hbp.d.idx] = (hbp.angle, hbp)
        return [hb[1] for hb in second_set.values()]

    def refine_pi_cation_laro(self, all_picat, stacks):
        """Just important for constellations with histidine involved. If the histidine ring is positioned in stacking
        position to an aromatic ring in the ligand, there is in most cases stacking and pi-cation interaction reported
        as histidine also carries a positive charge in the ring. For such cases, only report stacking.
        """
        i_set = []
        for picat in all_picat:
            exclude = False
            for stack in stacks:
                if whichrestype(stack.proteinring.atoms[0]) == 'HIS' and picat.ring.obj == stack.ligandring.obj:
                    exclude = True
            if not exclude:
                i_set.append(picat)
        return i_set

    def refine_water_bridges(self, wbridges, hbonds_ldon, hbonds_pdon):
        """A donor atom already forming a hydrogen bond is not allowed to form a water bridge. Each water molecule
        can only be donor for two water bridges, selecting the constellation with the omega angle closest to 110 deg."""
        donor_atoms_hbonds = [hb.d.idx for hb in hbonds_ldon + hbonds_pdon]
        wb_dict = {}
        wb_dict2 = {}
        omega = 110.0

        # Just one hydrogen bond per donor atom
        for wbridge in [wb for wb in wbridges if wb.d.idx not in donor_atoms_hbonds]:
            if (wbridge.water.idx, wbridge.a.idx) not in wb_dict:
                wb_dict[(wbridge.water.idx, wbridge.a.idx)] = wbridge
            else:
                if abs(omega - wb_dict[(wbridge.water.idx, wbridge.a.idx)].w_angle) < abs(omega - wbridge.w_angle):
                    wb_dict[(wbridge.water.idx, wbridge.a.idx)] = wbridge
        for wb_tuple in wb_dict:
            water, acceptor = wb_tuple
            if water not in wb_dict2:
                wb_dict2[water] = [(abs(omega - wb_dict[wb_tuple].w_angle), wb_dict[wb_tuple]), ]
            elif len(wb_dict2[water]) == 1:
                wb_dict2[water].append((abs(omega - wb_dict[wb_tuple].w_angle), wb_dict[wb_tuple]))
                wb_dict2[water] = sorted(wb_dict2[water])
            else:
                if wb_dict2[water][1][0] < abs(omega - wb_dict[wb_tuple].w_angle):
                    wb_dict2[water] = [wb_dict2[water][0], (wb_dict[wb_tuple].w_angle, wb_dict[wb_tuple])]

        filtered_wb = []
        for fwbridges in wb_dict2.values():
            [filtered_wb.append(fwb[1]) for fwb in fwbridges]
        return filtered_wb


class BindingSite(Mol):
    def __init__(self, atoms, protcomplex, cclass, altconf, min_dist):
        """Find all relevant parts which could take part in interactions"""
        Mol.__init__(self, altconf)
        self.complex = cclass
        self.full_mol = protcomplex
        self.all_atoms = atoms
        self.min_dist = min_dist # Minimum distance of bs res to ligand
        self.bs_res = list(set([''.join([str(whichresnumber(a)), whichchain(a)]) for a in self.all_atoms]))  # e.g. 47A
        self.rings = self.find_rings(self.full_mol, self.all_atoms)
        self.hydroph_atoms = self.hydrophobic_atoms(self.all_atoms)
        self.hbond_acc_atoms = self.find_hba(self.all_atoms)
        self.hbond_don_atom_pairs = self.find_hbd(self.all_atoms, self.hydroph_atoms)
        self.charged = self.find_charged(self.full_mol)
        self.halogenbond_acc = self.find_hal(self.all_atoms)
        self.metal_binding = self.find_metal_binding(self.full_mol)

    def find_hal(self, atoms):
        """Look for halogen bond acceptors (Y-{O|P|N|S}, with Y=C,P,S)"""
        data = namedtuple('hal_acceptor', 'o y')
        a_set = []
        # All oxygens, nitrogen, sulfurs with neighboring carbon, phosphor, nitrogen or sulfur
        for a in [at for at in atoms if at.atomicnum in [8, 7, 16]]:
            n_atoms = [na for na in pybel.ob.OBAtomAtomIter(a.OBAtom) if na.GetAtomicNum() in [6, 7, 15, 16]]
            if len(n_atoms) == 1:  # Proximal atom
                a_set.append(data(o=a, y=pybel.Atom(n_atoms[0])))
        return a_set

    def find_charged(self, mol):
        """Looks for positive charges in arginine, histidine or lysine, for negative in aspartic and glutamic acid."""
        data = namedtuple('pcharge', 'atoms type center restype resnr reschain')
        a_set = []
        for res in pybel.ob.OBResidueIter(mol.OBMol):
            a_contributing = []
            if res.GetName() in ('ARG', 'HIS', 'LYS'):  # Arginine, Histidine or Lysine have charged sidechains
                for a in pybel.ob.OBResidueAtomIter(res):
                    if a.GetType().startswith('N') and res.GetAtomProperty(a, 8) \
                            and not self.complex.idx_to_pdb_mapping[a.GetIdx()] in self.altconf:
                        a_contributing.append(pybel.Atom(a))
                if not len(a_contributing) == 0:
                    a_set.append(data(atoms=a_contributing,
                                      type='positive',
                                      center=centroid([ac.coords for ac in a_contributing]),
                                      restype=res.GetName(),
                                      resnr=res.GetNum(),
                                      reschain=res.GetChain()))
            if res.GetName() in ('GLU', 'ASP'):  # Aspartic or Glutamic Acid
                for a in pybel.ob.OBResidueAtomIter(res):
                    if a.GetType().startswith('O') and res.GetAtomProperty(a, 8) \
                            and not self.complex.idx_to_pdb_mapping[a.GetIdx()] in self.altconf:
                        a_contributing.append(pybel.Atom(a))
                if not len(a_contributing) == 0:
                    a_set.append(data(atoms=a_contributing,
                                      type='negative',
                                      center=centroid([ac.coords for ac in a_contributing]),
                                      restype=res.GetName(),
                                      resnr=res.GetNum(),
                                      reschain=res.GetChain()))
        return a_set

    def find_metal_binding(self, mol):
        """Looks for atoms that could possibly be involved in chelating a metal ion.
        This can be any main chain oxygen atom or oxygen, nitrogen and sulfur from specific amino acids"""
        data = namedtuple('metal_binding', 'atom type restype resnr reschain location')
        a_set = []
        for res in pybel.ob.OBResidueIter(mol.OBMol):
            restype, reschain, resnr = res.GetName().upper(), res.GetChain(), res.GetNum()
            if restype in ['ASP', 'GLU', 'SER', 'THR', 'TYR']:  # Look for oxygens here
                for a in pybel.ob.OBResidueAtomIter(res):
                    if a.GetType().startswith('O') and res.GetAtomProperty(a, 8) \
                            and not self.complex.idx_to_pdb_mapping[a.GetIdx()] in self.altconf:
                                a_set.append(data(atom=pybel.Atom(a), type='O', restype=restype,
                                                  resnr=resnr, reschain=reschain,
                                                  location='protein.sidechain'))
            if restype == 'HIS':  # Look for nitrogen here
                for a in pybel.ob.OBResidueAtomIter(res):
                    if a.GetType().startswith('N') and res.GetAtomProperty(a, 8) \
                            and not self.complex.idx_to_pdb_mapping[a.GetIdx()] in self.altconf:
                                a_set.append(data(atom=pybel.Atom(a), type='N', restype=restype,
                                                  resnr=resnr, reschain=reschain,
                                                  location='protein.sidechain'))
            if restype == 'CYS':  # Look for sulfur here
                for a in pybel.ob.OBResidueAtomIter(res):
                    if a.GetType().startswith('S') and res.GetAtomProperty(a, 8) \
                            and not self.complex.idx_to_pdb_mapping[a.GetIdx()] in self.altconf:
                                a_set.append(data(atom=pybel.Atom(a), type='S', restype=restype,
                                                  resnr=resnr, reschain=reschain,
                                                  location='protein.sidechain'))
            for a in pybel.ob.OBResidueAtomIter(res):  # All main chain oxygens
                if a.GetType().startswith('O') and res.GetAtomProperty(a, 2) \
                and not self.complex.idx_to_pdb_mapping[a.GetIdx()] in self.altconf and restype != 'HOH':
                    a_set.append(data(atom=pybel.Atom(a), type='O', restype=res.GetName(),
                                      resnr=res.GetNum(), reschain=res.GetChain(),
                                      location='protein.mainchain'))
        return a_set


class Ligand(Mol):
    def __init__(self, cclass, ligand):
        altconf = cclass.altconf
        Mol.__init__(self, altconf)
        mapping = ligand.mapping
        members = ligand.members
        water = ligand.water
        self.longname = ligand.longname
        self.type = ligand.type
        self.complex = cclass
        self.molecule = ligand.mol  # Pybel Molecule
        self.smile = self.molecule.write(format='smi')
        if not len(self.smile) == 0:
            self.smile = self.smile.split()[0]
        else:
            message('[Warning] Could not write smile for this ligand.\n', indent=True)
            self.smile = ''
        self.heavy_atoms = self.molecule.OBMol.NumHvyAtoms()  # Heavy atoms count
        self.name = self.molecule.title
        self.all_atoms = self.molecule.atoms
        self.atmdict = {l.idx: l for l in self.all_atoms}
        self.rings = self.find_rings(self.molecule, self.all_atoms)
        self.hydroph_atoms = self.hydrophobic_atoms(self.all_atoms)
        self.hbond_acc_atoms = self.find_hba(self.all_atoms)
        num_rings = len(self.rings)
        if num_rings != 0:
            message('Contains %i aromatic ring(s).\n' % num_rings, indent=True)

        ##########################################################
        # Special Case for hydrogen bond acceptor identification #
        ##########################################################

        self.inverse_mapping = {v: k for k, v in mapping.items()}
        self.pdb_to_idx_mapping = {v: k for k, v in cclass.idx_to_pdb_mapping.items()}
        self.hbond_don_atom_pairs = self.find_hbd(self.all_atoms, self.hydroph_atoms)

        ######
        donor_pairs = []
        data = namedtuple('hbonddonor', 'd h type')
        for donor in self.all_atoms:
            pdbidx = cclass.idx_to_pdb_mapping[mapping[donor.idx]]  # Work with protonated atoms for HBD search
            d = cclass.atoms[self.pdb_to_idx_mapping[pdbidx]]
            if d.OBAtom.IsHbondDonor():
                for adj_atom in [a for a in pybel.ob.OBAtomAtomIter(d.OBAtom) if a.IsHbondDonorH()]:
                    donor_pairs.append(data(d=donor, h=pybel.Atom(adj_atom), type='regular'))
        self.hbond_don_atom_pairs = donor_pairs
        #######

        self.charged = self.find_charged(self.all_atoms)
        self.centroid = centroid([a.coords for a in self.all_atoms])
        self.max_dist_to_center = max((euclidean3d(self.centroid, a.coords) for a in self.all_atoms))
        self.water = []
        self.members = members
        for hoh in water:
            oxy = None
            for at in pybel.ob.OBResidueAtomIter(hoh):
                if at.GetAtomicNum() == 8 and at.GetIdx() not in self.altconf:
                    oxy = pybel.Atom(at)
            # There are some cases where there is no oxygen in a water residue, ignore those
            if not set([at.GetAtomicNum() for at in pybel.ob.OBResidueAtomIter(hoh)]) == {1} and oxy is not None:
                if euclidean3d(self.centroid, oxy.coords) < self.max_dist_to_center + config.BS_DIST:
                    self.water.append(oxy)
        s = self.name.split('-')
        data = namedtuple('pymol_data', 'hetid chain resid maptopdb bs_id')
        self.pymol_data = data(hetid=s[0], chain=s[1], resid=s[2], maptopdb=mapping, bs_id=s)
        self.halogenbond_don = self.find_hal(self.all_atoms)
        self.metal_binding = self.find_metal_binding(self.all_atoms, self.water)
        # #@todo Update documentation
        self.metals = [a for a in self.all_atoms if a.type.upper() in config.METAL_IONS]

    def is_functional_group(self, atom, group):
        """Given a pybel atom, look up if it belongs to a function group"""
        n_atoms = [a_neighbor.GetAtomicNum() for a_neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom)]

        if group in ['quartamine', 'tertamine'] and atom.atomicnum == 7:  # Nitrogen
            # It's a nitrogen, so could be a protonated amine or quaternary ammonium
            if '1' not in n_atoms and len(n_atoms) == 4:
                return True if group == 'quartamine' else False  # It's a quat. ammonium (N with 4 residues != H)
            elif atom.OBAtom.GetHyb() == 3 and len(n_atoms) >= 3:
                return True if group == 'tertamine' else False  # It's sp3-hybridized, so could pick up an hydrogen
            else:
                return False

        if group in ['sulfonium', 'sulfonicacid', 'sulfate'] and atom.atomicnum == 16:  # Sulfur
            if '1' not in n_atoms and len(n_atoms) == 3:  # It's a sulfonium (S with 3 residues != H)
                return True if group == 'sulfonium' else False
            elif n_atoms.count(8) == 3:  # It's a sulfonate or sulfonic acid
                return True if group == 'sulfonicacid' else False
            elif n_atoms.count(8) == 4:  # It's a sulfate
                return True if group == 'sulfate' else False

        if group == 'phosphate' and atom.atomicnum == 15:  # Phosphor
            if set(n_atoms) == {8}:  # It's a phosphate
                return True

        if group in ['carboxylate', 'guanidine'] and atom.atomicnum == 6:  # It's a carbon atom
            if n_atoms.count(8) == 2 and n_atoms.count(6) == 1:  # It's a carboxylate group
                return True if group == 'carboxylate' else False
            elif n_atoms.count(7) == 3 and len(n_atoms) == 3:  # It's a guanidine group
                nitro_partners = []
                for nitro in pybel.ob.OBAtomAtomIter(atom.OBAtom):
                    nitro_partners.append(len([b_neighbor for b_neighbor in pybel.ob.OBAtomAtomIter(nitro)]))
                if min(nitro_partners) == 1:  # One nitrogen is only connected to the carbon, can pick up a H
                    return True if group == 'guanidine' else False

        if group == 'halocarbon' and atom.atomicnum in [9, 17, 35, 53]:  # Halogen atoms
            n_atoms = [na for na in pybel.ob.OBAtomAtomIter(atom.OBAtom) if na.GetAtomicNum() == 6]
            if len(n_atoms) == 1:  # Halocarbon
                return True
        else:
            return False

    def find_hal(self, atoms):
        """Look for halogen bond donors (X-C, with X=F, Cl, Br, I)"""
        data = namedtuple('hal_donor', 'x c')
        a_set = []
        for a in atoms:
            if self.is_functional_group(a, 'halocarbon'):
                n_atoms = [na for na in pybel.ob.OBAtomAtomIter(a.OBAtom) if na.GetAtomicNum() == 6]
                a_set.append(data(x=a, c=pybel.Atom(n_atoms[0])))
        if len(a_set) != 0:
            message('Ligand contains %i halogen atom(s).\n' % len(a_set), indent=True)
        return a_set

    def find_charged(self, all_atoms):
        """Identify all positively charged groups in a ligand. This search is not exhaustive, as the cases can be quite
        diverse. The typical cases seem to be protonated amines, quaternary ammoinium and sulfonium
        as mentioned in 'Cation-pi interactions in ligand recognition and catalysis' (Zacharias et al., 2002)).
        Identify negatively charged groups in the ligand.
        """
        data = namedtuple('lcharge', 'atoms type center fgroup')
        a_set = []
        for a in all_atoms:
            if self.is_functional_group(a, 'quartamine'):
                a_set.append(data(atoms=[a, ], type='positive', center=list(a.coords), fgroup='quartamine'))
            elif self.is_functional_group(a, 'tertamine'):
                a_set.append(data(atoms=[a, ], type='positive', center=list(a.coords), fgroup='tertamine'))
            if self.is_functional_group(a, 'sulfonium'):
                a_set.append(data(atoms=[a, ], type='positive', center=list(a.coords), fgroup='sulfonium'))
            if self.is_functional_group(a, 'phosphate'):
                a_contributing = [a, ]
                [a_contributing.append(pybel.Atom(neighbor)) for neighbor in pybel.ob.OBAtomAtomIter(a.OBAtom)]
                a_set.append(data(atoms=a_contributing, type='negative', center=a.coords, fgroup='phosphate'))
            if self.is_functional_group(a, 'sulfonicacid'):
                a_contributing = [a, ]
                [a_contributing.append(pybel.Atom(neighbor)) for neighbor in pybel.ob.OBAtomAtomIter(a.OBAtom) if
                 neighbor.GetAtomicNum() == 8]
                a_set.append(data(atoms=a_contributing, type='negative', center=a.coords, fgroup='sulfonicacid'))
            elif self.is_functional_group(a, 'sulfate'):
                a_contributing = [a, ]
                [a_contributing.append(pybel.Atom(neighbor)) for neighbor in pybel.ob.OBAtomAtomIter(a.OBAtom)]
                a_set.append(data(atoms=a_contributing, type='negative', center=a.coords, fgroup='sulfate'))
            if self.is_functional_group(a, 'carboxylate'):
                a_contributing = [pybel.Atom(neighbor) for neighbor in pybel.ob.OBAtomAtomIter(a.OBAtom)
                                  if neighbor.GetAtomicNum() == 8]
                a_set.append(data(atoms=a_contributing, type='negative',
                                  center=centroid([a.coords for a in a_contributing]), fgroup='carboxylate'))
            elif self.is_functional_group(a, 'guanidine'):
                a_contributing = [pybel.Atom(neighbor) for neighbor in pybel.ob.OBAtomAtomIter(a.OBAtom)
                                  if neighbor.GetAtomicNum() == 7]
                a_set.append(data(atoms=a_contributing, type='positive', center=a.coords, fgroup='guanidine'))
        return a_set

    def find_metal_binding(self, lig_atoms, water_oxygens):
        """Looks for atoms that could possibly be involved in binding a metal ion.
        This can be any water oxygen, as well as oxygen from carboxylate, phophoryl, phenolate, alcohol;
        nitrogen from imidazole; sulfur from thiolate.
        """
        a_set = []
        data = namedtuple('metal_binding', 'atom type fgroup restype resnr reschain location')
        for oxygen in water_oxygens:
            a_set.append(data(atom=oxygen, type='O', fgroup='water', restype=whichrestype(oxygen),
                              resnr=whichresnumber(oxygen), reschain=whichchain(oxygen), location='water'))
        # #@todo Check detection
        # #@todo Refactor code
        for a in lig_atoms:
            n_atoms = pybel.ob.OBAtomAtomIter(a.OBAtom)  # Neighboring atoms
            # All atomic numbers of neighboring atoms
            n_atoms_atomicnum = [n.GetAtomicNum() for n in pybel.ob.OBAtomAtomIter(a.OBAtom)]
            if a.atomicnum == 8:  # Oxygen
                if n_atoms_atomicnum.count('1') == 1 and len(n_atoms_atomicnum) == 2:  # Oxygen in alcohol (R-[O]-H)
                    a_set.append(data(atom=a, type='O', fgroup='alcohol', restype=whichrestype(a),
                                      resnr=whichresnumber(a), reschain=whichchain(a), location='ligand'))
                if True in [n.IsAromatic() for n in n_atoms] and not a.OBAtom.IsAromatic():  # Phenolate oxygen
                    a_set.append(data(atom=a, type='O', fgroup='phenolate', restype=whichrestype(a),
                                      resnr=whichresnumber(a), reschain=whichchain(a), location='ligand'))
            if a.atomicnum == 6:  # It's a carbon atom
                if n_atoms_atomicnum.count(8) == 2 and n_atoms_atomicnum.count(6) == 1:  # It's a carboxylate group
                    for neighbor in [n for n in n_atoms if n.GetAtomicNum() == 8]:
                        a_set.append(data(atom=pybel.Atom(neighbor), type='O', fgroup='carboxylate', restype=whichrestype(neighbor),
                                          resnr=whichresnumber(neighbor), reschain=whichchain(neighbor),
                                          location='ligand'))
            if a.atomicnum == 15:  # It's a phosphor atom
                if n_atoms_atomicnum.count(8) >= 3:  # It's a phosphoryl
                    for neighbor in [n for n in n_atoms if n.GetAtomicNum() == 8]:
                        a_set.append(data(atom=pybel.Atom(neighbor), type='O', fgroup='phosphoryl', restype=whichrestype(neighbor),
                                          resnr=whichresnumber(neighbor), reschain=whichchain(neighbor),
                                          location='ligand'))
                if n_atoms_atomicnum.count(8) == 2:  # It's another phosphor-containing group #@todo (correct name?)
                    for neighbor in [n for n in n_atoms if n.GetAtomicNum() == 8]:
                        a_set.append(data(atom=pybel.Atom(neighbor), type='O', fgroup='phosphor.other',
                                          restype=whichrestype(neighbor), resnr=whichresnumber(neighbor),
                                          reschain=whichchain(neighbor), location='ligand'))
            if a.atomicnum == 7:  # It's a nitrogen atom
                if n_atoms_atomicnum.count(6) == 2:  # It's imidazole/pyrrole or similar
                    a_set.append(data(atom=a, type='N', fgroup='imidazole/pyrrole', restype=whichrestype(a),
                                      resnr=whichresnumber(a), reschain=whichchain(a), location='ligand'))
            if a.atomicnum == 16:  # It's a sulfur atom
                if True in [n.IsAromatic() for n in n_atoms] and not a.OBAtom.IsAromatic():  # Thiolate
                    a_set.append(data(atom=a, type='S', fgroup='thiolate', restype=whichrestype(a),
                                      resnr=whichresnumber(a), reschain=whichchain(a), location='ligand'))
                if set(n_atoms_atomicnum) == {26}:  # Sulfur in Iron sulfur cluster
                    a_set.append(data(atom=a, type='S', fgroup='iron-sulfur.cluster', restype=whichrestype(a),
                                      resnr=whichresnumber(a), reschain=whichchain(a), location='ligand'))

        return a_set


class PDBComplex:
    """Contains a collection of objects associated with a PDB complex, i.e. one or several ligands and their binding
    sites as well as information about the pliprofiler between them. Provides functions to load and prepare input files
    such as PDB files.
    """

    def __init__(self):
        self.interaction_sets = {}  # Dictionary with site identifiers as keys and object as value
        self.protcomplex = None
        self.filetype = None
        self.atoms = {}  # Dictionary of Pybel atoms, accessible by their idx
        self.sourcefiles = {}
        self.output_path = '/tmp'
        self.pymol_name = None
        self.idx_to_pdb_mapping = {}
        self.modres = set()
        self.altconf = []  # Atom idx of atoms with alternate conformations
        self.covalent = []  # Covalent linkages between ligands and protein residues/other ligands
        self.excluded = []  # Excluded ligands

    def load_pdb(self, pdbpath):
        """Loads a pdb file with protein AND ligand(s), separates and prepares them."""
        self.sourcefiles['pdbcomplex'] = pdbpath
        self.protcomplex, self.filetype = read_pdb(pdbpath)
        message('PDB structure successfully read.\n')

        # Counting is different from PDB if TER records present
        self.idx_to_pdb_mapping, self.modres, self.covalent = parse_pdb(open(pdbpath).readlines())
        # #@todo Include this in the parse_pdb function, return named tuple?
        self.altconf = get_altconf_atoms(open(pdbpath).readlines())
        try:
            self.pymol_name = self.protcomplex.data['HEADER'][56:60].lower()  # Get name from HEADER data
        except KeyError:  # Extract the PDBID from the filename
            self.pymol_name = extract_pdbid(pdbpath.split('/')[-1])
        self.protcomplex.OBMol.AddPolarHydrogens()
        for atm in self.protcomplex:
            self.atoms[atm.idx] = atm

        # Extract and prepare ligands
        ligands, excluded = getligs(self.protcomplex, self.altconf, self.idx_to_pdb_mapping, self.modres, self.covalent)
        self.excluded = excluded
        if len(excluded) != 0:
            message("Excluded molecules as ligands: %s\n" % ','.join([lig for lig in excluded]))

        resis = [obres for obres in pybel.ob.OBResidueIter(self.protcomplex.OBMol) if obres.GetResidueProperty(0)]

        num_ligs = len(ligands)
        if num_ligs == 1:
            message("Analyzing one ligand...\n")
        elif num_ligs > 1:
            message("Analyzing %i ligands...\n" % num_ligs)
        else:
            message("Structure contains no ligands.\n")

        for ligand in ligands:
            single_sites = []
            for member in ligand.members:
                single_sites.append('-'.join([str(x) for x in member]))
            site = ' : '.join(single_sites)
            site = site if not len(site) > 20 else site[:20] + '...'
            longname = ligand.longname if not len(ligand.longname) > 20 else ligand.longname[:20] + '...'
            ligtype = 'Unspecified type' if ligand.type == 'UNSPECIFIED' else ligand.type
            ligtext = "\n%s [%s] -- %s" % (longname, ligtype, site)
            any_in_biolip = len(set([x[0] for x in ligand.members]).intersection(config.biolip_list)) != 0
            if ligtype not in ['POLYMER', 'DNA', 'ION', 'DNA+ION', 'RNA+ION', 'SMALLMOLECULE+ION'] and any_in_biolip:
                ligtext += ' (possible artifact/unspecific binder)'
            message(ligtext)
            message('\n' + '-' * len(ligtext) + '\n')

            lig_obj = Ligand(self, ligand)
            cutoff = lig_obj.max_dist_to_center + config.BS_DIST
            bs_res = self.extract_bs(cutoff, lig_obj.centroid, resis)
            # Get a list of all atoms belonging to the binding site, search by idx
            bs_atoms = [self.atoms[idx] for idx in [i for i in self.atoms.keys()
                                                    if self.atoms[i].OBAtom.GetResidue().GetIdx() in bs_res]
                        if idx in self.idx_to_pdb_mapping and self.idx_to_pdb_mapping[idx] not in self.altconf]
            bs_atoms_refined = []

            # Create hash with BSRES -> (MINDIST_TO_LIG, AA_TYPE)
            # and refine binding site atom selection with exact threshold
            min_dist = {}
            for r in bs_atoms:
                bs_res_id = ''.join([str(whichresnumber(r)), whichchain(r)])
                for l in ligand.mol.atoms:
                    distance = euclidean3d(r.coords, l.coords)
                    if bs_res_id not in min_dist:
                        min_dist[bs_res_id] = (distance, whichrestype(r))
                    elif min_dist[bs_res_id][0] > distance:
                        min_dist[bs_res_id] = (distance, whichrestype(r))
                    if distance <= config.BS_DIST and r not in bs_atoms_refined:
                        bs_atoms_refined.append(r)
            num_bs_atoms = len(bs_atoms_refined)
            message('Binding site atoms in vicinity (%.1f A max. dist: %i).\n' % (config.BS_DIST, num_bs_atoms), indent=True)

            bs_obj = BindingSite(bs_atoms_refined, self.protcomplex, self, self.altconf, min_dist)
            pli_obj = PLInteraction(lig_obj, bs_obj, self)
            self.interaction_sets[ligand.mol.title] = pli_obj

    def extract_bs(self, cutoff, ligcentroid, resis):
        """Return list of ids from residues belonging to the binding site"""
        return [obres.GetIdx() for obres in resis if self.res_belongs_to_bs(obres, cutoff, ligcentroid)]

    def res_belongs_to_bs(self, res, cutoff, ligcentroid):
        """Check for each residue if its centroid is within a certain distance to the ligand centroid."""
        rescentroid = centroid([(atm.x(), atm.y(), atm.z()) for atm in pybel.ob.OBResidueAtomIter(res)])
        return True if euclidean3d(rescentroid, ligcentroid) < cutoff else False

    def get_atom(self, idx):
        return self.atoms[idx]

    @property
    def output_path(self):
        return self.output_path

    @output_path.setter
    def output_path(self, path):
        self.output_path = tilde_expansion(path)
