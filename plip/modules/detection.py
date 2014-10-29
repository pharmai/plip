"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
Module for detection of non-covalent interactions.
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


# Own modules
import itertools

from supplemental import *


######################
#THRESHOLD DEFINITIONS
######################

# Some distance thresholds were extended (max. 1.0A) if too restrictive too account for low-quality structures
HYDROPH_DIST_MAX = 4.0  # Distance cutoff for detection of hydrophobic contacts
HBOND_DIST_MAX = 4.0  # Max. distance between hydrogen bond donor and acceptor (Hubbard & Haider, 2001) + 0.5 A
HBOND_DON_ANGLE_MIN = 90  # Min. angle at the hydrogen bond donor (Hubbard & Haider, 2001)
PISTACK_DIST_MAX = 7.5  # Max. distance for parallel or offset pistacking (McGaughey, 1998)
PISTACK_ANG_DEV = 30  # Max. Deviation from parallel or perpendicular orientation (in degrees)
PISTACK_OFFSET_MAX = 2.0  # Maximum offset of the two rings (corresponds to the radius of benzene + 0.5 A)
PICATION_DIST_MAX = 6.0  # Max. distance between charged atom and aromatic ring center (Gallivan and Dougherty, 1999)
SALTBRIDGE_DIST_MAX = 5.0  # Max. distance between centers of charge for salt bridges (Barlow and Thornton, 1983) + 1.0
HALOGEN_DIST_MAX = 4.0  # Max. distance between oxy. and halogen (Halogen bonds in biological molecules., Auffinger)+0.5
HALOGEN_ACC_ANGLE = 120  # Optimal acceptor angle (Halogen bonds in biological molecules., Auffinger)
HALOGEN_DON_ANGLE = 165  # Optimal donor angle (Halogen bonds in biological molecules., Auffinger)
HALOGEN_ANGLE_DEV = 30  # Max. deviation from optimal angle
WATER_BRIDGE_MINDIST = 2.5  # Min. distance between water oxygen and polar atom (Jiang et al., 2005) -0.1
WATER_BRIDGE_MAXDIST = 4.0  # Max. distance between water oxygen and polar atom (Jiang et al., 2005) +0.4
WATER_BRIDGE_OMEGA_MIN = 80  # Min. angle between acceptor, water oxygen and donor hydrogen (Jiang et al., 2005)
WATER_BRIDGE_OMEGA_MAX = 160  # Max. angle between acceptor, water oxygen and donor hydrogen (Jiang et al., 2005)
WATER_BRIDGE_THETA_MIN = 100  # Min. angle between water oxygen, donor hydrogen and donor atom (Jiang et al., 2005)


##################################################
# FUNCTIONS FOR DETECTION OF SPECIFIC INTERACTIONS
##################################################
#@todo Implement checks to avoid atoms of the ligand or protein being in the path of an interaction. (e.g. 1f0r)

def hydrophobic_interactions(atom_set_a, atom_set_b):
    """Detection of hydrophobic pliprofiler between atom_set_a (binding site) and atom_set_b (ligand).
    Definition: All pairs of qualified carbon atoms within a distance of HYDROPH_DIST_MAX
    """
    data = namedtuple('hydroph_interaction', 'bsatom ligatom distance restype resnr')
    i_set = []
    for a, b in itertools.product(atom_set_a.atoms, atom_set_b.atoms):
        e = euclidean3d(a.coords, b.coords)
        if e < HYDROPH_DIST_MAX:
            i_set.append(data(bsatom=a, ligatom=b, distance=e, restype=whichrestype(a), resnr=whichresnumber(a)))
    return i_set


def hbonds(acceptors, donor_pairs, protisdon, typ):
    """Detection of hydrogen bonds between sets of acceptors and donor pairs.
    Definition: All pairs of hydrogen bond acceptor and donors with
    donor hydrogens and acceptor showing a distance within HBOND DIST MIN and HBOND DIST MAX
    and donor angles above HBOND_DON_ANGLE_MIN
    """
    i_set = []
    data = namedtuple('hbond', 'a d h distance_ah distance_ad angle type protisdon resnr restype')
    for acc, don in itertools.product(acceptors, donor_pairs):
        if typ == 'strong':  # Regular (strong) hydrogen bonds
            dist_ah = euclidean3d(acc.a.coords, don.h.coords)
            dist_ad = euclidean3d(acc.a.coords, don.d.coords)
            if dist_ad < HBOND_DIST_MAX:
                vec1, vec2 = vector(don.h.coords, don.d.coords), vector(don.h.coords, acc.a.coords)
                v = vecangle(vec1, vec2)
                if v > HBOND_DON_ANGLE_MIN:
                    restype = whichrestype(don.d) if protisdon else whichrestype(acc.a)
                    resnr = whichresnumber(don.d)if protisdon else whichresnumber(acc.a)
                    i_set.append(data(a=acc.a, d=don.d, h=don.h, distance_ah=dist_ah, distance_ad=dist_ad, angle=v,
                                      type=typ, protisdon=protisdon, resnr=resnr, restype=restype))
    return i_set


def pistacking(rings_bs, rings_lig):
    """Return all pi-stackings between the given aromatic ring systems in receptor and ligand."""
    i_set = []
    data = namedtuple('pistack', 'proteinring ligandring distance angle offset type restype resnr')
    for r, l in itertools.product(rings_bs, rings_lig):
        # DISTANCE AND RING ANGLE CALCULATION
        d = euclidean3d(r.center, l.center)
        b = vecangle(r.normal, l.normal)
        a = min(b, 180-b if not 180-b < 0 else b)  # Smallest of two angles, depending on direction of normal

        # RING CENTER OFFSET CALCULATION (project each ring center into the other ring)
        proj1 = projection(l.normal, l.center, r.center)
        proj2 = projection(r.normal, r.center, l.center)
        offset = min(euclidean3d(proj1, l.center), euclidean3d(proj2, r.center))

        # RECEPTOR DATA
        resnr, restype = whichresnumber(r.atoms[0]), whichrestype(r.atoms[0])

        # SELECTION BY DISTANCE, ANGLE AND OFFSET
        if d < PISTACK_DIST_MAX:
            if 0 < a < PISTACK_ANG_DEV and offset < PISTACK_OFFSET_MAX:
                i_set.append(data(proteinring=r, ligandring=l, distance=d, angle=a, offset=offset,
                                  type='P', resnr=resnr, restype=restype))
            if 90-PISTACK_ANG_DEV < a < 90+PISTACK_ANG_DEV and offset < PISTACK_OFFSET_MAX:
                i_set.append(data(proteinring=r, ligandring=l, distance=d, angle=a, offset=offset,
                                  type='T', resnr=resnr, restype=restype))
    return i_set


def pication(rings, pos_charged, protcharged):
    """Return all pi-Cation interaction between aromatic rings and positively charged groups.
    For tertiary and quaternary amines, check also the angle between the ring and the nitrogen.
    """
    i_set = []
    data = namedtuple('pication', 'ring charge distance offset type restype resnr protcharged')
    if not len(rings) == 0 and not len(pos_charged) == 0:
        for ring in rings:
            c = ring.center
            for p in pos_charged:
                d = euclidean3d(c, p.center)
                # Project the center of charge into the ring and measure distance to ring center
                proj = projection(ring.normal, ring.center, p.center)
                offset = euclidean3d(proj, ring.center)
                if d < PICATION_DIST_MAX and offset < PISTACK_OFFSET_MAX:
                    if type(p).__name__ == 'lcharge' and p.fgroup == 'tertamine':
                        # Special case here if the ligand has a tertiary amine, check an additional angle
                        # Otherwise, we might have have a pi-cation interaction 'through' the ligand
                        n_atoms = [a_neighbor for a_neighbor in OBAtomAtomIter(p.atoms[0].OBAtom)]
                        n_atoms_coords = [(a.x(), a.y(), a.z()) for a in n_atoms]
                        amine_normal = np.cross(vector(n_atoms_coords[0], n_atoms_coords[1]),
                                                vector(n_atoms_coords[2], n_atoms_coords[0]))
                        b = vecangle(ring.normal, amine_normal)
                        # Smallest of two angles, depending on direction of normal
                        a = min(b, 180-b if not 180-b < 0 else b)
                        if not a > 30.0:
                            resnr, restype = whichresnumber(ring.atoms[0]), whichrestype(ring.atoms[0])
                            i_set.append(data(ring=ring, charge=p, distance=d,
                                              offset=offset, type='regular', restype=restype,
                                              resnr=resnr, protcharged=protcharged))
                        break
                    resnr = whichresnumber(p.atoms[0]) if protcharged else whichresnumber(ring.atoms[0])
                    restype = whichrestype(p.atoms[0]) if protcharged else whichrestype(ring.atoms[0])
                    i_set.append(data(ring=ring, charge=p, distance=d,
                                      offset=offset, type='regular', restype=restype,
                                      resnr=resnr, protcharged=protcharged))
    return i_set


def saltbridge(poscenter, negcenter, protispos):
    """Detect all salt bridges (pliprofiler between centers of positive and negative charge)"""
    data = namedtuple('saltbridge', 'positive negative distance protispos resnr restype')
    i_set = []
    for pc, nc in itertools.product(poscenter, negcenter):
        dists = []
        for pa, na in itertools.product(pc.atoms, [n for n in nc.atoms if not nc.atoms == []]):
            dists.append(euclidean3d(pa.coords, na.coords))
        if min(dists) < SALTBRIDGE_DIST_MAX:
            resnr = pc.resnr if protispos else nc.resnr
            restype = pc.restype if protispos else nc.restype
            i_set.append(data(positive=pc, negative=nc, distance=euclidean3d(pc.center, nc.center),
                              protispos=protispos, resnr=resnr, restype=restype))
    return i_set


def halogen(acceptor, donor):
    """Detect all halogen bonds of the type Y-O...X-C"""
    data = namedtuple('halogenbond', 'acc don distance don_angle acc_angle restype resnr donortype')
    i_set = []
    for acc, don in itertools.product(acceptor, donor):
        dist = euclidean3d(acc.o.coords, don.x.coords)
        if dist < HALOGEN_DIST_MAX:
            vec1, vec2 = vector(acc.o.coords, acc.y.coords), vector(acc.o.coords, don.x.coords)
            vec3, vec4 = vector(don.x.coords, acc.o.coords), vector(don.x.coords, don.c.coords)
            acc_angle, don_angle = vecangle(vec1, vec2), vecangle(vec3, vec4)
            if HALOGEN_ACC_ANGLE-HALOGEN_ANGLE_DEV < acc_angle < HALOGEN_ACC_ANGLE+HALOGEN_ANGLE_DEV:
                if HALOGEN_DON_ANGLE-HALOGEN_ANGLE_DEV < don_angle < HALOGEN_DON_ANGLE+HALOGEN_ANGLE_DEV:
                    i_set.append(data(acc=acc, don=don, distance=dist, don_angle=don_angle, acc_angle=acc_angle,
                                      restype=whichrestype(acc.o), resnr=whichresnumber(acc.o),
                                      donortype=don.x.OBAtom.GetType()))
    return i_set


def water_bridges(bs_hba, lig_hba, bs_hbd, lig_hbd, water):
    """Find water-bridged hydrogen bonds between ligand and protein. For now only considers bridged of first degree."""
    i_set = []
    data = namedtuple('waterbridge', 'a d h water distance_aw distance_dw d_angle w_angle type resnr restype protisdon')
    # First find all acceptor-water pairs with distance within d
    # and all donor-water pairs with distance within d and angle greater theta
    lig_aw, prot_aw, lig_dw, prot_hw = [], [], [], []
    for w in water:
        for acc1 in lig_hba:
            dist = euclidean3d(acc1.a.coords, w.coords)
            if WATER_BRIDGE_MINDIST <= dist <= WATER_BRIDGE_MAXDIST:
                lig_aw.append((acc1, w, dist))
        for acc2 in bs_hba:
            dist = euclidean3d(acc2.a.coords, w.coords)
            if WATER_BRIDGE_MINDIST <= dist <= WATER_BRIDGE_MAXDIST:
                prot_aw.append((acc2, w, dist))
        for don1 in lig_hbd:
            dist = euclidean3d(don1.d.coords, w.coords)
            d_angle = vecangle(vector(don1.h.coords, don1.d.coords), vector(don1.h.coords, w.coords))
            if WATER_BRIDGE_MINDIST <= dist <= WATER_BRIDGE_MAXDIST and d_angle > WATER_BRIDGE_THETA_MIN:
                lig_dw.append((don1, w, dist, d_angle))
        for don2 in bs_hbd:
            dist = euclidean3d(don2.d.coords, w.coords)
            d_angle = vecangle(vector(don2.h.coords, don2.d.coords), vector(don2.h.coords, w.coords))
            if WATER_BRIDGE_MINDIST <= dist <= WATER_BRIDGE_MAXDIST and d_angle > WATER_BRIDGE_THETA_MIN:
                prot_hw.append((don2, w, dist, d_angle))

    for l, p in itertools.product(lig_aw, prot_hw):
        acc, wl, distance_aw = l
        don, wd, distance_dw, d_angle = p
        if wl == wd:  # Same water molecule and angle within omega
            w_angle = vecangle(vector(acc.a.coords, wl.coords), vector(wl.coords, don.h.coords))
            if WATER_BRIDGE_OMEGA_MIN < w_angle < WATER_BRIDGE_OMEGA_MAX:

                i_set.append(data(a=acc.a, d=don.d, h=don.h, water=wl, distance_aw=distance_aw, distance_dw=distance_dw,
                                  d_angle=d_angle, w_angle=w_angle, type='first_deg', resnr=whichresnumber(don.d),
                                  restype=whichrestype(don.d), protisdon=True))

    for p, l in itertools.product(prot_aw, lig_dw):
        acc, wl, distance_aw = p
        don, wd, distance_dw, d_angle = l
        if wl == wd:  # Same water molecule and angle within omega
            w_angle = vecangle(vector(acc.a.coords, wl.coords), vector(wl.coords, don.h.coords))
            if WATER_BRIDGE_OMEGA_MIN < w_angle < WATER_BRIDGE_OMEGA_MAX:
                i_set.append(data(a=acc.a, d=don.d, h=don.h, water=wl, distance_aw=distance_aw, distance_dw=distance_dw,
                                  d_angle=d_angle, w_angle=w_angle, type='first_deg', resnr=whichresnumber(acc.a),
                                  restype=whichrestype(acc.a), protisdon=False))

    return i_set
