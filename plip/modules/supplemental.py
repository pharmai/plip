"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
supplemental.py - Supplemental functions for PLIP analysis.
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

# Compatibility
from __future__ import print_function

# PLIP Modules
import config

# Python standard library
import re
from collections import namedtuple
import os
if os.name != 'nt':  # Resource module not available for Windows
    import resource
import subprocess
import codecs
import gzip
import zipfile

# External libraries
import pybel
from pybel import *
from openbabel import *
import numpy as np
from pymol import cmd
from pymol import finish_launching


def is_lig(hetid):
    """Checks if a PDB compound can be excluded as a small molecule ligand"""
    h = hetid.upper()
    return not (h == 'HOH' or h in config.UNSUPPORTED)


def parse_pdb(fil):
    """Extracts additional information from PDB files.
    I. When reading in a PDB file, OpenBabel numbers ATOMS and HETATOMS continously.
    In PDB files, TER records are also counted, leading to a different numbering system.
    This functions reads in a PDB file and provides a mapping as a dictionary.
    II. Additionally, it returns a list of modified residues.
    III. Furthermore, covalent linkages between ligands and protein residues/other ligands are identified
    IV. Alternative conformations
    """
    # #@todo Also consider SSBOND entries here
    i, j = 0, 0  # idx and PDB numbering
    d = {}
    modres = set()
    covlinkage = namedtuple("covlinkage", "id1 chain1 pos1 conf1 id2 chain2 pos2 conf2")
    covalent = []
    alt = []
    previous_ter = False
    for line in fil:
        if line.startswith(("ATOM", "HETATM")):

            # Retrieve alternate conformations
            atomid, location = int(line[6:11]), line[16]
            location = 'A' if location == ' ' else location
            if location != 'A':
                alt.append(atomid)

            if not previous_ter:
                i += 1
                j += 1
            else:
                i += 1
                j += 2
            d[i] = j
            previous_ter = False
        # Numbering Changes at TER records
        if line.startswith("TER"):
            previous_ter = True
        # Get modified residues
        if line.startswith("MODRES"):
            modres.add(line[12:15].strip())
        # Get covalent linkages between ligands
        if line.startswith("LINK"):
            conf1, id1, chain1, pos1 = line[16].strip(), line[17:20].strip(), line[21].strip(), int(line[22:26])
            conf2, id2, chain2, pos2 = line[46].strip(), line[47:50].strip(), line[51].strip(), int(line[52:56])
            covalent.append(covlinkage(id1=id1, chain1=chain1, pos1=pos1, conf1=conf1,
                                       id2=id2, chain2=chain2, pos2=pos2, conf2=conf2))
    return d, modres, covalent, alt


def extract_pdbid(string):
    """Use regular expressions to get a PDB ID from a string"""
    p = re.compile("[0-9][0-9a-z]{3}")
    m = p.search(string.lower())
    try:
        return m.group()
    except AttributeError:
        return "UnknownProtein"


def whichrestype(atom):
    """Returns the residue name of an Pybel or OpenBabel atom."""
    if isinstance(atom, Atom):
        try:
            return atom.OBAtom.GetResidue().GetName() if atom.OBAtom.GetResidue() is not None else None
        except AttributeError:
            return None
    elif isinstance(atom, OBAtom):
            return atom.GetResidue().GetName() if atom.GetResidue() is not None else None
    else:
        return None


def whichresnumber(atom):
    """Returns the residue number of an Pybel or OpenBabel atom (numbering as in original PDB file)."""
    if isinstance(atom, Atom):
        return atom.OBAtom.GetResidue().GetNum() if atom.OBAtom.GetResidue() is not None else None
    elif isinstance(atom, OBAtom):
        return atom.GetResidue().GetNum() if atom.GetResidue() is not None else None
    else:
        return None


def whichchain(atom):
    """Returns the residue number of an PyBel or OpenBabel atom."""
    if isinstance(atom, Atom):
        return atom.OBAtom.GetResidue().GetChain() if atom.OBAtom.GetResidue() is not None else None
    elif isinstance(atom, OBAtom):
        return atom.GetResidue().GetChain() if atom.GetResidue() is not None else None
    else:
        return None

#########################
# Mathematical operations
#########################


def euclidean3d(v1, v2):
    """Faster implementation of euclidean distance for the 3D case."""
    if not len(v1) == 3 and len(v2) == 3:
        print("Vectors are not in 3D space. Returning None.")
        return None
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def vector(p1, p2):
    """Vector from p1 to p2.
    :param p1: coordinates of point p1
    :param p2: coordinates of point p2
    :returns : numpy array with vector coordinates
    """
    return None if len(p1) != len(p2) else np.array([p2[i] - p1[i] for i in xrange(len(p1))])


def vecangle(v1, v2, deg=True):
    """Calculate the angle between two vectors
    :param v1: coordinates of vector v1
    :param v2: coordinates of vector v2
    :returns : angle in degree or rad
    """
    if np.array_equal(v1, v2):
        return 0.0
    dm = np.dot(v1, v2)
    cm = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle = np.arccos(round(dm / cm, 10))  # Round here to prevent floating point errors
    return np.degrees([angle, ])[0] if deg else angle


def normalize_vector(v):
    """Take a vector and return the normalized vector
    :param v: a vector v
    :returns : normalized vector v
    """
    norm = np.linalg.norm(v)
    return v/norm if not norm == 0 else v


def centroid(coo):
    """Calculates the centroid from a 3D point cloud and returns the coordinates
    :param coo: Array of coordinate arrays
    :returns : centroid coordinates as list
    """
    return map(np.mean, (([c[0] for c in coo]), ([c[1] for c in coo]), ([c[2] for c in coo])))


def projection(pnormal1, ppoint, tpoint):
    """Calculates the centroid from a 3D point cloud and returns the coordinates
    :param pnormal1: normal of plane
    :param ppoint: coordinates of point in the plane
    :param tpoint: coordinates of point to be projected
    :returns : coordinates of point orthogonally projected on the plane
    """
    # Choose the plane normal pointing to the point to be projected
    pnormal2 = [coo*(-1) for coo in pnormal1]
    d1 = euclidean3d(tpoint, pnormal1 + ppoint)
    d2 = euclidean3d(tpoint, pnormal2 + ppoint)
    pnormal = pnormal1 if d1 < d2 else pnormal2
    # Calculate the projection of tpoint to the plane
    sn = -np.dot(pnormal, vector(ppoint, tpoint))
    sd = np.dot(pnormal, pnormal)
    sb = sn / sd
    return [c1 + c2 for c1, c2 in zip(tpoint, [sb * pn for pn in pnormal])]


def cluster_doubles(double_list):
    """Given a list of doubles, they are clustered if they share one element
    :param double_list: list of doubles
    :returns : list of clusters (tuples)
    """
    location = {}  # hashtable of which cluster each element is in
    clusters = []
    # Go through each double
    for t in double_list:
        a, b = t[0], t[1]
        # If they both are already in different clusters, merge the clusters
        if a in location and b in location:
            if location[a] != location[b]:
                if location[a] < location[b]:
                    clusters[location[a]] = clusters[location[a]].union(clusters[location[b]])  # Merge clusters
                    clusters = clusters[:location[b]] + clusters[location[b]+1:]
                else:
                    clusters[location[b]] = clusters[location[b]].union(clusters[location[a]])  # Merge clusters
                    clusters = clusters[:location[a]] + clusters[location[a]+1:]
                # Rebuild index of locations for each element as they have changed now
                location = {}
                for i, cluster in enumerate(clusters):
                    for c in cluster:
                        location[c] = i
        else:
            # If a is already in a cluster, add b to that cluster
            if a in location:
                clusters[location[a]].add(b)
                location[b] = location[a]
            # If b is already in a cluster, add a to that cluster
            if b in location:
                clusters[location[b]].add(a)
                location[a] = location[b]
            # If neither a nor b is in any cluster, create a new one with a and b
            if not (b in location and a in location):
                clusters.append(set(t))
                location[a] = len(clusters) - 1
                location[b] = len(clusters) - 1
    return map(tuple, clusters)


#################
# File operations
#################

def tilde_expansion(folder_path):
    """Tilde expansion, i.e. converts '~' in paths into <value of $HOME>."""
    return os.path.expanduser(folder_path) if '~' in folder_path else folder_path


def folder_exists(folder_path):
    """Checks if a folder exists"""
    return os.path.exists(folder_path)


def create_folder_if_not_exists(folder_path):
    """Creates a folder if it does not exists."""
    folder_path = tilde_expansion(folder_path)
    folder_path = "".join([folder_path, '/']) if not folder_path[-1] == '/' else folder_path
    direc = os.path.dirname(folder_path)
    if not folder_exists(direc):
        os.makedirs(direc)


def cmd_exists(c):
    return subprocess.call("type " + c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

################
# PyMOL-specific
################


def object_exists(object_name):
    """Checks if an object exists in the open PyMOL session."""
    return object_name in cmd.get_names("objects")


def initialize_pymol(options):
    """Initializes PyMOL"""
    # Pass standard arguments of function to prevent PyMOL from printing out PDB headers (workaround)
    finish_launching(args=['pymol', options, '-K'])
    cmd.reinitialize()


def start_pymol(quiet=False, options='-p', run=False):
    """Starts up PyMOL and sets general options. Quiet mode suppresses all PyMOL output.
    Command line options can be passed as the second argument."""
    import pymol
    pymol.pymol_argv = ['pymol', '%s' % options] + sys.argv[1:]
    if run:
        initialize_pymol(options)
    if quiet:
        cmd.feedback('disable', 'all', 'everything')


def standard_settings():
    """Sets up standard settings for a nice visualization."""
    cmd.set('bg_rgb', [1.0, 1.0, 1.0])  # White background
    cmd.set('depth_cue', 0)  # Turn off depth cueing (no fog)
    cmd.set('cartoon_side_chain_helper', 1)  # Improve combined visualization of sticks and cartoon
    cmd.set('cartoon_fancy_helices', 1)  # Nicer visualization of helices (using tapered ends)
    cmd.set('transparency_mode', 1)  # Turn on multilayer transparency
    cmd.set('dash_radius', 0.05)
    set_custom_colorset()


def set_custom_colorset():
    """Defines a colorset with matching colors. Provided by Joachim."""
    cmd.set_color('myorange', '[253, 174, 97]')
    cmd.set_color('mygreen', '[171, 221, 164]')
    cmd.set_color('myred', '[215, 25, 28]')
    cmd.set_color('myblue', '[43, 131, 186]')
    cmd.set_color('mylightblue', '[158, 202, 225]')
    cmd.set_color('mylightgreen', '[229, 245, 224]')


#############################################
# Following code adapted from Joachim Haupt #
#############################################


def getligs(mol, altconf, idx_to_pdb, modres, covalent):
    """Get all ligands from a PDB file. Adapted from Joachim's structTools"""
    #############################
    # Read in file and get name #
    #############################

    data = namedtuple('ligand', 'mol mapping water members longname type')
    ligands = []

    #########################
    # Filtering using lists #
    #########################

    all_res1 = [o for o in pybel.ob.OBResidueIter(mol.OBMol)
               if not (o.GetResidueProperty(9) or o.GetResidueProperty(0))]
    all_lignames = set([a.GetName() for a in all_res1])

    water = [o for o in pybel.ob.OBResidueIter(mol.OBMol) if o.GetResidueProperty(9)]
    all_res2 = [a for a in all_res1 if is_lig(a.GetName()) and a.GetName() not in modres]  # Filter out non-ligands

    ############################################
    # Filtering by counting and artifacts list #
    ############################################
    artifacts = []
    unique_ligs = set(a.GetName() for a in all_res2)
    for ulig in unique_ligs:
        # Discard if appearing 15 times or more and is possible artifact
        if ulig in config.biolip_list and [a.GetName() for a in all_res2].count(ulig) >= 15:
            artifacts.append(ulig)
    all_res3 = [a for a in all_res2 if a.GetName() not in artifacts]
    all_res_dict = {(a.GetName(), a.GetChain(), a.GetNum()): a for a in all_res3}

    #######################################
    # Basic support for RNA/DNA as ligand #
    #######################################
    nucleotides = ['A', 'C', 'T', 'G', 'U', 'DA', 'DC', 'DT', 'DG', 'DU']
    dna_rna = {}  # Dictionary of DNA/RNA residues by chain
    covlinkage = namedtuple("covlinkage", "id1 chain1 pos1 conf1 id2 chain2 pos2 conf2")
    # Create missing covlinkage entries for DNA/RNA
    for ligand in all_res_dict:
        resname, chain, pos = ligand
        if resname in nucleotides:
            if chain not in dna_rna:
                dna_rna[chain] = [(resname, pos), ]
            else:
                dna_rna[chain].append((resname, pos))
    for chain in dna_rna:
        nuc_list = dna_rna[chain]
        for i, nucleotide in enumerate(nuc_list):
            if not i == len(nuc_list) - 1:
                name, pos = nucleotide
                nextnucleotide = nuc_list[i + 1]
                nextname, nextpos = nextnucleotide
                newlink = covlinkage(id1=name, chain1=chain, pos1=pos, conf1='',
                                     id2=nextname, chain2=chain, pos2=nextpos, conf2='')
                covalent.append(newlink)

    #########################
    # Identify kmer ligands #
    #########################

    lignames = list(set([a.GetName() for a in all_res3]))
    # Remove all those not considered by ligands and pairings including alternate conformations

    ligdoubles = [[(link.id1, link.chain1, link.pos1),
                   (link.id2, link.chain2, link.pos2)] for link in
                  [c for c in covalent if c.id1 in lignames and c.id2 in lignames and
                   c.conf1 in ['A', ''] and c.conf2 in ['A', '']
                  and (c.id1, c.chain1, c.pos1) in all_res_dict and (c.id2, c.chain2, c.pos2) in all_res_dict]]
    kmers = cluster_doubles(ligdoubles)
    if not kmers:  # No ligand kmers, just normal independent ligands
        res_kmers = [[all_res_dict[res]] for res in all_res_dict]
    else:
        # res_kmers contains clusters of covalently bound ligand residues (kmer ligands)
        res_kmers = [[all_res_dict[res] for res in kmer] for kmer in kmers]

        # In this case, add other ligands which are not part of a kmer
        in_kmer = []
        for res_kmer in res_kmers:
            for res in res_kmer:
                in_kmer.append((res.GetName(), res.GetChain(), res.GetNum()))
        for res in all_res_dict:
            if res not in in_kmer:
                newres = [all_res_dict[res], ]
                res_kmers.append(newres)

    ###################
    # Extract ligands #
    ###################

    for kmer in res_kmers:  # iterate over all ligands
        members = [(res.GetName(), res.GetChain(), res.GetNum()) for res in kmer]
        rname, rchain, rnum = sorted(members)[0]  # representative name, chain, and number
        ordered_members = sorted(members, key=lambda x: (x[1], x[2]))
        names = [x[0] for x in ordered_members]
        longname = '-'.join([x[0] for x in ordered_members])
        if len(names) > 3:  # Polymer
            if len({'U', 'A', 'C', 'G'}.intersection(set(names))) != 0:
                ligtype = 'RNA'
            elif len({'DT', 'DA', 'DC', 'DG'}.intersection(set(names))) != 0:
                ligtype = 'DNA'
            else:
                ligtype = "POLYMER"

        else:
            ligtype = 'SMALLMOLECULE'

        for name in names:
            if name in config.METAL_IONS:
                if len(names) == 1:
                    ligtype = 'ION'
                else:
                    if "ION" not in ligtype:
                        ligtype += '+ION'
        hetatoms = set()
        for obresidue in kmer:
            hetatoms_res = set([(obatom.GetIdx(), obatom) for obatom in pybel.ob.OBResidueAtomIter(obresidue)
                        if not obatom.IsHydrogen()])

            hetatoms_res = set([atm for atm in hetatoms_res if not idx_to_pdb[atm[0]] in altconf])  # Remove alt. conformations
            hetatoms.update(hetatoms_res)
        if len(hetatoms) == 0:
            continue
        hetatoms = dict(hetatoms)  # make it a dict with idx as key and OBAtom as value
        lig = pybel.ob.OBMol()  # new ligand mol
        neighbours = dict()
        for obatom in hetatoms.values():  # iterate over atom objects
            idx = obatom.GetIdx()
            lig.AddAtom(obatom)
            # ids of all neighbours of obatom
            neighbours[idx] = set([neighbour_atom.GetIdx() for neighbour_atom
                                   in pybel.ob.OBAtomAtomIter(obatom)]) & set(hetatoms.keys())

        ##############################################################
        # map the old atom idx of OBMol to the new idx of the ligand #
        ##############################################################

        newidx = dict(zip(hetatoms.keys(), [obatom.GetIdx() for obatom in pybel.ob.OBMolAtomIter(lig)]))
        mapold = dict(zip(newidx.values(), newidx))
        # copy the bonds
        for obatom in hetatoms:
            for neighbour_atom in neighbours[obatom]:
                bond = hetatoms[obatom].GetBond(hetatoms[neighbour_atom])
                lig.AddBond(newidx[obatom], newidx[neighbour_atom], bond.GetBondOrder())
        lig = pybel.Molecule(lig)
        # For kmers, the representative ids are chosen (first residue of kmer)
        lig.data.update({'Name': rname,
                         'Chain': rchain,
                         'ResNr': rnum})
        lig.title = '-'.join((rname, rchain, str(rnum)))
        ligands.append(data(mol=lig, mapping=mapold, water=water, members=members, longname=longname, type=ligtype))
    excluded = sorted(list(all_lignames.difference(set(lignames))))
    return ligands, excluded


def read_pdb(pdbfname):
    """Reads a given PDB file and returns a Pybel Molecule."""
    pybel.ob.obErrorLog.StopLogging()  # Suppress all OpenBabel warnings
    if os.name != 'nt':  # Resource module not available for Windows
        maxsize = resource.getrlimit(resource.RLIMIT_STACK)[-1]
        resource.setrlimit(resource.RLIMIT_STACK, (min(2 ** 28, maxsize), maxsize))
    sys.setrecursionlimit(10 ** 5)  # increase Python recoursion limit
    return readmol(pdbfname)


def read(fil):
    """Returns a file handler and detects gzipped files."""
    if os.path.splitext(fil)[-1] == '.gz':
        return gzip.open(fil, 'rb')
    elif os.path.splitext(fil)[-1] == '.zip':
        zf = zipfile.ZipFile(fil, 'r')
        return zf.open(zf.infolist()[0].filename)
    else:
        try:
            codecs.open(fil, 'r', 'utf-8').read()
            return codecs.open(fil, 'r', 'utf-8')
        except UnicodeDecodeError:
            return open(fil, 'r')


def readmol(path):
    """Reads the given molecule file and returns the corresponding Pybel molecule as well as the input file type.
    In contrast to the standard Pybel implementation, the file is closed properly."""
    supported_formats = ['pdb', 'pdbqt']
    obc = pybel.ob.OBConversion()

    with read(path) as f:
        filestr = str(f.read())

    for sformat in supported_formats:
        obc.SetInFormat(sformat)
        mol = pybel.ob.OBMol()
        obc.ReadString(mol, filestr)
        if not mol.Empty():
            if sformat == 'pdbqt':
                message('[EXPERIMENTAL] Input is PDBQT file. Some features (especially visualization) might not '
                        'work as expected. Please consider using PDB format instead.\n')
            return pybel.Molecule(mol), sformat
    sysexit(4, 'No valid PDB or PDBQT file provided.')


def sysexit(code, msg):
    """Exit using an custom error message and error code."""
    sys.stderr.write(msg)
    sys.exit(code)


def message(msg, indent=False):
    """Writes messages in verbose mode"""
    if config.VERBOSE:
        if indent:
            msg = '  ' + msg
        sys.stdout.write(msg)
