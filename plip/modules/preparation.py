"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
preparation.py - Prepare PDB input files for processing.
"""


# Python Standard Library
from __future__ import absolute_import
from builtins import filter
from operator import itemgetter

# Own modules
from .detection import *
from .supplemental import *
from . import config


################
# MAIN CLASSES #
################

class PDBParser:
    def __init__(self, pdbpath, as_string):
        self.as_string = as_string
        self.pdbpath = pdbpath
        self.num_fixed_lines = 0
        self.covlinkage = namedtuple("covlinkage", "id1 chain1 pos1 conf1 id2 chain2 pos2 conf2")
        self.proteinmap, self.modres, self.covalent, self.altconformations, self.corrected_pdb = self.parse_pdb()


    def parse_pdb(self):
        """Extracts additional information from PDB files.
        I. When reading in a PDB file, OpenBabel numbers ATOMS and HETATOMS continously.
        In PDB files, TER records are also counted, leading to a different numbering system.
        This functions reads in a PDB file and provides a mapping as a dictionary.
        II. Additionally, it returns a list of modified residues.
        III. Furthermore, covalent linkages between ligands and protein residues/other ligands are identified
        IV. Alternative conformations
        """
        if self.as_string:
            fil = self.pdbpath.rstrip('\n').split('\n') # Removing trailing newline character
        else:
            f = read(self.pdbpath)
            fil = f.readlines()
            f.close()
        corrected_lines = []
        i, j = 0, 0  # idx and PDB numbering
        d = {}
        modres = set()
        covalent = []
        alt = []
        previous_ter = False

        # Standard without fixing
        if not config.NOFIX:
            if not config.PLUGIN_MODE:
                lastnum = 0 # Atom numbering (has to be consecutive)
                other_models = False
                for line in fil:
                    if not other_models: # Only consider the first model in an NRM structure
                        corrected_line, newnum = self.fix_pdbline(line, lastnum)
                        if corrected_line is not None:
                            if corrected_line.startswith('MODEL'):
                                try: # Get number of MODEL (1,2,3)
                                    model_num = int(corrected_line[10:14])
                                    if model_num > 1: # MODEL 2,3,4 etc.
                                        other_models = True
                                except ValueError:
                                    write_message("Ignoring invalid MODEL entry: %s\n" % corrected_line, mtype='debug')
                            corrected_lines.append(corrected_line)
                            lastnum = newnum
                corrected_pdb = ''.join(corrected_lines)
            else:
                corrected_pdb = self.pdbpath
                corrected_lines = fil
        else:
            corrected_pdb = self.pdbpath
            corrected_lines = fil

        for line in corrected_lines:
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
                covalent.append(self.get_linkage(line))
        return d, modres, covalent, alt, corrected_pdb

    def fix_pdbline(self, pdbline, lastnum):
        """Fix a PDB line if information is missing."""
        pdbqt_conversion = {
        "HD": "H", "HS": "H", "NA": "N",
        "NS": "N", "OA": "O", "OS": "O", "SA": "S"}
        fixed = False
        newnum = 0
        forbidden_characters = "[^a-zA-Z0-9_]"
        pdbline = pdbline.strip('\n')
        # Some MD / Docking tools produce empty lines, leading to segfaults
        if len(pdbline.strip()) == 0:
            self.num_fixed_lines += 1
            return None, lastnum
        if len(pdbline) > 100: # Should be 80 long
            self.num_fixed_lines += 1
            return None, lastnum
        # TER Entries also have continuing numbering, consider them as well
        if pdbline.startswith('TER'):
            newnum = lastnum + 1
        if pdbline.startswith('ATOM'):
            newnum = lastnum + 1
            currentnum = int(pdbline[6:11])
            resnum = pdbline[22:27].strip()
            resname = pdbline[17:21].strip()
            # Invalid residue number
            try:
                int(resnum)
            except ValueError:
                pdbline = pdbline[:22] + '   0 ' + pdbline[27:]
                fixed = True
            # Invalid characters in residue name
            if re.match(forbidden_characters, resname.strip()):
                pdbline = pdbline[:17] + 'UNK ' + pdbline[21:]
                fixed = True
            if lastnum + 1 != currentnum:
                pdbline = pdbline[:6] + (5 - len(str(newnum))) * ' ' + str(newnum) + ' ' + pdbline[12:]
                fixed = True
            # No chain assigned
            if pdbline[21] == ' ':
                pdbline = pdbline[:21] + 'A' + pdbline[22:]
                fixed = True
            if pdbline.endswith('H'):
                self.num_fixed_lines += 1
                return None, lastnum
            # Sometimes, converted PDB structures contain PDBQT atom types. Fix that.
            for pdbqttype in pdbqt_conversion:
                if pdbline.strip().endswith(pdbqttype):
                    pdbline = pdbline.strip()[:-2] + ' ' + pdbqt_conversion[pdbqttype] + '\n'
                    self.num_fixed_lines += 1
        if pdbline.startswith('HETATM'):
            newnum = lastnum + 1
            try:
                currentnum = int(pdbline[6:11])
            except ValueError:
                currentnum = None
                write_message("Invalid HETATM entry: %s\n" % pdbline, mtype='debug')
            if lastnum + 1 != currentnum:
                pdbline = pdbline[:6] + (5 - len(str(newnum))) * ' ' + str(newnum) + ' ' + pdbline[12:]
                fixed = True
            # No chain assigned or number assigned as chain
            if pdbline[21] == ' ':
                pdbline = pdbline[:21] + 'Z' + pdbline[22:]
                fixed = True
            # No residue number assigned
            if pdbline[23:26] == '   ':
                pdbline = pdbline[:23] + '999' + pdbline[26:]
                fixed = True
            # Non-standard Ligand Names
            ligname = pdbline[17:21].strip()
            if len(ligname) > 3:
                pdbline = pdbline[:17] + ligname[:3] + ' ' + pdbline[21:]
                fixed = True
            if re.match(forbidden_characters, ligname.strip()):
                pdbline = pdbline[:17] + 'LIG ' + pdbline[21:]
                fixed = True
            if len(ligname.strip()) == 0:
                pdbline = pdbline[:17] + 'LIG ' + pdbline[21:]
                fixed = True
            if pdbline.endswith('H'):
                self.num_fixed_lines += 1
                return None, lastnum
            # Sometimes, converted PDB structures contain PDBQT atom types. Fix that.
            for pdbqttype in pdbqt_conversion:
                if pdbline.strip().endswith(pdbqttype):
                    pdbline = pdbline.strip()[:-2] + ' ' + pdbqt_conversion[pdbqttype] + ' '
                    self.num_fixed_lines +=1
        self.num_fixed_lines += 1 if fixed else 0
        return pdbline + '\n', max(newnum, lastnum)

    def get_linkage(self, line):
        """Get the linkage information from a LINK entry PDB line."""
        conf1, id1, chain1, pos1 = line[16].strip(), line[17:20].strip(), line[21].strip(), int(line[22:26])
        conf2, id2, chain2, pos2 = line[46].strip(), line[47:50].strip(), line[51].strip(), int(line[52:56])
        return self.covlinkage(id1=id1, chain1=chain1, pos1=pos1, conf1=conf1,
                               id2=id2, chain2=chain2, pos2=pos2, conf2=conf2)


class LigandFinder:
    def __init__(self, proteincomplex, altconf, modres, covalent, mapper):
        self.lignames_all = None
        self.lignames_kept = None
        self.water = None
        self.proteincomplex = proteincomplex
        self.altconformations = altconf
        self.modresidues = modres
        self.covalent = covalent
        self.mapper = mapper
        self.ligands = self.getligs()
        self.excluded = sorted(list(self.lignames_all.difference(set(self.lignames_kept))))

    def getpeptides(self, chain):
        """If peptide ligand chains are defined via the command line options,
        try to extract the underlying ligand formed by all residues in the
        given chain without water
        """
        data = namedtuple('ligand', 'mol hetid chain position water members longname type atomorder can_to_pdb')
        all_from_chain = [o for o in pybel.ob.OBResidueIter(self.proteincomplex.OBMol) if o.GetChain() == chain] # All residues from chain
        if len(all_from_chain) == 0:
            return None
        else:
            water = [o for o in all_from_chain if o.GetResidueProperty(9)]
            non_water = water = [o for o in all_from_chain if not o.GetResidueProperty(9)]
            ligand = self.extract_ligand(non_water)
            return ligand


    def getligs(self):
        """Get all ligands from a PDB file and prepare them for analysis.
        Returns all non-empty ligands.
        """

        if config.PEPTIDES == [] and config.INTRA is None:
            # Extract small molecule ligands (default)
            ligands = []

            # Filter for ligands using lists
            ligand_residues, self.lignames_all, self.water = self.filter_for_ligands()

            all_res_dict = {(a.GetName(), a.GetChain(), a.GetNum()): a for a in ligand_residues}
            self.lignames_kept = list(set([a.GetName() for a in ligand_residues]))

            if not config.BREAKCOMPOSITE:
                #  Update register of covalent links with those between DNA/RNA subunits
                self.covalent += nucleotide_linkage(all_res_dict)
                #  Find fragment linked by covalent bonds
                res_kmers = self.identify_kmers(all_res_dict)
            else:
                res_kmers = [[a, ] for a in ligand_residues]
            write_message("{} ligand kmer(s) detected for closer inspection.\n".format(len(res_kmers)), mtype='debug')
            for kmer in res_kmers:  # iterate over all ligands and extract molecules + information
                if len(kmer) > config.MAX_COMPOSITE_LENGTH:
                    write_message("Ligand kmer(s) filtered out with a length of {} fragments ({} allowed).\n".format(len(kmer), config.MAX_COMPOSITE_LENGTH), mtype='debug')
                else:
                    ligands.append(self.extract_ligand(kmer))

        else:
            # Extract peptides from given chains
            self.water = [o for o in pybel.ob.OBResidueIter(self.proteincomplex.OBMol) if o.GetResidueProperty(9)]
            if config.PEPTIDES != []:
                peptide_ligands = [self.getpeptides(chain) for chain in config.PEPTIDES]
            elif config.INTRA is not None:
                peptide_ligands = [self.getpeptides(config.INTRA), ]

            ligands = [p for p in peptide_ligands if p is not None ]
            self.covalent, self.lignames_kept, self.lignames_all = [], [], set()

        return [lig for lig in ligands if len(lig.mol.atoms) != 0]

    def extract_ligand(self, kmer):
        """Extract the ligand by copying atoms and bonds and assign all information necessary for later steps."""
        data = namedtuple('ligand', 'mol hetid chain position water members longname type atomorder can_to_pdb')
        members = [(res.GetName(), res.GetChain(), int32_to_negative(res.GetNum())) for res in kmer]
        members = sort_members_by_importance(members)
        rname, rchain, rnum = members[0]
        write_message("Finalizing extraction for ligand %s:%s:%s with %i elements\n" % (rname, rchain, rnum, len(kmer)), mtype='debug')
        names = [x[0] for x in members]
        longname = '-'.join([x[0] for x in members])

        if config.PEPTIDES != []:
            ligtype = 'PEPTIDE'
        elif config.INTRA is not None:
            ligtype = 'INTRA'
        else:
            # Classify a ligand by its HETID(s)
            ligtype = classify_by_name(names)
        write_message("Ligand classified as {}\n".format(ligtype), mtype='debug')

        hetatoms = set()
        for obresidue in kmer:
            hetatoms_res = set([(obatom.GetIdx(), obatom) for obatom in pybel.ob.OBResidueAtomIter(obresidue)
                                if obatom.GetAtomicNum() != 1])

            if not config.ALTLOC:
                # Remove alternative conformations (standard -> True)
                hetatoms_res = set([atm for atm in hetatoms_res
                                    if not self.mapper.mapid(atm[0], mtype='protein',
                                                             to='internal') in self.altconformations])
            hetatoms.update(hetatoms_res)
        write_message("Hetero atoms determined (n={})\n".format(len(hetatoms)), mtype='debug')

        hetatoms = dict(hetatoms)  # make it a dict with idx as key and OBAtom as value
        lig = pybel.ob.OBMol()  # new ligand mol
        neighbours = dict()
        for obatom in hetatoms.values():  # iterate over atom objects
            idx = obatom.GetIdx()
            lig.AddAtom(obatom)
            # ids of all neighbours of obatom
            neighbours[idx] = set([neighbour_atom.GetIdx() for neighbour_atom
                                   in pybel.ob.OBAtomAtomIter(obatom)]) & set(hetatoms.keys())
        write_message("Atom neighbours mapped\n", mtype='debug')

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
        lig.data.update({'Name': rname, 'Chain': rchain, 'ResNr': rnum})

        # Check if a negative residue number is represented as a 32 bit integer
        if rnum > 10 ** 5:
            rnum = int32_to_negative(rnum)

        lig.title = ':'.join((rname, rchain, str(rnum)))
        self.mapper.ligandmaps[lig.title] = mapold

        write_message("Renumerated molecule generated\n", mtype='debug')

        if not config.NOPDBCANMAP:
            atomorder = canonicalize(lig)
        else:
            atomorder =  None

        can_to_pdb = {}
        if atomorder is not None:
            can_to_pdb = {atomorder[key-1]: mapold[key] for key in mapold}

        ligand = data(mol=lig, hetid=rname, chain=rchain, position=rnum, water=self.water,
                      members=members, longname=longname, type=ligtype, atomorder=atomorder,
                      can_to_pdb=can_to_pdb)
        return ligand

    def is_het_residue(self, obres):
        """Given an OBResidue, determines if the residue is indeed a possible ligand
        in the PDB file"""
        if not obres.GetResidueProperty(0):
            # If the residue is NOT amino (0)
            # It can be amino_nucleo, coenzme, ion, nucleo, protein, purine, pyrimidine, solvent
            # In these cases, it is a ligand candidate
            return True
        else:
            # Here, the residue is classified as amino
            # Amino acids can still be ligands, so we check for HETATM entries
            # Only residues with at least one HETATM entry are processed as ligands
            het_atoms = []
            for atm in pybel.ob.OBResidueAtomIter(obres):
                het_atoms.append(obres.IsHetAtom(atm))
            if True in het_atoms:
                return True
        return False


    def filter_for_ligands(self):
        """Given an OpenBabel Molecule, get all ligands, their names, and water"""

        candidates1 = [o for o in pybel.ob.OBResidueIter(self.proteincomplex.OBMol) if not o.GetResidueProperty(9) and self.is_het_residue(o)]

        if config.DNARECEPTOR: # If DNA is the receptor, don't consider DNA as a ligand
            candidates1 = [res for res in candidates1 if res.GetName() not in config.DNA+config.RNA]
        all_lignames = set([a.GetName() for a in candidates1])

        water = [o for o in pybel.ob.OBResidueIter(self.proteincomplex.OBMol) if o.GetResidueProperty(9)]
        # Filter out non-ligands
        if not config.KEEPMOD:  # Keep modified residues as ligands
            candidates2 = [a for a in candidates1 if is_lig(a.GetName()) and a.GetName() not in self.modresidues]
        else:
            candidates2 = [a for a in candidates1 if is_lig(a.GetName())]
        write_message("%i ligand(s) after first filtering step.\n" % len(candidates2), mtype='debug')

        ############################################
        # Filtering by counting and artifacts list #
        ############################################
        artifacts = []
        unique_ligs = set(a.GetName() for a in candidates2)
        for ulig in unique_ligs:
            # Discard if appearing 15 times or more and is possible artifact
            if ulig in config.biolip_list and [a.GetName() for a in candidates2].count(ulig) >= 15:
                artifacts.append(ulig)

        selected_ligands = [a for a in candidates2 if a.GetName() not in artifacts]

        return selected_ligands, all_lignames, water

    def identify_kmers(self, residues):
        """Using the covalent linkage information, find out which fragments/subunits form a ligand."""

        # Remove all those not considered by ligands and pairings including alternate conformations
        ligdoubles = [[(link.id1, link.chain1, link.pos1),
                       (link.id2, link.chain2, link.pos2)] for link in
                      [c for c in self.covalent if c.id1 in self.lignames_kept and c.id2 in self.lignames_kept and
                       c.conf1 in ['A', ''] and c.conf2 in ['A', '']
                      and (c.id1, c.chain1, c.pos1) in residues
                      and (c.id2, c.chain2, c.pos2) in residues]]
        kmers = cluster_doubles(ligdoubles)
        if not kmers:  # No ligand kmers, just normal independent ligands
            return [[residues[res]] for res in residues]

        else:
            # res_kmers contains clusters of covalently bound ligand residues (kmer ligands)
            res_kmers = [[residues[res] for res in kmer] for kmer in kmers]

            # In this case, add other ligands which are not part of a kmer
            in_kmer = []
            for res_kmer in res_kmers:
                for res in res_kmer:
                    in_kmer.append((res.GetName(), res.GetChain(), res.GetNum()))
            for res in residues:
                if res not in in_kmer:
                    newres = [residues[res], ]
                    res_kmers.append(newres)
            return res_kmers


class Mapper:
    """Provides functions for mapping atom IDs in the correct way"""
    def __init__(self):
        self.proteinmap = None  # Map internal atom IDs of protein residues to original PDB Atom IDs
        self.ligandmaps = {}  # Map IDs of new ligand molecules to internal IDs (or PDB IDs?)
        self.original_structure = None

    def mapid(self, idx, mtype, bsid=None, to='original'):  # Mapping to original IDs is standard for ligands
        if mtype == 'reversed': # Needed to map internal ID back to original protein ID
            return self.reversed_proteinmap[idx]
        if mtype == 'protein':
            return self.proteinmap[idx]
        elif mtype == 'ligand':
            if to == 'internal':
                return self.ligandmaps[bsid][idx]
            elif to == 'original':
                return self.proteinmap[self.ligandmaps[bsid][idx]]

    def id_to_atom(self, idx):
        """Returns the atom for a given original ligand ID.
        To do this, the ID is mapped to the protein first and then the atom returned.
        """
        mapped_idx = self.mapid(idx, 'reversed')
        return pybel.Atom(self.original_structure.GetAtom(mapped_idx))



class Mol:
    def __init__(self, altconf, mapper, mtype, bsid):
        self.mtype = mtype
        self.bsid = bsid
        self.rings = None
        self.hydroph_atoms = None
        self.charged = None
        self.hbond_don_atom_pairs = None
        self.hbond_acc_atoms = None
        self.altconf = altconf
        self.Mapper = mapper

    def hydrophobic_atoms(self, all_atoms):
        """Select all carbon atoms which have only carbons and/or hydrogens as direct neighbors."""
        atom_set = []
        data = namedtuple('hydrophobic', 'atom orig_atom orig_idx')
        atm = [a for a in all_atoms if a.atomicnum == 6 and set([natom.GetAtomicNum() for natom
                                                                 in pybel.ob.OBAtomAtomIter(a.OBAtom)]).issubset(
            {1, 6})]
        for atom in atm:
            orig_idx = self.Mapper.mapid(atom.idx, mtype=self.mtype, bsid=self.bsid)
            orig_atom = self.Mapper.id_to_atom(orig_idx)
            if atom.idx not in self.altconf:
                atom_set.append(data(atom=atom, orig_atom=orig_atom, orig_idx=orig_idx))
        return atom_set

    def find_hba(self, all_atoms):
        """Find all possible hydrogen bond acceptors"""
        data = namedtuple('hbondacceptor', 'a a_orig_atom a_orig_idx type')
        a_set = []
        for atom in filter(lambda at: at.OBAtom.IsHbondAcceptor(), all_atoms):
            if atom.atomicnum not in [9, 17, 35, 53] and atom.idx not in self.altconf:  # Exclude halogen atoms
                a_orig_idx = self.Mapper.mapid(atom.idx, mtype=self.mtype, bsid=self.bsid)
                a_orig_atom = self.Mapper.id_to_atom(a_orig_idx)
                a_set.append(data(a=atom, a_orig_atom=a_orig_atom, a_orig_idx=a_orig_idx, type='regular'))
        return a_set

    def find_hbd(self, all_atoms, hydroph_atoms):
        """Find all possible strong and weak hydrogen bonds donors (all hydrophobic C-H pairings)"""
        donor_pairs = []
        data = namedtuple('hbonddonor', 'd d_orig_atom d_orig_idx h type')
        for donor in [a for a in all_atoms if a.OBAtom.IsHbondDonor() and a.idx not in self.altconf]:
            in_ring = False
            if not in_ring:
                for adj_atom in [a for a in pybel.ob.OBAtomAtomIter(donor.OBAtom) if a.IsHbondDonorH()]:
                    d_orig_idx = self.Mapper.mapid(donor.idx, mtype=self.mtype, bsid=self.bsid)
                    d_orig_atom = self.Mapper.id_to_atom(d_orig_idx)
                    donor_pairs.append(data(d=donor, d_orig_atom=d_orig_atom, d_orig_idx=d_orig_idx, h=pybel.Atom(adj_atom), type='regular'))
        for carbon in hydroph_atoms:
            for adj_atom in [a for a in pybel.ob.OBAtomAtomIter(carbon.atom.OBAtom) if a.GetAtomicNum() == 1]:
                d_orig_idx = self.Mapper.mapid(carbon.atom.idx, mtype=self.mtype, bsid=self.bsid)
                d_orig_atom = self.Mapper.id_to_atom(d_orig_idx)
                donor_pairs.append(data(d=carbon, d_orig_atom=d_orig_atom, d_orig_idx=d_orig_idx, h=pybel.Atom(adj_atom), type='weak'))
        return donor_pairs


    def find_rings(self, mol, all_atoms):
        """Find rings and return only aromatic.
        Rings have to be sufficiently planar OR be detected by OpenBabel as aromatic."""
        data = namedtuple('aromatic_ring', 'atoms orig_atoms atoms_orig_idx normal obj center type')
        rings = []
        aromatic_amino = ['TYR', 'TRP', 'HIS', 'PHE']
        ring_candidates = mol.OBMol.GetSSSR()
        write_message("Number of aromatic ring candidates: %i\n" % len(ring_candidates), mtype="debug")
        # Check here first for ligand rings not being detected as aromatic by Babel and check for planarity
        for ring in ring_candidates:
            r_atoms = [a for a in all_atoms if ring.IsMember(a.OBAtom)]
            if 4 < len(r_atoms) <= 6:
                res = list(set([whichrestype(a) for a in r_atoms]))
                if ring.IsAromatic() or res[0] in aromatic_amino or ring_is_planar(ring, r_atoms):
                    # Causes segfault with OpenBabel 2.3.2, so deactivated
                    #typ = ring.GetType() if not ring.GetType() == '' else 'unknown'
                    # Alternative typing
                    typ = '%s-membered' % len(r_atoms)
                    ring_atms = [r_atoms[a].coords for a in [0, 2, 4]]  # Probe atoms for normals, assuming planarity
                    ringv1 = vector(ring_atms[0], ring_atms[1])
                    ringv2 = vector(ring_atms[2], ring_atms[0])
                    atoms_orig_idx = [self.Mapper.mapid(r_atom.idx, mtype=self.mtype, bsid=self.bsid) for r_atom in r_atoms]
                    orig_atoms = [self.Mapper.id_to_atom(idx) for idx in atoms_orig_idx]
                    rings.append(data(atoms=r_atoms,
                                  orig_atoms=orig_atoms,
                                  atoms_orig_idx=atoms_orig_idx,
                                  normal=normalize_vector(np.cross(ringv1, ringv2)),
                                  obj=ring,
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
        self.lig_members = lig_obj.members
        self.pdbid = protcomplex.pymol_name
        self.bindingsite = bs_obj
        self.Mapper = protcomplex.Mapper
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

        self.all_itypes = self.saltbridge_lneg + self.saltbridge_pneg + self.hbonds_pdon
        self.all_itypes = self.all_itypes + self.hbonds_ldon + self.pistacking + self.pication_laro + self.pication_paro
        self.all_itypes = self.all_itypes + self.hydrophobic_contacts + self.halogen_bonds + self.water_bridges
        self.all_itypes = self.all_itypes + self.metal_complexes

        self.no_interactions = all(len(i) == 0 for i in self.all_itypes)
        self.unpaired_hba, self.unpaired_hbd, self.unpaired_hal = self.find_unpaired_ligand()
        self.unpaired_hba_orig_idx = [self.Mapper.mapid(atom.idx, mtype='ligand', bsid=self.ligand.bsid)
                                      for atom in self.unpaired_hba]
        self.unpaired_hbd_orig_idx = [self.Mapper.mapid(atom.idx, mtype='ligand', bsid=self.ligand.bsid)
                                      for atom in self.unpaired_hbd]
        self.unpaired_hal_orig_idx = [self.Mapper.mapid(atom.idx, mtype='ligand', bsid=self.ligand.bsid)
                                      for atom in self.unpaired_hal]
        self.num_unpaired_hba, self.num_unpaired_hbd = len(self.unpaired_hba), len(self.unpaired_hbd)
        self.num_unpaired_hal = len(self.unpaired_hal)

        # Exclude empty chains (coming from ligand as a target, from metal complexes)
        self.interacting_chains = sorted(list(set([i.reschain for i in self.all_itypes
                                                   if i.reschain not in [' ', None]])))

        # Get all interacting residues, excluding ligand and water molecules
        self.interacting_res = list(set([''.join([str(i.resnr), str(i.reschain)]) for i in self.all_itypes
                                         if i.restype not in ['LIG', 'HOH']]))
        if len(self.interacting_res) != 0:
            write_message('Ligand interacts with %i binding site residue(s) in chain(s) %s.\n'
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
                write_message('Complex uses %s.\n' % ', '.join(interactions_list), indent=True)
        else:
            write_message('No interactions for this ligand.\n', indent=True)

    def find_unpaired_ligand(self):
        """Identify unpaired functional in groups in ligands, involving H-Bond donors, acceptors, halogen bond donors.
        """
        unpaired_hba, unpaired_hbd, unpaired_hal = [], [], []
        # Unpaired hydrogen bond acceptors/donors in ligand (not used for hydrogen bonds/water, salt bridges/mcomplex)
        involved_atoms = [hbond.a.idx for hbond in self.hbonds_pdon] + [hbond.d.idx for hbond in self.hbonds_ldon]
        [[involved_atoms.append(atom.idx) for atom in sb.negative.atoms] for sb in self.saltbridge_lneg]
        [[involved_atoms.append(atom.idx) for atom in sb.positive.atoms] for sb in self.saltbridge_pneg]
        [involved_atoms.append(wb.a.idx) for wb in self.water_bridges if wb.protisdon]
        [involved_atoms.append(wb.d.idx) for wb in self.water_bridges if not wb.protisdon]
        [involved_atoms.append(mcomplex.target.atom.idx) for mcomplex in self.metal_complexes
         if mcomplex.location == 'ligand']
        for atom in [hba.a for hba in self.ligand.get_hba()]:
            if atom.idx not in involved_atoms:
                unpaired_hba.append(atom)
        for atom in [hbd.d for hbd in self.ligand.get_hbd()]:
            if atom.idx not in involved_atoms:
                unpaired_hbd.append(atom)

        # unpaired halogen bond donors in ligand (not used for the previous + halogen bonds)
        [involved_atoms.append(atom.don.x.idx) for atom in self.halogen_bonds]
        for atom in [haldon.x for haldon in self.ligand.halogenbond_don]:
            if atom.idx not in involved_atoms:
                unpaired_hal.append(atom)

        return unpaired_hba, unpaired_hbd, unpaired_hal

    def refine_hydrophobic(self, all_h, pistacks):
        """Apply several rules to reduce the number of hydrophobic interactions."""
        sel = {}
        #  1. Rings interacting via stacking can't have additional hydrophobic contacts between each other.
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
            write_message('Reduced number of hydrophobic contacts from %i to %i.\n' % (before, reduced), indent=True)
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
    def __init__(self, atoms, protcomplex, cclass, altconf, min_dist, mapper):
        """Find all relevant parts which could take part in interactions"""
        Mol.__init__(self, altconf, mapper, mtype='protein', bsid=None)
        self.complex = cclass
        self.full_mol = protcomplex
        self.all_atoms = atoms
        self.min_dist = min_dist  # Minimum distance of bs res to ligand
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
        data = namedtuple('hal_acceptor', 'o o_orig_idx y y_orig_idx')
        a_set = []
        # All oxygens, nitrogen, sulfurs with neighboring carbon, phosphor, nitrogen or sulfur
        for a in [at for at in atoms if at.atomicnum in [8, 7, 16]]:
            n_atoms = [na for na in pybel.ob.OBAtomAtomIter(a.OBAtom) if na.GetAtomicNum() in [6, 7, 15, 16]]
            if len(n_atoms) == 1:  # Proximal atom
                o_orig_idx = self.Mapper.mapid(a.idx, mtype=self.mtype, bsid=self.bsid)
                y_orig_idx = self.Mapper.mapid(n_atoms[0].GetIdx(), mtype=self.mtype, bsid=self.bsid)
                a_set.append(data(o=a, o_orig_idx=o_orig_idx, y=pybel.Atom(n_atoms[0]), y_orig_idx=y_orig_idx))
        return a_set

    def find_charged(self, mol):
        """Looks for positive charges in arginine, histidine or lysine, for negative in aspartic and glutamic acid."""
        data = namedtuple('pcharge', 'atoms atoms_orig_idx type center restype resnr reschain')
        a_set = []
        # Iterate through all residue, exclude those in chains defined as peptides
        for res in [r for r in pybel.ob.OBResidueIter(mol.OBMol) if not r.GetChain() in config.PEPTIDES]:
            if config.INTRA is not None:
                if res.GetChain() != config.INTRA: continue
            a_contributing = []
            a_contributing_orig_idx = []
            if res.GetName() in ('ARG', 'HIS', 'LYS'):  # Arginine, Histidine or Lysine have charged sidechains
                for a in pybel.ob.OBResidueAtomIter(res):
                    if a.GetType().startswith('N') and res.GetAtomProperty(a, 8) \
                            and not self.Mapper.mapid(a.GetIdx(), mtype='protein') in self.altconf:
                        a_contributing.append(pybel.Atom(a))
                        a_contributing_orig_idx.append(self.Mapper.mapid(a.GetIdx(), mtype='protein'))
                if not len(a_contributing) == 0:
                    a_set.append(data(atoms=a_contributing,
                                      atoms_orig_idx=a_contributing_orig_idx,
                                      type='positive',
                                      center=centroid([ac.coords for ac in a_contributing]),
                                      restype=res.GetName(),
                                      resnr=res.GetNum(),
                                      reschain=res.GetChain()))
            if res.GetName() in ('GLU', 'ASP'):  # Aspartic or Glutamic Acid
                for a in pybel.ob.OBResidueAtomIter(res):
                    if a.GetType().startswith('O') and res.GetAtomProperty(a, 8) \
                            and not self.Mapper.mapid(a.GetIdx(), mtype='protein') in self.altconf:
                        a_contributing.append(pybel.Atom(a))
                        a_contributing_orig_idx.append(self.Mapper.mapid(a.GetIdx(), mtype='protein'))
                if not len(a_contributing) == 0:
                    a_set.append(data(atoms=a_contributing,
                                      atoms_orig_idx=a_contributing_orig_idx,
                                      type='negative',
                                      center=centroid([ac.coords for ac in a_contributing]),
                                      restype=res.GetName(),
                                      resnr=res.GetNum(),
                                      reschain=res.GetChain()))
        return a_set

    def find_metal_binding(self, mol):
        """Looks for atoms that could possibly be involved in chelating a metal ion.
        This can be any main chain oxygen atom or oxygen, nitrogen and sulfur from specific amino acids"""
        data = namedtuple('metal_binding', 'atom atom_orig_idx type restype resnr reschain location')
        a_set = []
        for res in pybel.ob.OBResidueIter(mol.OBMol):
            restype, reschain, resnr = res.GetName().upper(), res.GetChain(), res.GetNum()
            if restype in ['ASP', 'GLU', 'SER', 'THR', 'TYR']:  # Look for oxygens here
                for a in pybel.ob.OBResidueAtomIter(res):
                    if a.GetType().startswith('O') and res.GetAtomProperty(a, 8) \
                            and not self.Mapper.mapid(a.GetIdx(), mtype='protein') in self.altconf:
                        atom_orig_idx = self.Mapper.mapid(a.GetIdx(), mtype=self.mtype, bsid=self.bsid)
                        a_set.append(data(atom=pybel.Atom(a), atom_orig_idx=atom_orig_idx, type='O', restype=restype,
                                          resnr=resnr, reschain=reschain,
                                          location='protein.sidechain'))
            if restype == 'HIS':  # Look for nitrogen here
                for a in pybel.ob.OBResidueAtomIter(res):
                    if a.GetType().startswith('N') and res.GetAtomProperty(a, 8) \
                            and not self.Mapper.mapid(a.GetIdx(), mtype='protein') in self.altconf:
                        atom_orig_idx = self.Mapper.mapid(a.GetIdx(), mtype=self.mtype, bsid=self.bsid)
                        a_set.append(data(atom=pybel.Atom(a), atom_orig_idx=atom_orig_idx, type='N', restype=restype,
                                          resnr=resnr, reschain=reschain,
                                          location='protein.sidechain'))
            if restype == 'CYS':  # Look for sulfur here
                for a in pybel.ob.OBResidueAtomIter(res):
                    if a.GetType().startswith('S') and res.GetAtomProperty(a, 8) \
                            and not self.Mapper.mapid(a.GetIdx(), mtype='protein') in self.altconf:
                        atom_orig_idx = self.Mapper.mapid(a.GetIdx(), mtype=self.mtype, bsid=self.bsid)
                        a_set.append(data(atom=pybel.Atom(a), atom_orig_idx=atom_orig_idx, type='S', restype=restype,
                                          resnr=resnr, reschain=reschain,
                                          location='protein.sidechain'))
            for a in pybel.ob.OBResidueAtomIter(res):  # All main chain oxygens
                if a.GetType().startswith('O') and res.GetAtomProperty(a, 2) \
                        and not self.Mapper.mapid(a.GetIdx(), mtype='protein') in self.altconf and restype != 'HOH':
                    atom_orig_idx = self.Mapper.mapid(a.GetIdx(), mtype=self.mtype, bsid=self.bsid)
                    a_set.append(data(atom=pybel.Atom(a), atom_orig_idx=atom_orig_idx, type='O', restype=res.GetName(),
                                      resnr=res.GetNum(), reschain=res.GetChain(),
                                      location='protein.mainchain'))
        return a_set


class Ligand(Mol):
    def __init__(self, cclass, ligand):
        altconf = cclass.altconf
        self.hetid, self.chain, self.position = ligand.hetid, ligand.chain, ligand.position
        self.bsid = ':'.join([self.hetid, self.chain, str(self.position)])
        Mol.__init__(self, altconf, cclass.Mapper, mtype='ligand', bsid=self.bsid)
        self.members = ligand.members
        self.longname = ligand.longname
        self.type = ligand.type
        self.complex = cclass
        self.molecule = ligand.mol  # Pybel Molecule
        self.smiles = self.molecule.write(format='can')  # SMILES String
        self.inchikey = self.molecule.write(format='inchikey')
        self.can_to_pdb = ligand.can_to_pdb
        if not len(self.smiles) == 0:
            self.smiles = self.smiles.split()[0]
        else:
            write_message('Could not write SMILES for this ligand.\n', indent=True, mtype='warning')
            self.smiles = ''
        self.heavy_atoms = self.molecule.OBMol.NumHvyAtoms()  # Heavy atoms count
        self.all_atoms = self.molecule.atoms
        self.atmdict = {l.idx: l for l in self.all_atoms}
        self.rings = self.find_rings(self.molecule, self.all_atoms)
        self.hydroph_atoms = self.hydrophobic_atoms(self.all_atoms)
        self.hbond_acc_atoms = self.find_hba(self.all_atoms)
        self.num_rings = len(self.rings)
        if self.num_rings != 0:
            write_message('Contains %i aromatic ring(s).\n' % self.num_rings, indent=True)
        descvalues = self.molecule.calcdesc()
        self.molweight, self.logp = float(descvalues['MW']), float(descvalues['logP'])
        self.num_rot_bonds = int(self.molecule.OBMol.NumRotors())
        self.atomorder = ligand.atomorder

        ##########################################################
        # Special Case for hydrogen bond acceptor identification #
        ##########################################################

        self.inverse_mapping = {v: k for k, v in self.Mapper.ligandmaps[self.bsid].items()}
        self.pdb_to_idx_mapping = {v: k for k, v in self.Mapper.proteinmap.items()}
        self.hbond_don_atom_pairs = self.find_hbd(self.all_atoms, self.hydroph_atoms)

        ######
        donor_pairs = []
        data = namedtuple('hbonddonor', 'd d_orig_atom d_orig_idx h type')
        for donor in self.all_atoms:
            pdbidx = self.Mapper.mapid(donor.idx, mtype='ligand', bsid=self.bsid, to='original')
            d = cclass.atoms[self.pdb_to_idx_mapping[pdbidx]]
            if d.OBAtom.IsHbondDonor():
                for adj_atom in [a for a in pybel.ob.OBAtomAtomIter(d.OBAtom) if a.IsHbondDonorH()]:
                    d_orig_atom = self.Mapper.id_to_atom(pdbidx)
                    donor_pairs.append(data(d=donor, d_orig_atom=d_orig_atom, d_orig_idx=pdbidx, h=pybel.Atom(adj_atom), type='regular'))
        self.hbond_don_atom_pairs = donor_pairs
        #######

        self.charged = self.find_charged(self.all_atoms)
        self.centroid = centroid([a.coords for a in self.all_atoms])
        self.max_dist_to_center = max((euclidean3d(self.centroid, a.coords) for a in self.all_atoms))
        self.water = []
        data = namedtuple('water', 'oxy oxy_orig_idx')
        for hoh in ligand.water:
            oxy = None
            for at in pybel.ob.OBResidueAtomIter(hoh):
                if at.GetAtomicNum() == 8 and at.GetIdx() not in self.altconf:
                    oxy = pybel.Atom(at)
            # There are some cases where there is no oxygen in a water residue, ignore those
            if not set([at.GetAtomicNum() for at in pybel.ob.OBResidueAtomIter(hoh)]) == {1} and oxy is not None:
                if euclidean3d(self.centroid, oxy.coords) < self.max_dist_to_center + config.BS_DIST:
                    oxy_orig_idx = self.Mapper.mapid(oxy.idx, mtype='protein')
                    self.water.append(data(oxy=oxy, oxy_orig_idx=oxy_orig_idx))
        self.halogenbond_don = self.find_hal(self.all_atoms)
        self.metal_binding = self.find_metal_binding(self.all_atoms, self.water)
        self.metals = []
        data = namedtuple('metal', 'm orig_m m_orig_idx')
        for a in [a for a in self.all_atoms if a.type.upper() in config.METAL_IONS]:
            m_orig_idx = self.Mapper.mapid(a.idx, mtype=self.mtype, bsid=self.bsid)
            orig_m = self.Mapper.id_to_atom(m_orig_idx)
            self.metals.append(data(m=a, m_orig_idx=m_orig_idx, orig_m=orig_m))
        self.num_hba, self.num_hbd = len(self.hbond_acc_atoms), len(self.hbond_don_atom_pairs)
        self.num_hal = len(self.halogenbond_don)

    def get_canonical_num(self, atomnum):
        """Converts internal atom ID into canonical atom ID. Agrees with Canonical SMILES in XML."""
        return self.atomorder[atomnum-1]


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
        data = namedtuple('hal_donor', 'x orig_x x_orig_idx c c_orig_idx')
        a_set = []
        for a in atoms:
            if self.is_functional_group(a, 'halocarbon'):
                n_atoms = [na for na in pybel.ob.OBAtomAtomIter(a.OBAtom) if na.GetAtomicNum() == 6]
                x_orig_idx = self.Mapper.mapid(a.idx, mtype=self.mtype, bsid=self.bsid)
                orig_x = self.Mapper.id_to_atom(x_orig_idx)
                c_orig_idx = [self.Mapper.mapid(na.GetIdx(), mtype=self.mtype, bsid=self.bsid) for na in n_atoms]
                a_set.append(data(x=a, orig_x=orig_x, x_orig_idx=x_orig_idx, c=pybel.Atom(n_atoms[0]), c_orig_idx=c_orig_idx))
        if len(a_set) != 0:
            write_message('Ligand contains %i halogen atom(s).\n' % len(a_set), indent=True)
        return a_set

    def find_charged(self, all_atoms):
        """Identify all positively charged groups in a ligand. This search is not exhaustive, as the cases can be quite
        diverse. The typical cases seem to be protonated amines, quaternary ammoinium and sulfonium
        as mentioned in 'Cation-pi interactions in ligand recognition and catalysis' (Zacharias et al., 2002)).
        Identify negatively charged groups in the ligand.
        """
        data = namedtuple('lcharge', 'atoms orig_atoms atoms_orig_idx type center fgroup')
        a_set = []
        for a in all_atoms:
            a_orig_idx = self.Mapper.mapid(a.idx, mtype=self.mtype, bsid=self.bsid)
            a_orig = self.Mapper.id_to_atom(a_orig_idx)
            if self.is_functional_group(a, 'quartamine'):
                a_set.append(data(atoms=[a, ], orig_atoms = [a_orig, ], atoms_orig_idx=[a_orig_idx, ], type='positive',
                                  center=list(a.coords), fgroup='quartamine'))
            elif self.is_functional_group(a, 'tertamine'):
                a_set.append(data(atoms=[a, ], orig_atoms = [a_orig, ], atoms_orig_idx=[a_orig_idx, ], type='positive', center=list(a.coords),
                                  fgroup='tertamine'))
            if self.is_functional_group(a, 'sulfonium'):
                a_set.append(data(atoms=[a, ], orig_atoms = [a_orig, ], atoms_orig_idx=[a_orig_idx, ], type='positive', center=list(a.coords),
                                  fgroup='sulfonium'))
            if self.is_functional_group(a, 'phosphate'):
                a_contributing = [a, ]
                a_contributing_orig_idx = [a_orig_idx, ]
                [a_contributing.append(pybel.Atom(neighbor)) for neighbor in pybel.ob.OBAtomAtomIter(a.OBAtom)]
                [a_contributing_orig_idx.append(self.Mapper.mapid(neighbor.idx, mtype=self.mtype, bsid=self.bsid))
                 for neighbor in a_contributing]
                orig_contributing = [self.Mapper.id_to_atom(idx) for idx in a_contributing_orig_idx]
                a_set.append(data(atoms=a_contributing, orig_atoms=orig_contributing, atoms_orig_idx=a_contributing_orig_idx, type='negative',
                                  center=a.coords, fgroup='phosphate'))
            if self.is_functional_group(a, 'sulfonicacid'):
                a_contributing = [a, ]
                a_contributing_orig_idx = [a_orig_idx, ]
                [a_contributing.append(pybel.Atom(neighbor)) for neighbor in pybel.ob.OBAtomAtomIter(a.OBAtom) if
                 neighbor.GetAtomicNum() == 8]
                [a_contributing_orig_idx.append(self.Mapper.mapid(neighbor.idx, mtype=self.mtype, bsid=self.bsid))
                 for neighbor in a_contributing]
                orig_contributing = [self.Mapper.id_to_atom(idx) for idx in a_contributing_orig_idx]
                a_set.append(data(atoms=a_contributing, orig_atoms=orig_contributing, atoms_orig_idx=a_contributing_orig_idx, type='negative',
                                  center=a.coords, fgroup='sulfonicacid'))
            elif self.is_functional_group(a, 'sulfate'):
                a_contributing = [a, ]
                a_contributing_orig_idx = [a_orig_idx, ]
                [a_contributing_orig_idx.append(self.Mapper.mapid(neighbor.idx, mtype=self.mtype, bsid=self.bsid))
                 for neighbor in a_contributing]
                [a_contributing.append(pybel.Atom(neighbor)) for neighbor in pybel.ob.OBAtomAtomIter(a.OBAtom)]
                orig_contributing = [self.Mapper.id_to_atom(idx) for idx in a_contributing_orig_idx]
                a_set.append(data(atoms=a_contributing, orig_atoms=orig_contributing, atoms_orig_idx=a_contributing_orig_idx, type='negative',
                                  center=a.coords, fgroup='sulfate'))
            if self.is_functional_group(a, 'carboxylate'):
                a_contributing = [pybel.Atom(neighbor) for neighbor in pybel.ob.OBAtomAtomIter(a.OBAtom)
                                  if neighbor.GetAtomicNum() == 8]
                a_contributing_orig_idx = [self.Mapper.mapid(neighbor.idx, mtype=self.mtype, bsid=self.bsid)
                                           for neighbor in a_contributing]
                orig_contributing = [self.Mapper.id_to_atom(idx) for idx in a_contributing_orig_idx]
                a_set.append(data(atoms=a_contributing, orig_atoms=orig_contributing, atoms_orig_idx=a_contributing_orig_idx, type='negative',
                                  center=centroid([a.coords for a in a_contributing]), fgroup='carboxylate'))
            elif self.is_functional_group(a, 'guanidine'):
                a_contributing = [pybel.Atom(neighbor) for neighbor in pybel.ob.OBAtomAtomIter(a.OBAtom)
                                  if neighbor.GetAtomicNum() == 7]
                a_contributing_orig_idx = [self.Mapper.mapid(neighbor.idx, mtype=self.mtype, bsid=self.bsid)
                                           for neighbor in a_contributing]
                orig_contributing = [self.Mapper.id_to_atom(idx) for idx in a_contributing_orig_idx]
                a_set.append(data(atoms=a_contributing, orig_atoms=orig_contributing, atoms_orig_idx=a_contributing_orig_idx, type='positive',
                                  center=a.coords, fgroup='guanidine'))
        return a_set

    def find_metal_binding(self, lig_atoms, water_oxygens):
        """Looks for atoms that could possibly be involved in binding a metal ion.
        This can be any water oxygen, as well as oxygen from carboxylate, phophoryl, phenolate, alcohol;
        nitrogen from imidazole; sulfur from thiolate.
        """
        a_set = []
        data = namedtuple('metal_binding', 'atom orig_atom atom_orig_idx type fgroup restype resnr reschain location')
        for oxygen in water_oxygens:
            a_set.append(data(atom=oxygen.oxy, atom_orig_idx=oxygen.oxy_orig_idx, type='O', fgroup='water',
                              restype=whichrestype(oxygen.oxy), resnr=whichresnumber(oxygen.oxy),
                              reschain=whichchain(oxygen.oxy), location='water', orig_atom=self.Mapper.id_to_atom(oxygen.oxy_orig_idx)))
        # #@todo Refactor code
        for a in lig_atoms:
            a_orig_idx = self.Mapper.mapid(a.idx, mtype='ligand', bsid=self.bsid)
            n_atoms = pybel.ob.OBAtomAtomIter(a.OBAtom)  # Neighboring atoms
            # All atomic numbers of neighboring atoms
            n_atoms_atomicnum = [n.GetAtomicNum() for n in pybel.ob.OBAtomAtomIter(a.OBAtom)]
            if a.atomicnum == 8:  # Oxygen
                if n_atoms_atomicnum.count('1') == 1 and len(n_atoms_atomicnum) == 2:  # Oxygen in alcohol (R-[O]-H)
                    a_set.append(data(atom=a, atom_orig_idx=a_orig_idx, type='O', fgroup='alcohol',
                                      restype=self.hetid, resnr=self.position, reschain=self.chain,
                                      location='ligand', orig_atom=self.Mapper.id_to_atom(a_orig_idx)))
                if True in [n.IsAromatic() for n in n_atoms] and not a.OBAtom.IsAromatic():  # Phenolate oxygen
                    a_set.append(data(atom=a, atom_orig_idx=a_orig_idx, type='O', fgroup='phenolate',
                                      restype=self.hetid, resnr=self.position, reschain=self.chain,
                                      location='ligand', orig_atom=self.Mapper.id_to_atom(a_orig_idx)))
            if a.atomicnum == 6:  # It's a carbon atom
                if n_atoms_atomicnum.count(8) == 2 and n_atoms_atomicnum.count(6) == 1:  # It's a carboxylate group
                    for neighbor in [n for n in n_atoms if n.GetAtomicNum() == 8]:
                        neighbor_orig_idx = self.Mapper.mapid(neighbor.GetIdx(), mtype='ligand', bsid=self.bsid)
                        a_set.append(data(atom=pybel.Atom(neighbor), atom_orig_idx=neighbor_orig_idx, type='O',
                                          fgroup='carboxylate',
                                          restype=self.hetid,
                                          resnr=self.position, reschain=self.chain,
                                          location='ligand', orig_atom=self.Mapper.id_to_atom(a_orig_idx)))
            if a.atomicnum == 15:  # It's a phosphor atom
                if n_atoms_atomicnum.count(8) >= 3:  # It's a phosphoryl
                    for neighbor in [n for n in n_atoms if n.GetAtomicNum() == 8]:
                        neighbor_orig_idx = self.Mapper.mapid(neighbor.GetIdx(), mtype='ligand', bsid=self.bsid)
                        a_set.append(data(atom=pybel.Atom(neighbor), atom_orig_idx=neighbor_orig_idx, type='O',
                                          fgroup='phosphoryl',
                                          restype=self.hetid,
                                          resnr=self.position, reschain=self.chain,
                                          location='ligand', orig_atom=self.Mapper.id_to_atom(a_orig_idx)))
                if n_atoms_atomicnum.count(8) == 2:  # It's another phosphor-containing group #@todo (correct name?)
                    for neighbor in [n for n in n_atoms if n.GetAtomicNum() == 8]:
                        neighbor_orig_idx = self.Mapper.mapid(neighbor.GetIdx(), mtype='ligand', bsid=self.bsid)
                        a_set.append(data(atom=pybel.Atom(neighbor), atom_orig_idx=neighbor_orig_idx, type='O',
                                          fgroup='phosphor.other', restype=self.hetid,
                                          resnr=self.position,
                                          reschain=self.chain, location='ligand', orig_atom=self.Mapper.id_to_atom(a_orig_idx)))
            if a.atomicnum == 7:  # It's a nitrogen atom
                if n_atoms_atomicnum.count(6) == 2:  # It's imidazole/pyrrole or similar
                    a_set.append(data(atom=a, atom_orig_idx=a_orig_idx, type='N', fgroup='imidazole/pyrrole',
                                      restype=self.hetid, resnr=self.position, reschain=self.chain,
                                      location='ligand', orig_atom=self.Mapper.id_to_atom(a_orig_idx)))
            if a.atomicnum == 16:  # It's a sulfur atom
                if True in [n.IsAromatic() for n in n_atoms] and not a.OBAtom.IsAromatic():  # Thiolate
                    a_set.append(data(atom=a, atom_orig_idx=a_orig_idx, type='S', fgroup='thiolate',
                                      restype=self.hetid, resnr=self.position, reschain=self.chain,
                                      location='ligand', orig_atom=self.Mapper.id_to_atom(a_orig_idx)))
                if set(n_atoms_atomicnum) == {26}:  # Sulfur in Iron sulfur cluster
                    a_set.append(data(atom=a, atom_orig_idx=a_orig_idx, type='S', fgroup='iron-sulfur.cluster',
                                      restype=self.hetid, resnr=self.position, reschain=self.chain,
                                      location='ligand', orig_atom=self.Mapper.id_to_atom(a_orig_idx)))

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
        self.information = {}
        self.corrected_pdb = ''
        self._output_path = '/tmp'
        self.pymol_name = None
        self.modres = set()
        self.resis = []
        self.altconf = []  # Atom idx of atoms with alternate conformations
        self.covalent = []  # Covalent linkages between ligands and protein residues/other ligands
        self.excluded = []  # Excluded ligands
        self.Mapper = Mapper()
        self.ligands = []

    def __str__(self):
        formatted_lig_names = [":".join([x.hetid, x.chain, str(x.position)]) for x in self.ligands]
        return "Protein structure %s with ligands:\n" % (self.pymol_name) + "\n".join([lig for lig in formatted_lig_names])

    def load_pdb(self, pdbpath, as_string=False):
        """Loads a pdb file with protein AND ligand(s), separates and prepares them.
        If specified 'as_string', the input is a PDB string instead of a path."""
        if as_string:
            self.sourcefiles['pdbcomplex.original'] = None
            self.sourcefiles['pdbcomplex'] = None
            self.sourcefiles['pdbstring'] = pdbpath
        else:
            self.sourcefiles['pdbcomplex.original'] = pdbpath
            self.sourcefiles['pdbcomplex'] = pdbpath
        self.information['pdbfixes'] = False
        pdbparser = PDBParser(pdbpath, as_string=as_string)  # Parse PDB file to find errors and get additonal data
        # #@todo Refactor and rename here
        self.Mapper.proteinmap = pdbparser.proteinmap
        self.Mapper.reversed_proteinmap = inv_map = {v: k for k, v in self.Mapper.proteinmap.items()}
        self.modres = pdbparser.modres
        self.covalent = pdbparser.covalent
        self.altconf = pdbparser.altconformations
        self.corrected_pdb = pdbparser.corrected_pdb

        if not config.PLUGIN_MODE:
            if pdbparser.num_fixed_lines > 0:
                write_message('%i lines automatically fixed in PDB input file.\n' % pdbparser.num_fixed_lines)
                # Save modified PDB file
                if not as_string:
                    basename = os.path.basename(pdbpath).split('.')[0]
                else:
                    basename = "from_stdin"
                pdbpath_fixed = tmpfile(prefix='plipfixed.' + basename + '_', direc=self.output_path)
                create_folder_if_not_exists(self.output_path)
                self.sourcefiles['pdbcomplex'] = pdbpath_fixed
                self.corrected_pdb = re.sub(r'[^\x00-\x7F]+', ' ', self.corrected_pdb)  # Strip non-unicode chars
                if not config.NOFIXFILE: # Only write to file if this option is not activated
                    with open(pdbpath_fixed, 'w') as f:
                        f.write(self.corrected_pdb)
                self.information['pdbfixes'] = True


        if not as_string:
            self.sourcefiles['filename'] = os.path.basename(self.sourcefiles['pdbcomplex'])
        self.protcomplex, self.filetype = read_pdb(self.corrected_pdb, as_string=True)

        # Update the model in the Mapper class instance
        self.Mapper.original_structure = self.protcomplex.OBMol
        write_message('PDB structure successfully read.\n')


        # Determine (temporary) PyMOL Name from Filename
        self.pymol_name = pdbpath.split('/')[-1].split('.')[0] + '-Protein'
        # Replace characters causing problems in PyMOL
        self.pymol_name = self.pymol_name.replace(' ', '').replace('(', '').replace(')', '').replace('-','_')
        # But if possible, name it after PDBID in Header
        if 'HEADER' in self.protcomplex.data:  # If the PDB file has a proper header
            potential_name = self.protcomplex.data['HEADER'][56:60].lower()
            if extract_pdbid(potential_name) != 'UnknownProtein':
                self.pymol_name = potential_name
        write_message("Pymol Name set as: '%s'\n" % self.pymol_name, mtype='debug')

        # Extract and prepare ligands
        ligandfinder = LigandFinder(self.protcomplex, self.altconf, self.modres, self.covalent, self.Mapper)
        self.ligands = ligandfinder.ligands
        self.excluded = ligandfinder.excluded

        # Add polar hydrogens
        self.protcomplex.OBMol.AddPolarHydrogens()
        for atm in self.protcomplex:
            self.atoms[atm.idx] = atm
        write_message("Assigned polar hydrogens\n", mtype='debug')

        if len(self.excluded) != 0:
            write_message("Excluded molecules as ligands: %s\n" % ','.join([lig for lig in self.excluded]))

        if config.DNARECEPTOR:
            self.resis = [obres for obres in pybel.ob.OBResidueIter(self.protcomplex.OBMol) if obres.GetName() in config.DNA+config.RNA]
        else:
            self.resis = [obres for obres in pybel.ob.OBResidueIter(self.protcomplex.OBMol) if obres.GetResidueProperty(0)]

        num_ligs = len(self.ligands)
        if num_ligs == 1:
            write_message("Analyzing one ligand...\n")
        elif num_ligs > 1:
            write_message("Analyzing %i ligands...\n" % num_ligs)
        else:
            write_message("Structure contains no ligands.\n\n")

    def analyze(self):
        """Triggers analysis of all complexes in structure"""
        for ligand in self.ligands:
            self.characterize_complex(ligand)

    def characterize_complex(self, ligand):
        """Handles all basic functions for characterizing the interactions for one ligand"""

        single_sites = []
        for member in ligand.members:
            single_sites.append(':'.join([str(x) for x in member]))
        site = ' + '.join(single_sites)
        site = site if not len(site) > 20 else site[:20] + '...'
        longname = ligand.longname if not len(ligand.longname) > 20 else ligand.longname[:20] + '...'
        ligtype = 'Unspecified type' if ligand.type == 'UNSPECIFIED' else ligand.type
        ligtext = "\n%s [%s] -- %s" % (longname, ligtype, site)
        if ligtype == 'PEPTIDE':
            ligtext = '\n Chain %s [PEPTIDE / INTER-CHAIN]' % ligand.chain
        if ligtype == 'INTRA':
            ligtext = "\n Chain %s [INTRA-CHAIN]" % ligand.chain
        any_in_biolip = len(set([x[0] for x in ligand.members]).intersection(config.biolip_list)) != 0
        write_message(ligtext)
        write_message('\n' + '-' * len(ligtext) + '\n')

        if ligtype not in ['POLYMER', 'DNA', 'ION', 'DNA+ION', 'RNA+ION', 'SMALLMOLECULE+ION'] and any_in_biolip:
            write_message('may be biologically irrelevant\n', mtype='info', indent=True)

        lig_obj = Ligand(self, ligand)
        cutoff = lig_obj.max_dist_to_center + config.BS_DIST
        bs_res = self.extract_bs(cutoff, lig_obj.centroid, self.resis)
        # Get a list of all atoms belonging to the binding site, search by idx
        bs_atoms = [self.atoms[idx] for idx in [i for i in self.atoms.keys()
                                                if self.atoms[i].OBAtom.GetResidue().GetIdx() in bs_res]
                    if idx in self.Mapper.proteinmap and self.Mapper.mapid(idx, mtype='protein') not in self.altconf]
        if ligand.type == 'PEPTIDE':
            # If peptide, don't consider the peptide chain as part of the protein binding site
            bs_atoms = [a for a in bs_atoms if a.OBAtom.GetResidue().GetChain() != lig_obj.chain]
        if ligand.type == 'INTRA':
            # Interactions within the chain
            bs_atoms = [a for a in bs_atoms if a.OBAtom.GetResidue().GetChain() == lig_obj.chain]
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
        write_message('Binding site atoms in vicinity (%.1f A max. dist: %i).\n' % (config.BS_DIST, num_bs_atoms),
                indent=True)

        bs_obj = BindingSite(bs_atoms_refined, self.protcomplex, self, self.altconf, min_dist, self.Mapper)
        pli_obj = PLInteraction(lig_obj, bs_obj, self)
        self.interaction_sets[ligand.mol.title] = pli_obj

    def extract_bs(self, cutoff, ligcentroid, resis):
        """Return list of ids from residues belonging to the binding site"""
        return [obres.GetIdx() for obres in resis if self.res_belongs_to_bs(obres, cutoff, ligcentroid)]

    def res_belongs_to_bs(self, res, cutoff, ligcentroid):
        """Check for each residue if its centroid is within a certain distance to the ligand centroid.
        Additionally checks if a residue belongs to a chain restricted by the user (e.g. by defining a peptide chain)"""
        rescentroid = centroid([(atm.x(), atm.y(), atm.z()) for atm in pybel.ob.OBResidueAtomIter(res)])
        # Check geometry
        near_enough = True if euclidean3d(rescentroid, ligcentroid) < cutoff else False
        # Check chain membership
        restricted_chain = True if res.GetChain() in config.PEPTIDES else False
        return (near_enough and not restricted_chain)

    def get_atom(self, idx):
        return self.atoms[idx]

    @property
    def output_path(self):
        return self._output_path

    @output_path.setter
    def output_path(self, path):
        self._output_path = tilde_expansion(path)
