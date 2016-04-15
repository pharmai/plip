"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
visualize.py - Visualization of PLIP results using PyMOL.
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


# Own modules
#from supplemental import *
from collections import namedtuple
import config
from pymol import cmd
from pymolplip import PyMOLVisualizer

# Python Standard Library
import json

###### from supplemental.py

from pymol import finish_launching
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

def message(msg, indent=False, mtype='standard', caption=False):
    """Writes messages in verbose mode"""
    if caption:
        msg = msg + '\n' + '-'*len(msg)
    if mtype == 'warning':
        msg = colorlog('Warning:  ' + msg, 'yellow')
    if mtype == 'error':
        msg = colorlog('Error:  ' + msg, 'red')
    if mtype == 'debug':
        msg = colorlog('Debug:  ' + msg, 'pink')
    if mtype == 'info':
        msg = colorlog('Info:  ' + msg, 'green')
    if indent:
        msg = '  ' + msg
    if mtype in ['error', 'warning']:
        sys.stderr.write(msg)
    else:
        sys.stdout.write(msg)

def write_message(msg, indent=False, mtype='standard', caption=False):
    """Writes message if verbose mode is set."""
    if (mtype=='debug' and config.DEBUG) or (mtype !='debug' and config.VERBOSE) or mtype=='error':
        message(msg, indent=indent, mtype=mtype, caption=caption)

#################



hbonds_info = namedtuple('hbonds_info', 'ldon_id lig_don_id prot_acc_id pdon_id prot_don_id lig_acc_id')
hydrophobic_info = namedtuple('hydrophobic_info', 'bs_ids lig_ids pairs_ids')
halogen_info = namedtuple('halogen_info', 'don_id acc_id')
pistack_info = namedtuple('pistack_info', 'proteinring_atoms, proteinring_center ligandring_atoms '
                                          'ligandring_center type')
pication_info = namedtuple('pication_info', 'ring_center charge_center ring_atoms charge_atoms, protcharged')
sbridge_info = namedtuple('sbridge_info', 'positive_atoms negative_atoms positive_center negative_center protispos')
wbridge_info = namedtuple('wbridge_info', 'don_id acc_id water_id protisdon')
metal_info = namedtuple('metal_info', 'metal_id, target_id location')


def select_by_ids(selname, idlist, selection_exists=False, chunksize=20, restrict=None):
    """Selection with a large number of ids concatenated into a selection
    list can cause buffer overflow in PyMOL. This function takes a selection
    name and and list of IDs (list of integers) as input and makes a careful
    step-by-step selection (packages of 20 by default)"""
    idlist = list(set(idlist))  # Remove duplicates
    if not selection_exists:
        cmd.select(selname, 'None')  # Empty selection first
    idchunks = [idlist[i:i+chunksize] for i in xrange(0, len(idlist), chunksize)]
    for idchunk in idchunks:
        cmd.select(selname, '%s or (id %s)' % (selname, '+'.join(map(str, idchunk))))
    if restrict is not None:
        cmd.select(selname, '%s and %s' % (selname, restrict))

class VisualizerData:
    """Contains all information on a complex relevant for visualization. Can be pickled"""
    def __init__(self, mol, site):
        pcomp = mol
        pli = mol.interaction_sets[site]
        ligand = pli.ligand

        # General Information
        self.lig_members = sorted(pli.ligand.members)
        self.sourcefile = pcomp.sourcefiles['pdbcomplex']
        self.corrected_pdb = pcomp.corrected_pdb
        self.pdbid = mol.pymol_name
        self.hetid = ligand.hetid
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
        self.hydrophobic_contacts = hydrophobic_info(bs_ids=[hp[0] for hp in hydroph_pairs_id],
                                                     lig_ids=[hp[1] for hp in hydroph_pairs_id],
                                                     pairs_ids=hydroph_pairs_id)

        # Hydrogen Bonds
        # #@todo Don't use indices, simplify this code here
        hbonds_ldon, hbonds_pdon = pli.hbonds_ldon, pli.hbonds_pdon
        hbonds_ldon_id = [(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_ldon]
        hbonds_pdon_id = [(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_pdon]
        self.hbonds = hbonds_info(ldon_id=[(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_ldon],
                                  lig_don_id=[hb[1] for hb in hbonds_ldon_id],
                                  prot_acc_id=[hb[0] for hb in hbonds_ldon_id],
                                  pdon_id=[(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_pdon],
                                  prot_don_id=[hb[1] for hb in hbonds_pdon_id],
                                  lig_acc_id=[hb[0] for hb in hbonds_pdon_id])

        # Halogen Bonds
        self.halogen_bonds = [halogen_info(don_id=h.don_orig_idx, acc_id=h.acc_orig_idx)
                              for h in pli.halogen_bonds]

        # Pistacking
        self.pistacking = [pistack_info(proteinring_atoms=pistack.proteinring.atoms_orig_idx,
                                        proteinring_center=pistack.proteinring.center,
                                        ligandring_atoms=pistack.ligandring.atoms_orig_idx,
                                        ligandring_center=pistack.ligandring.center,
                                        type=pistack.type) for pistack in pli.pistacking]

        # Pi-cation interactions
        self.pication = [pication_info(ring_center=picat.ring.center,
                                       charge_center=picat.charge.center,
                                       ring_atoms=picat.ring.atoms_orig_idx,
                                       charge_atoms=picat.charge.atoms_orig_idx,
                                       protcharged=picat.protcharged)
                         for picat in pli.pication_paro+pli.pication_laro]

        # Salt Bridges
        self.saltbridges = [sbridge_info(positive_atoms=sbridge.positive.atoms_orig_idx,
                                         negative_atoms=sbridge.negative.atoms_orig_idx,
                                         positive_center=sbridge.positive.center,
                                         negative_center=sbridge.negative.center,
                                         protispos=sbridge.protispos)
                            for sbridge in pli.saltbridge_lneg+pli.saltbridge_pneg]

        # Water Bridgese('wbridge_info', 'don_id acc_id water_id protisdon')
        self.waterbridges = [wbridge_info(don_id=wbridge.d_orig_idx,
                                          acc_id=wbridge.a_orig_idx,
                                          water_id=wbridge.water_orig_idx,
                                          protisdon=wbridge.protisdon) for wbridge in pli.water_bridges]

        # Metal Complexes
        self.metal_complexes = [metal_info(metal_id=metalc.metal_orig_idx,
                                           target_id=metalc.target_orig_idx,
                                           location=metalc.location) for metalc in pli.metal_complexes]

    def to_json(self):
        """Generates a JSON dump of the class without the contained PDB string
        in self.corrected_pdb"""
        mindict = self.__dict__
        del mindict['corrected_pdb']
        return json.dumps(mindict)



def visualize_in_pymol(plcomplex):
    """Visualizes the protein-ligand pliprofiler at one site in PyMOL."""

    vis = PyMOLVisualizer(plcomplex)



    #####################
    # Set everything up #
    #####################

    pdbid = plcomplex.pdbid
    lig_members = plcomplex.lig_members
    chain = plcomplex.chain
    ligname = plcomplex.hetid
    metal_ids = plcomplex.metal_ids
    metal_ids_str = '+'.join([str(i) for i in metal_ids])

    ########################
    # Basic visualizations #
    ########################

    start_pymol(run=True, options='-pcq', quiet=not config.DEBUG)
    vis.set_initial_representations()

    cmd.load(plcomplex.sourcefile)
    current_name = cmd.get_object_list(selection='(all)')[0]
    write_message('Setting current_name to "%s" and pdbid to "%s\n"' % (current_name, pdbid), mtype='debug')
    cmd.set_name(current_name, pdbid)
    cmd.hide('everything', 'all')
    cmd.select(ligname, 'resn %s and chain %s and resi %s*' % (ligname, chain, plcomplex.position))

    # Visualize and color metal ions if there are any
    if not len(metal_ids) == 0:
        vis.select_by_ids(ligname, metal_ids, selection_exists=True)
        cmd.show('spheres', 'id %s and %s' % (metal_ids_str, pdbid))

    # Additionally, select all members of composite ligands
    for member in lig_members:
        resid, chain, resnr = member[0], member[1], str(member[2])
        cmd.select(ligname, '%s or (resn %s and chain %s and resi %s)' % (ligname, resid, chain, resnr))
    cmd.show('sticks', ligname)
    cmd.color('myblue')
    cmd.color('myorange', ligname)
    cmd.util.cnc('all')
    if not len(metal_ids) == 0:
        cmd.color('hotpink', 'id %s' % metal_ids_str)
        cmd.hide('sticks', 'id %s' % metal_ids_str)
        cmd.set('sphere_scale', 0.3, ligname)
    cmd.deselect()

    vis.make_initial_selections()

    vis.show_hydrophobic()  # Hydrophobic Contacts
    vis.show_hbonds()  # Hydrogen Bonds
    vis.show_halogen()  # Halogen Bonds
    vis.show_stacking()  # pi-Stacking Interactions
    vis.show_cationpi()  # pi-Cation Interactions
    vis.show_sbridges()  # Salt Bridges
    vis.show_wbridges()  # Water Bridges
    vis.show_metal()  # Metal Coordination

    vis.refinements()


    vis.zoom_to_ligand()

    vis.selections_cleanup()
    vis.selections_group()
    vis.additional_cleanup()
    if config.PYMOL:
        vis.save_session(config.OUTPATH)
    if config.PICS:
        filename = '%s_%s' % (pdbid.upper(), "_".join([ligname, plcomplex.chain, plcomplex.position]))
        vis.save_picture(config.OUTPATH, filename)
