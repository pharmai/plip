"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
visualize.py - Visualization of PLIP results using PyMOL.
"""

# Python standard library
from __future__ import absolute_import

# Own modules
from .supplemental import initialize_pymol, start_pymol, write_message, colorlog, sysexit
from . import config
from .pymolplip import PyMOLVisualizer
from .plipremote import VisualizerData

# Python Standard Library
import json
import sys

# Special imports
from pymol import cmd
import pymol


def select_by_ids(selname, idlist, selection_exists=False, chunksize=20, restrict=None):
    """Selection with a large number of ids concatenated into a selection
    list can cause buffer overflow in PyMOL. This function takes a selection
    name and and list of IDs (list of integers) as input and makes a careful
    step-by-step selection (packages of 20 by default)"""
    idlist = list(set(idlist))  # Remove duplicates
    if not selection_exists:
        cmd.select(selname, 'None')  # Empty selection first
    idchunks = [idlist[i:i+chunksize] for i in range(0, len(idlist), chunksize)]
    for idchunk in idchunks:
        cmd.select(selname, '%s or (id %s)' % (selname, '+'.join(map(str, idchunk))))
    if restrict is not None:
        cmd.select(selname, '%s and %s' % (selname, restrict))


def visualize_in_pymol(plcomplex):
    """Visualizes the protein-ligand pliprofiler at one site in PyMOL."""

    vis = PyMOLVisualizer(plcomplex)



    #####################
    # Set everything up #
    #####################

    pdbid = plcomplex.pdbid
    lig_members = plcomplex.lig_members
    chain = plcomplex.chain
    if config.PEPTIDES != []:
        vis.ligname = 'PeptideChain%s' % plcomplex.chain
    if config.INTRA is not None:
        vis.ligname = 'Intra%s' % plcomplex.chain

    ligname = vis.ligname
    hetid = plcomplex.hetid

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
    if config.PEPTIDES != []:
        cmd.select(ligname, 'chain %s and not resn HOH' % plcomplex.chain)
    else:
        cmd.select(ligname, 'resn %s and chain %s and resi %s*' % (hetid, chain, plcomplex.position))
    write_message("Selecting ligand for PDBID %s and ligand name %s with: " % (pdbid, ligname), mtype='debug')
    write_message('resn %s and chain %s and resi %s*' % (hetid, chain, plcomplex.position), mtype='debug')

    # Visualize and color metal ions if there are any
    if not len(metal_ids) == 0:
        vis.select_by_ids(ligname, metal_ids, selection_exists=True)
        cmd.show('spheres', 'id %s and %s' % (metal_ids_str, pdbid))

    # Additionally, select all members of composite ligands
    if len(lig_members) > 1:
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
    if config.DNARECEPTOR:
        # Rename Cartoon selection to Line selection and change repr.
        cmd.set_name('%sCartoon' % plcomplex.pdbid, '%sLines' % plcomplex.pdbid)
        cmd.hide('cartoon', '%sLines' % plcomplex.pdbid)
        cmd.show('lines', '%sLines' % plcomplex.pdbid)

    if config.PEPTIDES != []:
        filename = "%s_PeptideChain%s" % (pdbid.upper(), plcomplex.chain)
        if config.PYMOL:
            vis.save_session(config.OUTPATH, override=filename)
    elif config.INTRA is not None:
        filename = "%s_IntraChain%s" % (pdbid.upper(), plcomplex.chain)
        if config.PYMOL:
            vis.save_session(config.OUTPATH, override=filename)
    else:
        filename = '%s_%s' % (pdbid.upper(), "_".join([hetid, plcomplex.chain, plcomplex.position]))
        if config.PYMOL:
            vis.save_session(config.OUTPATH)
    if config.PICS:
        vis.save_picture(config.OUTPATH, filename)
