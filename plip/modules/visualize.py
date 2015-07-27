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
from supplemental import *
from time import sleep


def set_fancy_ray():
    """Give the molecule a flat, modern look."""
    cmd.set('light_count', 6)
    cmd.set('spec_count', 1.5)
    cmd.set('shininess', 4)
    cmd.set('specular', 0.3)
    cmd.set('reflect', 1.6)
    cmd.set('ambient', 0)
    cmd.set('direct', 0)
    cmd.set('ray_shadow', 0)  # Gives the molecules a flat, modern look
    cmd.set('ambient_occlusion_mode', 1)


def png_workaround(filepath, width=1200, height=800):
    """Workaround for (a) severe bug(s) in PyMOL preventing ray-traced images to be produced in command-line mode.
    Use this function in case neither cmd.ray() or cmd.png() work.
    """
    sys.stdout = sys.__stdout__
    cmd.feedback('disable', 'movie', 'everything')
    cmd.viewport(width, height)
    cmd.zoom('visible', 1.5)  # Adapt the zoom to the viewport
    cmd.set('ray_trace_frames', 1)  # Frames are raytraced before saving an image.
    cmd.mpng(filepath, 1, 1)  # Use batch png mode with 1 frame only
    cmd.mplay()  # cmd.mpng needs the animation to 'run'
    cmd.refresh()
    originalfile = "".join([filepath, '0001.png'])
    newfile = "".join([filepath, '.png'])

    #################################################
    # Wait for file for max. 1 second and rename it #
    #################################################

    attempts = 0
    while not os.path.isfile(originalfile) and attempts <= 10:
        sleep(0.1)
        attempts += 1
    if os.name == 'nt':  # In Windows, make sure there is no file of the same name, cannot be overwritten as in Unix
        if os.path.isfile(newfile):
            os.remove(newfile)
    os.rename(originalfile, newfile)  # Remove frame number in filename

    #  Check if imagemagick is available and crop + resize the images
    if cmd_exists('convert'):
        attempts, ecode = 0, 1
        # Check if file is truncated and wait if that's the case
        while ecode != 0 and attempts <= 10:
            ecode = subprocess.call(['convert', newfile, '/dev/null'], stdout=open('/dev/null', 'w'),
                                    stderr=subprocess.STDOUT)
            sleep(0.1)
            attempts += 1
        trim = 'convert -trim ' + newfile + ' -bordercolor White -border 20x20 ' + newfile + ';'  # Trim the image
        os.system(trim)
        getwidth = 'w=`convert ' + newfile + ' -ping -format "%w" info:`;'  # Get the width of the new image
        getheight = 'h=`convert ' + newfile + ' -ping -format "%h" info:`;'  # Get the hight of the new image
        newres = 'if [ "$w" -gt "$h" ]; then newr="${w%.*}x$w"; else newr="${h%.*}x$h"; fi;'  # Set quadratic ratio
        quadratic = 'convert ' + newfile + ' -gravity center -extent "$newr" ' + newfile  # Fill with whitespace
        os.system(getwidth + getheight + newres + quadratic)
    else:
        sys.stderr.write('Imagemagick not available. Images will not be resized or cropped.')


def visualize_in_pymol(protcomplex_class, pli_site):
    """Visualizes the protein-ligand pliprofiler at one site in PyMOL."""

    #####################
    # Set everything up #
    #####################

    pcomp = protcomplex_class
    pdbid = pcomp.pymol_name
    pli = pcomp.interaction_sets[pli_site]  # Select the interaction class corresponding to the selection
    ligand = pli.ligand
    lig_members = sorted(pli.ligand.members)
    save_to = protcomplex_class.output_path
    chain = ligand.chain if not ligand.chain == "0" else ""
    ligname = ligand.hetid
    metal_ids = [x.m_orig_idx for x in pli.ligand.metals]
    metal_ids_str = '+'.join([str(i) for i in metal_ids])

    ########################
    # Basic visualizations #
    ########################

    start_pymol(run=True, options='-pcq', quiet=True)
    standard_settings()
    cmd.set('dash_gap', 0)  # Show not dashes, but lines for the pliprofiler
    cmd.set('ray_shadow', 0)  # Turn on ray shadows for clearer ray-traced images
    cmd.set('cartoon_color', 'mylightblue')
    cmd.load(pcomp.sourcefiles['pdbcomplex'])
    current_name = cmd.get_object_list(selection='(all)')[0]
    cmd.set_name(current_name, pdbid)
    cmd.hide('everything', 'all')
    cmd.select(ligname, 'resn %s and chain %s and resi %s' % (ligand.hetid, chain, ligand.position))

    # Visualize and color metal ions if there are any
    if not len(metal_ids) == 0:
        cmd.select(ligname, '%s or id %s' % (ligname, metal_ids_str))
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

    ###########################
    # Create empty selections #
    ###########################

    for group in ['Hydrophobic-P', 'Hydrophobic-L', 'HBondDonor-P', 'HBondDonor-L', 'HBondAccept-P', 'HBondAccept-L',
                  'HalogenAccept', 'HalogenDonor', 'Water', 'MetalIons', 'StackRings-P', 'PosCharge-P', 'PosCharge-L',
                  'NegCharge-P', 'NegCharge-L', 'PiCatRing-P', 'StackRings-L', 'PiCatRing-L', 'Metal-M', 'Metal-P',
                  'Metal-W', 'Metal-L', 'Unpaired-HBA', 'Unpaired-HBD', 'Unpaired-HAL', 'Unpaired-RINGS']:
        cmd.select(group, 'None')

    ######################################
    # Visualize hydrophobic interactions #
    ######################################

    # #@todo Dont use indices!
    hydroph = pli.hydrophobic_contacts
    print hydroph
    hydroph_pairs_id = [(h.bsatom_orig_idx, h.ligatom_orig_idx) for h in hydroph]
    hydroph_bs_id, hydroph_lig_id = [hp[0] for hp in hydroph_pairs_id], [hp[1] for hp in hydroph_pairs_id]
    if not len(hydroph_bs_id) == 0:
        for h in [['Hydrophobic-P', hydroph_bs_id], ['Hydrophobic-L', hydroph_lig_id]]:
            cmd.select(h[0], 'id %s' % '+'.join(map(str, h[1])))
        for i in hydroph_pairs_id:
            cmd.select('tmp_bs', 'id %i' % i[0])
            cmd.select('tmp_lig', 'id %i' % i[1])
            cmd.distance('Hydrophobic', 'tmp_bs', 'tmp_lig')
        if object_exists('Hydrophobic'):
            cmd.set('dash_gap', 0.5, 'Hydrophobic')
            cmd.set('dash_color', 'grey50', 'Hydrophobic')
    else:
        cmd.select('Hydrophobic-P', 'None')

    #####################
    # Visualize H-Bonds #
    #####################

    hbonds_ldon, hbonds_pdon = pli.hbonds_ldon, pli.hbonds_pdon
    hbonds_ldon_id = [(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_ldon]
    hbonds_lig_don_id, hbonds_prot_acc_id = [hb[1] for hb in hbonds_ldon_id], [hb[0] for hb in hbonds_ldon_id]
    hbonds_pdon_id = [(hb.a_orig_idx, hb.d_orig_idx) for hb in hbonds_pdon]
    hbonds_prot_don_id, hbonds_lig_acc_id = [hb[1] for hb in hbonds_pdon_id], [hb[0] for hb in hbonds_pdon_id]
    for group in [['HBondDonor-L', hbonds_lig_don_id], ['HBondDonor-P', hbonds_prot_don_id],
                  ['HBondAccept-L', hbonds_lig_acc_id], ['HBondAccept-P', hbonds_prot_acc_id]]:
        if not len(group[1]) == 0:
            cmd.select(group[0], 'id %s' % '+'.join(map(str, group[1])))
    for i in hbonds_ldon_id:
        cmd.select('tmp_bs', 'id %i' % i[0])
        cmd.select('tmp_lig', 'id %i' % i[1])
        cmd.distance('HBonds', 'tmp_bs', 'tmp_lig')
    for i in hbonds_pdon_id:
        cmd.select('tmp_bs', 'id %i' % i[1])
        cmd.select('tmp_lig', 'id %i' % i[0])
        cmd.distance('HBonds', 'tmp_bs', 'tmp_lig')
    if object_exists('HBonds'):
        cmd.set('dash_color', 'blue', 'HBonds')

    ###########################
    # Visualize Halogen Bonds #
    ###########################

    halogen = pli.halogen_bonds
    all_don_x, all_acc_o = [], []
    for h in halogen:
        all_don_x.append(h.don_orig_idx)
        all_acc_o.append(h.acc_orig_idx)
        for group in [['tmp_bs', h.acc_orig_idx], ['tmp_lig', h.don_orig_idx]]:
            cmd.select(group[0], 'id %i' % group[1])
        cmd.distance('HalogenBonds', 'tmp_bs', 'tmp_lig')
    if not len(all_acc_o) == 0:
        cmd.select('HalogenAccept', 'id %s' % '+'.join(map(str, all_acc_o)))
        cmd.select('HalogenDonor', 'id %s' % '+'.join(map(str, all_don_x)))
    if object_exists('HalogenBonds'):
        cmd.set('dash_color', 'greencyan', 'HalogenBonds')

    #########################
    # Visualize Pi-Stacking #
    #########################

    stacks = pli.pistacking
    for i, stack in enumerate(stacks):
        pires_ids = '+'.join(map(str, stack.proteinring.atoms_orig_idx))
        pilig_ids = '+'.join(map(str, stack.ligandring.atoms_orig_idx))
        cmd.select('StackRings-P', 'StackRings-P or id %s' % pires_ids)
        cmd.select('StackRings-L', 'StackRings-L or id %s' % pilig_ids)
        cmd.select('StackRings-P', 'byres StackRings-P')
        cmd.show('sticks', 'StackRings-P')

        cmd.pseudoatom('ps-pistack-1-%i' % i, pos=stack.proteinring.center)
        cmd.pseudoatom('ps-pistack-2-%i' % i, pos=stack.ligandring.center)
        cmd.pseudoatom('Centroids-P', pos=stack.proteinring.center)
        cmd.pseudoatom('Centroids-L', pos=stack.ligandring.center)

        if stack.type == 'P':
            cmd.distance('PiStackingP', 'ps-pistack-1-%i' % i, 'ps-pistack-2-%i' % i)
        if stack.type == 'T':
            cmd.distance('PiStackingT', 'ps-pistack-1-%i' % i, 'ps-pistack-2-%i' % i)
    if object_exists('PiStackingP'):
        cmd.set('dash_color', 'green', 'PiStackingP')
        cmd.set('dash_gap', 0.3, 'PiStackingP')
        cmd.set('dash_length', 0.6, 'PiStackingP')
    if object_exists('PiStackingT'):
        cmd.set('dash_color', 'smudge', 'PiStackingT')
        cmd.set('dash_gap', 0.3, 'PiStackingT')
        cmd.set('dash_length', 0.6, 'PiStackingT')

    ####################################
    # Visualize Cation-pi interactions #
    ####################################

    for i, p in enumerate(pli.pication_paro+pli.pication_laro):
        cmd.pseudoatom('ps-picat-1-%i' % i, pos=p.ring.center)
        cmd.pseudoatom('ps-picat-2-%i' % i, pos=p.charge.center)
        if p.protcharged:
            cmd.pseudoatom('Chargecenter-P', pos=p.charge.center)
            cmd.pseudoatom('Centroids-L', pos=p.ring.center)
            pilig_ids = '+'.join(map(str, p.ring.atoms_orig_idx))
            cmd.select('PiCatRing-L', 'PiCatRing-L or id %s' % pilig_ids)
            for a in p.charge.atoms_orig_idx:
                cmd.select('PosCharge-P', 'PosCharge-P or id %i' % a)
        else:
            cmd.pseudoatom('Chargecenter-L', pos=p.charge.center)
            cmd.pseudoatom('Centroids-P', pos=p.ring.center)
            pires_ids = '+'.join(map(str, p.ring.atoms_orig_idx))
            cmd.select('PiCatRing-P', 'PiCatRing-P or id %s' % pires_ids)
            for a in p.charge.atoms_orig_idx:
                cmd.select('PosCharge-L', 'PosCharge-L or id %i' % a)
        cmd.distance('PiCation', 'ps-picat-1-%i' % i, 'ps-picat-2-%i' % i)
    if object_exists('PiCation'):
        cmd.set('dash_color', 'orange', 'PiCation')
        cmd.set('dash_gap', 0.3, 'PiCation')
        cmd.set('dash_length', 0.6, 'PiCation')

    ##########################
    # Visualize salt bridges #
    ##########################

    saltbridge_lneg = pli.saltbridge_lneg
    saltbridge_pneg = pli.saltbridge_pneg

    for i, saltb in enumerate(saltbridge_lneg):
        for patom in saltb.positive.atoms_orig_idx:
            cmd.select('PosCharge-P', 'PosCharge-P or id %i' % patom)
        for latom in saltb.negative.atoms_orig_idx:
            cmd.select('NegCharge-L', 'NegCharge-L or id %i' % latom)
        for sbgroup in [['ps-sbl-1-%i' % i, 'Chargecenter-P', saltb.positive.center],
                        ['ps-sbl-2-%i' % i, 'Chargecenter-L', saltb.negative.center]]:
            cmd.pseudoatom(sbgroup[0], pos=sbgroup[2])
            cmd.pseudoatom(sbgroup[1], pos=sbgroup[2])
        cmd.distance('Saltbridges', 'ps-sbl-1-%i' % i, 'ps-sbl-2-%i' % i)

    for i, saltb in enumerate(saltbridge_pneg):
        for patom in saltb.negative.atoms_orig_idx:
            cmd.select('NegCharge-P', 'NegCharge-P or id %i' % patom)
        for latom in saltb.positive.atoms_orig_idx:
            cmd.select('PosCharge-L', 'PosCharge-L or id %i' % latom)
        for sbgroup in [['ps-sbp-1-%i' % i, 'Chargecenter-P', saltb.negative.center],
                        ['ps-sbp-2-%i' % i, 'Chargecenter-L', saltb.positive.center]]:
            cmd.pseudoatom(sbgroup[0], pos=sbgroup[2])
            cmd.pseudoatom(sbgroup[1], pos=sbgroup[2])
        cmd.distance('Saltbridges', 'ps-sbp-1-%i' % i, 'ps-sbp-2-%i' % i)

    if object_exists('Saltbridges'):
        cmd.set('dash_color', 'yellow', 'Saltbridges')
        cmd.set('dash_gap', 0.5, 'Saltbridges')

    ########################################
    # Water-bridged H-Bonds (first degree) #
    ########################################

    for bridge in pli.water_bridges:
        if bridge.protisdon:
            cmd.select('HBondDonor-P', 'HBondDonor-P or id %i' % bridge.d_orig_idx)
            cmd.select('HBondAccept-L', 'HBondAccept-L or id %i' % bridge.a_orig_idx)
        else:
            cmd.select('HBondDonor-L', 'HBondDonor-L or id %i' % bridge.d_orig_idx)
            cmd.select('HBondAccept-P', 'HBondAccept-P or id %i' % bridge.a_orig_idx)
        cmd.select('Water', 'Water or id %i' % bridge.water_orig_idx)
        cmd.select('tmp_don', 'id %i' % bridge.d_orig_idx)
        cmd.select('tmp_water', 'id %i' % bridge.water_orig_idx)
        cmd.select('tmp_acc', 'id %i' % bridge.a_orig_idx)
        cmd.distance('WaterBridges', 'tmp_acc', 'tmp_water')
        cmd.distance('WaterBridges', 'tmp_don', 'tmp_water')
    if object_exists('WaterBridges'):
        cmd.set('dash_color', 'lightblue', 'WaterBridges')
    cmd.delete('tmp_water or tmp_acc or tmp_don')
    cmd.color('lightblue', 'Water')
    cmd.show('spheres', 'Water')

    ###################
    # Metal Complexes #
    ###################

    if not len(pli.metal_complexes) == 0:
        cmd.select('Metal-M', 'id %s' % metal_ids_str)
        for metal_complex in pli.metal_complexes:
            cmd.select('tmp_m', 'id %i' % metal_complex.metal_orig_idx)
            cmd.select('tmp_t', 'id %i' % metal_complex.target_orig_idx)
            if metal_complex.location == 'water':
                cmd.select('Metal-W', 'Metal-W or id %s' % metal_complex.target_orig_idx)
            if metal_complex.location.startswith('protein'):
                cmd.select('Metal-P', 'Metal-P or id %s' % metal_complex.target_orig_idx)
            if metal_complex.location == 'ligand':
                cmd.select('Metal-L', 'Metal-L or id %s' % metal_complex.target_orig_idx)
            cmd.distance('MetalComplexes', 'tmp_m', 'tmp_t')
            cmd.delete('tmp_m or tmp_t')
    if object_exists('MetalComplexes'):
        cmd.set('dash_color', 'violetpurple', 'MetalComplexes')
        cmd.set('dash_gap', 0.5, 'MetalComplexes')
        # Show water molecules for metal complexes
        cmd.show('spheres', 'Metal-W')
        cmd.color('lightblue', 'Metal-W')

    ######################
    # Visualize the rest #
    ######################

    # Show sticks for all residues interacing with the ligand
    cmd.select('AllBSRes', 'byres (Hydrophobic-P or HBondDonor-P or HBondAccept-P or PosCharge-P or NegCharge-P or '
                           'StackRings-P or PiCatRing-P or HalogenAcc or Metal-P)')
    cmd.show('sticks', 'AllBSRes')
    # Show spheres for the ring centroids
    cmd.hide('everything', 'centroids*')
    cmd.show('nb_spheres', 'centroids*')
    # Show spheres for centers of charge
    if object_exists('Chargecenter-P') or object_exists('Chargecenter-L'):
        cmd.hide('nonbonded', 'chargecenter*')
        cmd.show('spheres', 'chargecenter*')
        cmd.set('sphere_scale', 0.4, 'chargecenter*')
        cmd.color('yellow', 'chargecenter*')

    ####################
    # Last refinements #
    ####################

    cmd.set('valence', 1)  # Show bond valency (e.g. double bonds)
    # Optional cartoon representation of the protein
    cmd.copy('%sCartoon' % pdbid, pdbid)
    cmd.show('cartoon', '%sCartoon' % pdbid)
    cmd.show('sticks', '%sCartoon' % pdbid)
    cmd.set('stick_transparency', 1, '%sCartoon' % pdbid)
    # Set view. Zoom on the ligand (and its pliprofiler)
    if 'AllBSRes' in cmd.get_names("selections"):
        cmd.zoom('%s or AllBSRes' % ligname)
    else:
        if object_exists(ligname):
            cmd.zoom(ligname, 3)
    cmd.set('sphere_scale', 0.2, 'resn HOH')  # Needs to be done here because of the copy made
    cmd.set('sphere_transparency', 0.4, '!resn HOH')
    cmd.origin(ligname)
    if 'Centroids*' in cmd.get_names("selections"):
        cmd.color('grey80', 'Centroids*')
    cmd.hide('spheres', '%sCartoon' % pdbid)
    cmd.hide('cartoon', '%sCartoon and resn DA+DG+DC+DU+DT+A+G+C+U+T' % pdbid)  # Hide DNA/RNA Cartoon
    if ligname == 'SF4':  # Special case for iron-sulfur clusters, can't be visualized with sticks
        cmd.show('spheres', '%s' % ligname)

    ##################################
    # Selections for unpaired groups #
    ##################################
    if not len(pli.unpaired_hba) == 0:
        cmd.select('Unpaired-HBA', 'Unpaired-HBA or id %s' % '+'.join(str(idx) for idx in pli.unpaired_hba_orig_idx))
    if not len(pli.unpaired_hbd) == 0:
        cmd.select('Unpaired-HBD', 'Unpaired-HBD or id %s' % '+'.join(str(idx) for idx in pli.unpaired_hbd_orig_idx))
    if not len(pli.unpaired_hal) == 0:
        cmd.select('Unpaired-HAL', 'Unpaired-HAL or id %s' % '+'.join(str(idx) for idx in pli.unpaired_hal_orig_idx))

    ##############################
    # Organization of selections #
    ##############################

    # Delete all empty and temporary selections
    selections = cmd.get_names("selections")
    for selection in selections:
        if len(cmd.get_model(selection).atom) == 0:
            cmd.delete(selection)
    cmd.deselect()
    cmd.delete('tmp*')
    cmd.delete('ps-*')

    # Group non-empty selections
    cmd.group('Structures', '%s %s %sCartoon' % (pdbid, ligname, pdbid))
    cmd.group('Interactions', 'Hydrophobic HBonds HalogenBonds WaterBridges PiCation PiStackingP PiStackingT '
                              'Saltbridges MetalComplexes')
    cmd.group('Atoms', '')
    cmd.group('Atoms.Protein', 'Hydrophobic-P HBondAccept-P HBondDonor-P HalogenAccept Centroids-P PiCatRing-P '
                               'StackRings-P PosCharge-P NegCharge-P AllBSRes Chargecenter-P  Metal-P')
    cmd.group('Atoms.Ligand', 'Hydrophobic-L HBondAccept-L HBondDonor-L HalogenDonor Centroids-L NegCharge-L '
                              'PosCharge-L NegCharge-L ChargeCenter-L StackRings-L PiCatRing-L Metal-L Metal-M '
                              'Unpaired-HBA Unpaired-HBD Unpaired-HAL Unpaired-RINGS')
    cmd.group('Atoms.Other', 'Water Metal-W')
    cmd.order('*', 'y')

    ###############################################
    # Remove atoms with alternative conformations #
    ###############################################

    cmd.remove('not alt ""+A')

    ########################################
    # Clean up and save PyMOL session file #
    ########################################

    cmd.hide('labels', 'Interactions')
    cmd.disable('%sCartoon' % pdbid)
    cmd.hide('everything', 'hydrogens')

    filename = '%s_%s' % (pdbid.upper(), "_".join([ligand.hetid, ligand.chain, str(ligand.position)]))
    if config.PYMOL:
        cmd.save("".join([save_to, "%s.pse" % filename]))

    # Create output pictures (experimental)
    set_fancy_ray()
    if config.PICS:
        png_workaround("".join([save_to, filename]))
