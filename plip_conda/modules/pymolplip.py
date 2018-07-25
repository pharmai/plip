"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
pymolplip.py - Visualization class for PyMOL.
"""

from pymol import cmd
from time import sleep
import sys
import os
import subprocess

class PyMOLVisualizer:

    def __init__(self, plcomplex):
        if not plcomplex is None:
            self.plcomplex = plcomplex
            self.protname = plcomplex.pdbid  # Name of protein with binding site
            self.hetid = plcomplex.hetid
            self.ligandtype = plcomplex.ligandtype
            self.ligname = "Ligand_" + self.hetid  # Name of ligand
            self.metal_ids = plcomplex.metal_ids

    def set_initial_representations(self):
        """General settings for PyMOL"""
        self.standard_settings()
        cmd.set('dash_gap', 0)  # Show not dashes, but lines for the pliprofiler
        cmd.set('ray_shadow', 0)  # Turn on ray shadows for clearer ray-traced images
        cmd.set('cartoon_color', 'mylightblue')

        # Set clipping planes for full view
        cmd.clip('far', -1000)
        cmd.clip('near', 1000)

    def make_initial_selections(self):
        """Make empty selections for structures and interactions"""
        for group in ['Hydrophobic-P', 'Hydrophobic-L', 'HBondDonor-P',
        'HBondDonor-L', 'HBondAccept-P', 'HBondAccept-L',
        'HalogenAccept', 'HalogenDonor', 'Water', 'MetalIons', 'StackRings-P',
        'PosCharge-P', 'PosCharge-L', 'NegCharge-P', 'NegCharge-L',
        'PiCatRing-P', 'StackRings-L', 'PiCatRing-L', 'Metal-M', 'Metal-P',
        'Metal-W', 'Metal-L', 'Unpaired-HBA', 'Unpaired-HBD', 'Unpaired-HAL',
        'Unpaired-RINGS']:
            cmd.select(group, 'None')

    def standard_settings(self):
        """Sets up standard settings for a nice visualization."""
        cmd.set('bg_rgb', [1.0, 1.0, 1.0])  # White background
        cmd.set('depth_cue', 0)  # Turn off depth cueing (no fog)
        cmd.set('cartoon_side_chain_helper', 1)  # Improve combined visualization of sticks and cartoon
        cmd.set('cartoon_fancy_helices', 1)  # Nicer visualization of helices (using tapered ends)
        cmd.set('transparency_mode', 1)  # Turn on multilayer transparency
        cmd.set('dash_radius', 0.05)
        self.set_custom_colorset()

    def set_custom_colorset(self):
        """Defines a colorset with matching colors. Provided by Joachim."""
        cmd.set_color('myorange', '[253, 174, 97]')
        cmd.set_color('mygreen', '[171, 221, 164]')
        cmd.set_color('myred', '[215, 25, 28]')
        cmd.set_color('myblue', '[43, 131, 186]')
        cmd.set_color('mylightblue', '[158, 202, 225]')
        cmd.set_color('mylightgreen', '[229, 245, 224]')

    def select_by_ids(self, selname, idlist, selection_exists=False, chunksize=20, restrict=None):
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


    def object_exists(self, object_name):
        """Checks if an object exists in the open PyMOL session."""
        return object_name in cmd.get_names("objects")

    def show_hydrophobic(self):
        """Visualizes hydrophobic contacts."""
        hydroph = self.plcomplex.hydrophobic_contacts
        if not len(hydroph.bs_ids) == 0:
            self.select_by_ids('Hydrophobic-P', hydroph.bs_ids, restrict=self.protname)
            self.select_by_ids('Hydrophobic-L', hydroph.lig_ids, restrict=self.ligname)

            for i in hydroph.pairs_ids:
                cmd.select('tmp_bs', 'id %i & %s' % (i[0], self.protname))
                cmd.select('tmp_lig', 'id %i & %s' % (i[1], self.ligname))
                cmd.distance('Hydrophobic', 'tmp_bs', 'tmp_lig')
            if self.object_exists('Hydrophobic'):
                cmd.set('dash_gap', 0.5, 'Hydrophobic')
                cmd.set('dash_color', 'grey50', 'Hydrophobic')
        else:
            cmd.select('Hydrophobic-P', 'None')

    def show_hbonds(self):
        """Visualizes hydrogen bonds."""
        hbonds = self.plcomplex.hbonds
        for group in [['HBondDonor-P', hbonds.prot_don_id],
        ['HBondAccept-P', hbonds.prot_acc_id]]:
            if not len(group[1]) == 0:
                self.select_by_ids(group[0], group[1], restrict=self.protname)
        for group in [['HBondDonor-L', hbonds.lig_don_id],
        ['HBondAccept-L', hbonds.lig_acc_id]]:
            if not len(group[1]) == 0:
                self.select_by_ids(group[0], group[1], restrict=self.ligname)
        for i in hbonds.ldon_id:
            cmd.select('tmp_bs', 'id %i & %s' % (i[0], self.protname))
            cmd.select('tmp_lig', 'id %i & %s' % (i[1], self.ligname))
            cmd.distance('HBonds', 'tmp_bs', 'tmp_lig')
        for i in hbonds.pdon_id:
            cmd.select('tmp_bs', 'id %i & %s' % (i[1], self.protname))
            cmd.select('tmp_lig', 'id %i & %s' % (i[0], self.ligname))
            cmd.distance('HBonds', 'tmp_bs', 'tmp_lig')
        if self.object_exists('HBonds'):
            cmd.set('dash_color', 'blue', 'HBonds')

    def show_halogen(self):
        """Visualize halogen bonds."""
        halogen = self.plcomplex.halogen_bonds
        all_don_x, all_acc_o = [], []
        for h in halogen:
            all_don_x.append(h.don_id)
            all_acc_o.append(h.acc_id)
            cmd.select('tmp_bs', 'id %i & %s' % (h.acc_id, self.protname))
            cmd.select('tmp_lig', 'id %i & %s' % (h.don_id, self.ligname))

            cmd.distance('HalogenBonds', 'tmp_bs', 'tmp_lig')
        if not len(all_acc_o) == 0:
            self.select_by_ids('HalogenAccept', all_acc_o, restrict=self.protname)
            self.select_by_ids('HalogenDonor', all_don_x, restrict=self.ligname)
        if self.object_exists('HalogenBonds'):
            cmd.set('dash_color', 'greencyan', 'HalogenBonds')

    def show_stacking(self):
        """Visualize pi-stacking interactions."""
        stacks = self.plcomplex.pistacking
        for i, stack in enumerate(stacks):
            pires_ids = '+'.join(map(str, stack.proteinring_atoms))
            pilig_ids = '+'.join(map(str, stack.ligandring_atoms))
            cmd.select('StackRings-P', 'StackRings-P or (id %s & %s)' % (pires_ids, self.protname))
            cmd.select('StackRings-L', 'StackRings-L or (id %s & %s)' % (pilig_ids, self.ligname))
            cmd.select('StackRings-P', 'byres StackRings-P')
            cmd.show('sticks', 'StackRings-P')

            cmd.pseudoatom('ps-pistack-1-%i' % i, pos=stack.proteinring_center)
            cmd.pseudoatom('ps-pistack-2-%i' % i, pos=stack.ligandring_center)
            cmd.pseudoatom('Centroids-P', pos=stack.proteinring_center)
            cmd.pseudoatom('Centroids-L', pos=stack.ligandring_center)

            if stack.type == 'P':
                cmd.distance('PiStackingP', 'ps-pistack-1-%i' % i, 'ps-pistack-2-%i' % i)
            if stack.type == 'T':
                cmd.distance('PiStackingT', 'ps-pistack-1-%i' % i, 'ps-pistack-2-%i' % i)
        if self.object_exists('PiStackingP'):
            cmd.set('dash_color', 'green', 'PiStackingP')
            cmd.set('dash_gap', 0.3, 'PiStackingP')
            cmd.set('dash_length', 0.6, 'PiStackingP')
        if self.object_exists('PiStackingT'):
            cmd.set('dash_color', 'smudge', 'PiStackingT')
            cmd.set('dash_gap', 0.3, 'PiStackingT')
            cmd.set('dash_length', 0.6, 'PiStackingT')

    def show_cationpi(self):
        """Visualize cation-pi interactions."""
        for i, p in enumerate(self.plcomplex.pication):
            cmd.pseudoatom('ps-picat-1-%i' % i, pos=p.ring_center)
            cmd.pseudoatom('ps-picat-2-%i' % i, pos=p.charge_center)
            if p.protcharged:
                cmd.pseudoatom('Chargecenter-P', pos=p.charge_center)
                cmd.pseudoatom('Centroids-L', pos=p.ring_center)
                pilig_ids = '+'.join(map(str, p.ring_atoms))
                cmd.select('PiCatRing-L', 'PiCatRing-L or (id %s & %s)' % (pilig_ids, self.ligname))
                for a in p.charge_atoms:
                    cmd.select('PosCharge-P', 'PosCharge-P or (id %i & %s)' % (a, self.protname))
            else:
                cmd.pseudoatom('Chargecenter-L', pos=p.charge_center)
                cmd.pseudoatom('Centroids-P', pos=p.ring_center)
                pires_ids = '+'.join(map(str, p.ring_atoms))
                cmd.select('PiCatRing-P', 'PiCatRing-P or (id %s & %s)' % (pires_ids, self.protname))
                for a in p.charge_atoms:
                    cmd.select('PosCharge-L', 'PosCharge-L or (id %i & %s)' % (a, self.ligname))
            cmd.distance('PiCation', 'ps-picat-1-%i' % i, 'ps-picat-2-%i' % i)
        if self.object_exists('PiCation'):
            cmd.set('dash_color', 'orange', 'PiCation')
            cmd.set('dash_gap', 0.3, 'PiCation')
            cmd.set('dash_length', 0.6, 'PiCation')

    def show_sbridges(self):
        """Visualize salt bridges."""
        for i, saltb in enumerate(self.plcomplex.saltbridges):
            if saltb.protispos:
                for patom in saltb.positive_atoms:
                    cmd.select('PosCharge-P', 'PosCharge-P or (id %i & %s)' % (patom, self.protname))
                for latom in saltb.negative_atoms:
                    cmd.select('NegCharge-L', 'NegCharge-L or (id %i & %s)' % (latom, self.ligname))
                for sbgroup in [['ps-sbl-1-%i' % i, 'Chargecenter-P', saltb.positive_center],
                                ['ps-sbl-2-%i' % i, 'Chargecenter-L', saltb.negative_center]]:
                    cmd.pseudoatom(sbgroup[0], pos=sbgroup[2])
                    cmd.pseudoatom(sbgroup[1], pos=sbgroup[2])
                cmd.distance('Saltbridges', 'ps-sbl-1-%i' % i, 'ps-sbl-2-%i' % i)
            else:
                for patom in saltb.negative_atoms:
                    cmd.select('NegCharge-P', 'NegCharge-P or (id %i & %s)' % (patom, self.protname))
                for latom in saltb.positive_atoms:
                    cmd.select('PosCharge-L', 'PosCharge-L or (id %i & %s)' % (latom, self.ligname))
                for sbgroup in [['ps-sbp-1-%i' % i, 'Chargecenter-P', saltb.negative_center],
                                ['ps-sbp-2-%i' % i, 'Chargecenter-L', saltb.positive_center]]:
                    cmd.pseudoatom(sbgroup[0], pos=sbgroup[2])
                    cmd.pseudoatom(sbgroup[1], pos=sbgroup[2])
                cmd.distance('Saltbridges', 'ps-sbp-1-%i' % i, 'ps-sbp-2-%i' % i)

        if self.object_exists('Saltbridges'):
            cmd.set('dash_color', 'yellow', 'Saltbridges')
            cmd.set('dash_gap', 0.5, 'Saltbridges')

    def show_wbridges(self):
        """Visualize water bridges."""
        for bridge in self.plcomplex.waterbridges:
            if bridge.protisdon:
                cmd.select('HBondDonor-P', 'HBondDonor-P or (id %i & %s)' % (bridge.don_id, self.protname))
                cmd.select('HBondAccept-L', 'HBondAccept-L or (id %i & %s)' % (bridge.acc_id, self.ligname))
                cmd.select('tmp_don', 'id %i & %s' % (bridge.don_id, self.protname))
                cmd.select('tmp_acc', 'id %i & %s' % (bridge.acc_id, self.ligname))
            else:
                cmd.select('HBondDonor-L', 'HBondDonor-L or (id %i & %s)' % (bridge.don_id, self.ligname))
                cmd.select('HBondAccept-P', 'HBondAccept-P or (id %i & %s)' % (bridge.acc_id, self.protname))
                cmd.select('tmp_don', 'id %i & %s' % (bridge.don_id, self.ligname))
                cmd.select('tmp_acc', 'id %i & %s' % (bridge.acc_id, self.protname))
            cmd.select('Water', 'Water or (id %i & resn HOH)' % bridge.water_id)
            cmd.select('tmp_water', 'id %i & resn HOH' % bridge.water_id)
            cmd.distance('WaterBridges', 'tmp_acc', 'tmp_water')
            cmd.distance('WaterBridges', 'tmp_don', 'tmp_water')
        if self.object_exists('WaterBridges'):
            cmd.set('dash_color', 'lightblue', 'WaterBridges')
        cmd.delete('tmp_water or tmp_acc or tmp_don')
        cmd.color('lightblue', 'Water')
        cmd.show('spheres', 'Water')

    def show_metal(self):
        """Visualize metal coordination."""
        metal_complexes = self.plcomplex.metal_complexes
        if not len(metal_complexes) == 0:
            self.select_by_ids('Metal-M', self.metal_ids)
            for metal_complex in metal_complexes:
                cmd.select('tmp_m', 'id %i' % metal_complex.metal_id)
                cmd.select('tmp_t', 'id %i' % metal_complex.target_id)
                if metal_complex.location == 'water':
                    cmd.select('Metal-W', 'Metal-W or id %s' % metal_complex.target_id)
                if metal_complex.location.startswith('protein'):
                    cmd.select('tmp_t', 'tmp_t & %s' % self.protname)
                    cmd.select('Metal-P', 'Metal-P or (id %s & %s)' % (metal_complex.target_id, self.protname))
                if metal_complex.location == 'ligand':
                    cmd.select('tmp_t', 'tmp_t & %s' % self.ligname)
                    cmd.select('Metal-L', 'Metal-L or (id %s & %s)' % (metal_complex.target_id, self.ligname))
                cmd.distance('MetalComplexes', 'tmp_m', 'tmp_t')
                cmd.delete('tmp_m or tmp_t')
        if self.object_exists('MetalComplexes'):
            cmd.set('dash_color', 'violetpurple', 'MetalComplexes')
            cmd.set('dash_gap', 0.5, 'MetalComplexes')
            # Show water molecules for metal complexes
            cmd.show('spheres', 'Metal-W')
            cmd.color('lightblue', 'Metal-W')



    def selections_cleanup(self):
        """Cleans up non-used selections"""

        if not len(self.plcomplex.unpaired_hba_idx) == 0:
            self.select_by_ids('Unpaired-HBA', self.plcomplex.unpaired_hba_idx, selection_exists=True)
        if not len(self.plcomplex.unpaired_hbd_idx) == 0:
            self.select_by_ids('Unpaired-HBD', self.plcomplex.unpaired_hbd_idx, selection_exists=True)
        if not len(self.plcomplex.unpaired_hal_idx) == 0:
            self.select_by_ids('Unpaired-HAL', self.plcomplex.unpaired_hal_idx, selection_exists=True)

        selections = cmd.get_names("selections")
        for selection in selections:
            try:
                empty = len(cmd.get_model(selection).atom) == 0
            except:
                empty = True
            if empty:
                cmd.delete(selection)
        cmd.deselect()
        cmd.delete('tmp*')
        cmd.delete('ps-*')

    def selections_group(self):
        """Group all selections"""
        cmd.group('Structures', '%s %s %sCartoon' % (self.protname, self.ligname, self.protname))
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

    def additional_cleanup(self):
        """Cleanup of various representations"""

        cmd.remove('not alt ""+A')  # Remove alternate conformations
        cmd.hide('labels', 'Interactions')  # Hide labels of lines
        cmd.disable('%sCartoon' % self.protname)
        cmd.hide('everything', 'hydrogens')

    def zoom_to_ligand(self):
        """Zoom in too ligand and its interactions."""
        cmd.center(self.ligname)
        cmd.orient(self.ligname)
        cmd.turn('x', 110)  # If the ligand is aligned with the longest axis, aromatic rings are hidden
        if 'AllBSRes' in cmd.get_names("selections"):
            cmd.zoom('%s or AllBSRes' % self.ligname, 3)
        else:
            if self.object_exists(self.ligname):
                cmd.zoom(self.ligname, 3)
        cmd.origin(self.ligname)

    def save_session(self, outfolder, override=None):
        """Saves a PyMOL session file."""
        filename = '%s_%s' % (self.protname.upper(), "_".join([self.hetid, self.plcomplex.chain, self.plcomplex.position]))
        if override is not None:
            filename = override
        cmd.save("/".join([outfolder, "%s.pse" % filename]))

    def png_workaround(self, filepath, width=1200, height=800):
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
        if subprocess.call("type convert", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0:
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

    def save_picture(self, outfolder, filename):
        """Saves a picture"""
        self.set_fancy_ray()
        self.png_workaround("/".join([outfolder, filename]))

    def set_fancy_ray(self):
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
        cmd.set('ray_opaque_background', 0) # Transparent background

    def adapt_for_peptides(self):
        """Adapt visualization for peptide ligands and interchain contacts"""
        cmd.hide('sticks', self.ligname)
        cmd.set('cartoon_color', 'lightorange', self.ligname)
        cmd.show('cartoon', self.ligname)
        cmd.show('sticks', "byres *-L")
        cmd.util.cnc(self.ligname)
        cmd.remove('%sCartoon and chain %s' % (self.protname, self.plcomplex.chain))
        cmd.set('cartoon_side_chain_helper', 0)

    def adapt_for_intra(self):
        """Adapt visualization for intra-protein interactions"""





    def refinements(self):
        """Refinements for the visualization"""

        # Show sticks for all residues interacing with the ligand
        cmd.select('AllBSRes', 'byres (Hydrophobic-P or HBondDonor-P or HBondAccept-P or PosCharge-P or NegCharge-P or '
                               'StackRings-P or PiCatRing-P or HalogenAcc or Metal-P)')
        cmd.show('sticks', 'AllBSRes')
        # Show spheres for the ring centroids
        cmd.hide('everything', 'centroids*')
        cmd.show('nb_spheres', 'centroids*')
        # Show spheres for centers of charge
        if self.object_exists('Chargecenter-P') or self.object_exists('Chargecenter-L'):
            cmd.hide('nonbonded', 'chargecenter*')
            cmd.show('spheres', 'chargecenter*')
            cmd.set('sphere_scale', 0.4, 'chargecenter*')
            cmd.color('yellow', 'chargecenter*')

        cmd.set('valence', 1)  # Show bond valency (e.g. double bonds)
        # Optional cartoon representation of the protein
        cmd.copy('%sCartoon' % self.protname, self.protname)
        cmd.show('cartoon', '%sCartoon' % self.protname)
        cmd.show('sticks', '%sCartoon' % self.protname)
        cmd.set('stick_transparency', 1, '%sCartoon' % self.protname)


        # Resize water molecules. Sometimes they are not heteroatoms HOH, but part of the protein
        cmd.set('sphere_scale', 0.2, 'resn HOH or Water')  # Needs to be done here because of the copy made
        cmd.set('sphere_transparency', 0.4, '!(resn HOH or Water)')

        if 'Centroids*' in cmd.get_names("selections"):
            cmd.color('grey80', 'Centroids*')
        cmd.hide('spheres', '%sCartoon' % self.protname)
        cmd.hide('cartoon', '%sCartoon and resn DA+DG+DC+DU+DT+A+G+C+U+T' % self.protname)  # Hide DNA/RNA Cartoon
        if self.ligname == 'SF4':  # Special case for iron-sulfur clusters, can't be visualized with sticks
            cmd.show('spheres', '%s' % self.ligname)

        cmd.hide('everything', 'resn HOH &!Water')  # Hide all non-interacting water molecules
        cmd.hide('sticks', '%s and !%s and !AllBSRes' % (self.protname, self.ligname))  # Hide all non-interacting residues

        if self.ligandtype in ['PEPTIDE', 'INTRA']:
            self.adapt_for_peptides()

        if self.ligandtype == 'INTRA':
            self.adapt_for_intra()
