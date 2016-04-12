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
from collections import namedtuple

# Python Standard Library
import json

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


class PyMOLVisualizer:

    def __init__(self, plcomplex):
        if not plcomplex is None:
            self.plcomplex = plcomplex
            self.protname = plcomplex.pdbid  # Name of protein with binding site
            self.ligname = plcomplex.hetid  # Name of ligand
            self.metal_ids = plcomplex.metal_ids

    def set_initial_representations(self):
        """General settings for PyMOL"""
        standard_settings()
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

    def show_hydrophobic(self):
        """Visualizes hydrophobic contacts."""
        hydroph = self.plcomplex.hydrophobic_contacts
        if not len(hydroph.bs_ids) == 0:
            select_by_ids('Hydrophobic-P', hydroph.bs_ids, restrict=self.protname)
            select_by_ids('Hydrophobic-L', hydroph.lig_ids, restrict=self.ligname)
            for i in hydroph.pairs_ids:
                cmd.select('tmp_bs', 'id %i & %s' % (i[0], self.protname))
                cmd.select('tmp_lig', 'id %i & %s' % (i[1], self.ligname))
                cmd.distance('Hydrophobic', 'tmp_bs', 'tmp_lig')
            if object_exists('Hydrophobic'):
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
                select_by_ids(group[0], group[1], restrict=self.protname)
        for group in [['HBondDonor-L', hbonds.lig_don_id],
        ['HBondAccept-L', hbonds.lig_acc_id]]:
            if not len(group[1]) == 0:
                select_by_ids(group[0], group[1], restrict=self.ligname)
        for i in hbonds.ldon_id:
            cmd.select('tmp_bs', 'id %i & %s' % (i[0], self.protname))
            cmd.select('tmp_lig', 'id %i & %s' % (i[1], self.ligname))
            cmd.distance('HBonds', 'tmp_bs', 'tmp_lig')
        for i in hbonds.pdon_id:
            cmd.select('tmp_bs', 'id %i & %s' % (i[1], self.protname))
            cmd.select('tmp_lig', 'id %i & %s' % (i[0], self.ligname))
            cmd.distance('HBonds', 'tmp_bs', 'tmp_lig')
        if object_exists('HBonds'):
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
            select_by_ids('HalogenAccept', all_acc_o, restrict=self.protname)
            select_by_ids('HalogenDonor', all_don_x, restrict=self.ligname)
        if object_exists('HalogenBonds'):
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
        if object_exists('PiStackingP'):
            cmd.set('dash_color', 'green', 'PiStackingP')
            cmd.set('dash_gap', 0.3, 'PiStackingP')
            cmd.set('dash_length', 0.6, 'PiStackingP')
        if object_exists('PiStackingT'):
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
        if object_exists('PiCation'):
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

        if object_exists('Saltbridges'):
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
        if object_exists('WaterBridges'):
            cmd.set('dash_color', 'lightblue', 'WaterBridges')
        cmd.delete('tmp_water or tmp_acc or tmp_don')
        cmd.color('lightblue', 'Water')
        cmd.show('spheres', 'Water')

    def show_metal(self):
        """Visualize metal coordination."""
        metal_complexes = self.plcomplex.metal_complexes
        if not len(metal_complexes) == 0:
            select_by_ids('Metal-M', self.metal_ids)
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
        if object_exists('MetalComplexes'):
            cmd.set('dash_color', 'violetpurple', 'MetalComplexes')
            cmd.set('dash_gap', 0.5, 'MetalComplexes')
            # Show water molecules for metal complexes
            cmd.show('spheres', 'Metal-W')
            cmd.color('lightblue', 'Metal-W')



    def selections_cleanup(self):
        """Cleans up non-used selections"""

        if not len(self.plcomplex.unpaired_hba_idx) == 0:
            select_by_ids('Unpaired-HBA', self.plcomplex.unpaired_hba_idx, selection_exists=True)
        if not len(self.plcomplex.unpaired_hbd_idx) == 0:
            select_by_ids('Unpaired-HBD', self.plcomplex.unpaired_hbd_idx, selection_exists=True)
        if not len(self.plcomplex.unpaired_hal_idx) == 0:
            select_by_ids('Unpaired-HAL', self.plcomplex.unpaired_hal_idx, selection_exists=True)

        selections = cmd.get_names("selections")
        for selection in selections:
            if len(cmd.get_model(selection).atom) == 0:
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
            if object_exists(self.ligname):
                cmd.zoom(self.ligname, 3)
        cmd.origin(self.ligname)

    def generate_output(self):
        """Saves PyMOL session files and images."""
        filename = '%s_%s' % (self.protname.upper(), "_".join([self.ligname, self.plcomplex.chain, self.plcomplex.position]))
        if config.PYMOL:
            cmd.save("/".join([config.OUTPATH, "%s.pse" % filename]))

        # Create output pictures (experimental)
        set_fancy_ray()
        if config.PICS:
            png_workaround("/".join([config.OUTPATH, filename]))

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
        if object_exists('Chargecenter-P') or object_exists('Chargecenter-L'):
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


class ChimeraVisualizer():
    """Provides visualization for Chimera."""
    def __init__(self, plcomplex, chimera_module):
        self.chimera = chimera_module
        Molecule = self.chimera.Molecule
        self.colorbyname = self.chimera.colorTable.getColorByName
        self.rc = self.chimera.runCommand
        self.getPseudoBondGroup = self.chimera.misc.getPseudoBondGroup
        if not plcomplex is None:
            self.plcomplex = plcomplex
            self.protname = plcomplex.pdbid  # Name of protein with binding site
            self.ligname = plcomplex.hetid  # Name of ligand
            self.metal_ids = plcomplex.metal_ids
            self.water_ids = []
            self.bs_res_ids = []

            self.models = self.chimera.openModels
            self.model = self.models.list()[0]  # Model of interest (ref protein)
            self.helper_mol = Molecule()  # Chimera Molecule
            self.models.add([self.helper_mol])

            self.atoms = self.atom_by_serialnumber()

    def set_initial_representations(self):
        self.rc("background solid white")

    def make_initial_selections(self):
        pass

    def atom_by_serialnumber(self):
        """Provides a dictionary mapping serial numbers to their atom objects."""
        atm_by_snum = {}
        for atom in self.model.atoms:
            atm_by_snum[atom.serialNumber] =  atom
        return atm_by_snum

    def show_hydrophobic(self):
        """Visualizes hydrophobic contacts."""
        grp = self.getPseudoBondGroup("Hydrophobic Interactions", associateWith=[self.model])
        grp.lineType = self.chimera.Dash
        grp.lineWidth = 3
        grp.color = self.colorbyname('gray')
        for i in self.plcomplex.hydrophobic_contacts.pairs_ids:
            b = grp.newPseudoBond(self.atoms[i[0]], self.atoms[i[1]])
            self.bs_res_ids.append(i[0])

    def show_hbonds(self):
        """Visualizes hydrogen bonds."""
        grp = self.getPseudoBondGroup("Hydrogen Bonds", associateWith=[self.model])
        grp.lineWidth = 3
        for i in self.plcomplex.hbonds.ldon_id:
            b = grp.newPseudoBond(self.atoms[i[0]], self.atoms[i[1]])
            b.color = self.colorbyname('blue')
            self.bs_res_ids.append(i[0])
        for i in self.plcomplex.hbonds.pdon_id:
            b = grp.newPseudoBond(self.atoms[i[0]], self.atoms[i[1]])
            b.color = self.colorbyname('blue')
            self.bs_res_ids.append(i[1])


    def show_halogen(self):
        """Visualizes halogen bonds."""
        grp = self.getPseudoBondGroup("HalogenBonds", associateWith=[self.model])
        grp.lineWidth = 3
        for i in self.plcomplex.halogen_bonds:
            b = grp.newPseudoBond(self.atoms[i[0]], self.atoms[i[1]])
            b.color = self.colorbyname('turquoise')

            self.bs_res_ids.append(i.acc_id)

    def show_stacking(self):
        """Visualizes pi-stacking interactions."""
        grp = self.getPseudoBondGroup("pi-Stacking", associateWith=[self.model])
        grp.lineWidth = 3
        grp.lineType = self.chimera.Dash
        for i, stack in enumerate(self.plcomplex.pistacking):
            #@TODO Only one new model for centroids, chargecenters etc.

            m = self.helper_mol
            r = m.newResidue("Centroids", " ", 1, " ")
            centroid_prot = m.newAtom("CENTROID", Element("CENTROID"))
            x, y, z = stack.proteinring_center
            centroid_prot.setCoord(Coord(x, y, z))
            r.addAtom(centroid_prot)

            centroid_lig = m.newAtom("CENTROID", Element("CENTROID"))
            x, y, z = stack.ligandring_center
            centroid_lig.setCoord(Coord(x, y, z))
            r.addAtom(centroid_lig)


            #models.add([m])
            b = grp.newPseudoBond(centroid_lig, centroid_prot)
            b.color = self.colorbyname('forest green')

            self.bs_res_ids += stack.proteinring_atoms

    def show_cationpi(self):
        """Visualizes cation-pi interactions"""
        grp = self.getPseudoBondGroup("Cation-Pi", associateWith=[self.model])
        grp.lineWidth = 3
        grp.lineType = self.chimera.Dash
        for i, cat in enumerate(self.plcomplex.pication):

            m = self.helper_mol
            r = m.newResidue("Centroids2", " ", 1, " ")
            chargecenter = m.newAtom("CHARGE", Element("CHARGE"))
            x, y, z = cat.charge_center
            chargecenter.setCoord(Coord(x, y, z))
            r.addAtom(chargecenter)

            centroid = m.newAtom("CENTROID", Element("CENTROID"))
            x, y, z = cat.ring_center
            centroid.setCoord(Coord(x, y, z))
            r.addAtom(centroid)

            b = grp.newPseudoBond(centroid, chargecenter)
            b.color = self.colorbyname('orange')

            if cat.protcharged:
                self.bs_res_ids += cat.charge_atoms
            else:
                self.bs_res_ids += cat.ring_atoms

    def show_sbridges(self):
        """Visualizes salt bridges."""
        # Salt Bridges
        grp = self.getPseudoBondGroup("Salt Bridges", associateWith=[self.model])
        grp.lineWidth = 3
        grp.lineType = self.chimera.Dash
        for i, sbridge in enumerate(self.plcomplex.saltbridges):

            m = self.helper_mol
            r = m.newResidue("Centroids3", " ", 1, " ")
            chargecenter1 = m.newAtom("CHARGE", Element("CHARGE"))
            x, y, z = sbridge.positive_center
            chargecenter1.setCoord(Coord(x, y, z))
            r.addAtom(chargecenter1)

            chargecenter2 = m.newAtom("CHARGE", Element("CHARGE"))
            x, y, z = sbridge.negative_center
            chargecenter2.setCoord(Coord(x, y, z))
            r.addAtom(chargecenter2)

            b = grp.newPseudoBond(chargecenter1, chargecenter2)
            b.color = self.colorbyname('yellow')

            if sbridge.protispos:
                self.bs_res_ids += sbridge.positive_atoms
            else:
                self.bs_res_ids += sbridge_negative_atoms

    def show_wbridges(self):
        """Visualizes water bridges"""
        grp = self.getPseudoBondGroup("Water Bridges", associateWith=[self.model])
        grp.lineWidth = 3
        for i, wbridge in enumerate(self.plcomplex.waterbridges):
            c = grp.newPseudoBond(self.atoms[wbridge.water_id], self.atoms[wbridge.acc_id])
            c.color = self.colorbyname('cornflower blue')
            self.water_ids.append(wbridge.water_id)
            b = grp.newPseudoBond(self.atoms[wbridge.don_id], self.atoms[wbridge.water_id])
            b.color = self.colorbyname('cornflower blue')
            self.water_ids.append(wbridge.water_id)
            if wbridge.protisdon:
                self.bs_res_ids.append(wbridge.don_id)
            else:
                self.bs_res_ids.append(wbridge.acc_id)

    def show_metal(self):
        """Visualizes metal coordination."""
        # Metal Coordination
        grp = self.getPseudoBondGroup("Metal Coordination", associateWith=[self.model])
        grp.lineWidth = 3
        for i, metal in enumerate(self.plcomplex.metal_complexes):
            c = grp.newPseudoBond(self.atoms[metal.metal_id], self.atoms[metal.target_id])
            c.color = self.colorbyname('magenta')

            if metal.location == 'water':
                self.water_ids.append(metal.target_id)

            if metal.location.startswith('protein'):
                self.bs_res_ids.append(metal.target_id)

    def cleanup(self):
        """Clean up the visualization."""

        # Hide all non-interacting water molecules
        water_selection = []
        for wid in self.water_ids:
            water_selection.append('serialNumber=%i' % wid)
        self.rc("~display :HOH & ~:@/%s" % " or ".join(water_selection))

        # Show all interacting binding site residues
        self.rc("display :%s" % ",".join([str(self.atoms[bsid].residue.id) for bsid in self.bs_res_ids]))


    def zoom_to_ligand(self):
        pass

    def refinements(self):
        """Details for the visualization."""
        self.rc("setattr a color gray @CENTROID")
        self.rc("setattr a radius 0.3 @CENTROID")
        self.rc("setattr a color orange @CHARGE")
        self.rc("setattr a radius 0.5 @CHARGE")



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
        select_by_ids(ligname, metal_ids, selection_exists=True)
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

    vis.generate_output()
