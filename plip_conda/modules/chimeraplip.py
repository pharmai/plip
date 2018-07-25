"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
chimeraplip.py - Visualization class for Chimera.
"""

class ChimeraVisualizer():
    """Provides visualization for Chimera."""
    def __init__(self, plcomplex, chimera_module, tid):
        self.chimera = chimera_module
        self.tid = tid
        self.uid = plcomplex.uid
        self.plipname = 'PLIP-%i' % self.tid
        self.hetid, self.chain, self.pos = self.uid.split(':')
        self.pos = int(self.pos)
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

            for md in self.models.list():
                if md.name == self.plipname:
                    self.model = md

            self.atoms = self.atom_by_serialnumber()


    def set_initial_representations(self):
        """Set the initial representations"""
        self.update_model_dict()
        self.rc("background solid white")
        self.rc("setattr g display 0")  # Hide all pseudobonds
        self.rc("~display #%i & :/isHet & ~:%s" % (self.model_dict[self.plipname],self.hetid))

    def update_model_dict(self):
        """Updates the model dictionary"""
        dct = {}
        models = self.chimera.openModels
        for md in models.list():
            dct[md.name] = md.id
        self.model_dict = dct

    def atom_by_serialnumber(self):
        """Provides a dictionary mapping serial numbers to their atom objects."""
        atm_by_snum = {}
        for atom in self.model.atoms:
            atm_by_snum[atom.serialNumber] =  atom
        return atm_by_snum

    def show_hydrophobic(self):
        """Visualizes hydrophobic contacts."""
        grp = self.getPseudoBondGroup("Hydrophobic Interactions-%i" % self.tid, associateWith=[self.model])
        grp.lineType = self.chimera.Dash
        grp.lineWidth = 3
        grp.color = self.colorbyname('gray')
        for i in self.plcomplex.hydrophobic_contacts.pairs_ids:
            b = grp.newPseudoBond(self.atoms[i[0]], self.atoms[i[1]])
            self.bs_res_ids.append(i[0])

    def show_hbonds(self):
        """Visualizes hydrogen bonds."""
        grp = self.getPseudoBondGroup("Hydrogen Bonds-%i" % self.tid, associateWith=[self.model])
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
        grp = self.getPseudoBondGroup("HalogenBonds-%i" % self.tid, associateWith=[self.model])
        grp.lineWidth = 3
        for i in self.plcomplex.halogen_bonds:
            b = grp.newPseudoBond(self.atoms[i[0]], self.atoms[i[1]])
            b.color = self.colorbyname('turquoise')

            self.bs_res_ids.append(i.acc_id)

    def show_stacking(self):
        """Visualizes pi-stacking interactions."""
        grp = self.getPseudoBondGroup("pi-Stacking-%i" % self.tid, associateWith=[self.model])
        grp.lineWidth = 3
        grp.lineType = self.chimera.Dash
        for i, stack in enumerate(self.plcomplex.pistacking):

            m = self.model
            r = m.newResidue("pseudoatoms", " ", 1, " ")
            centroid_prot = m.newAtom("CENTROID", self.chimera.Element("CENTROID"))
            x, y, z = stack.proteinring_center
            centroid_prot.setCoord(self.chimera.Coord(x, y, z))
            r.addAtom(centroid_prot)

            centroid_lig = m.newAtom("CENTROID", self.chimera.Element("CENTROID"))
            x, y, z = stack.ligandring_center
            centroid_lig.setCoord(self.chimera.Coord(x, y, z))
            r.addAtom(centroid_lig)


            b = grp.newPseudoBond(centroid_lig, centroid_prot)
            b.color = self.colorbyname('forest green')

            self.bs_res_ids += stack.proteinring_atoms

    def show_cationpi(self):
        """Visualizes cation-pi interactions"""
        grp = self.getPseudoBondGroup("Cation-Pi-%i" % self.tid, associateWith=[self.model])
        grp.lineWidth = 3
        grp.lineType = self.chimera.Dash
        for i, cat in enumerate(self.plcomplex.pication):

            m = self.model
            r = m.newResidue("pseudoatoms", " ", 1, " ")
            chargecenter = m.newAtom("CHARGE", self.chimera.Element("CHARGE"))
            x, y, z = cat.charge_center
            chargecenter.setCoord(self.chimera.Coord(x, y, z))
            r.addAtom(chargecenter)

            centroid = m.newAtom("CENTROID", self.chimera.Element("CENTROID"))
            x, y, z = cat.ring_center
            centroid.setCoord(self.chimera.Coord(x, y, z))
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
        grp = self.getPseudoBondGroup("Salt Bridges-%i" % self.tid, associateWith=[self.model])
        grp.lineWidth = 3
        grp.lineType = self.chimera.Dash
        for i, sbridge in enumerate(self.plcomplex.saltbridges):

            m = self.model
            r = m.newResidue("pseudoatoms", " ", 1, " ")
            chargecenter1 = m.newAtom("CHARGE", self.chimera.Element("CHARGE"))
            x, y, z = sbridge.positive_center
            chargecenter1.setCoord(self.chimera.Coord(x, y, z))
            r.addAtom(chargecenter1)

            chargecenter2 = m.newAtom("CHARGE", self.chimera.Element("CHARGE"))
            x, y, z = sbridge.negative_center
            chargecenter2.setCoord(self.chimera.Coord(x, y, z))
            r.addAtom(chargecenter2)

            b = grp.newPseudoBond(chargecenter1, chargecenter2)
            b.color = self.colorbyname('yellow')

            if sbridge.protispos:
                self.bs_res_ids += sbridge.positive_atoms
            else:
                self.bs_res_ids += sbridge.negative_atoms

    def show_wbridges(self):
        """Visualizes water bridges"""
        grp = self.getPseudoBondGroup("Water Bridges-%i" % self.tid, associateWith=[self.model])
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
        grp = self.getPseudoBondGroup("Metal Coordination-%i" % self.tid, associateWith=[self.model])
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

        if not len(self.water_ids) == 0:
            # Hide all non-interacting water molecules
            water_selection = []
            for wid in self.water_ids:
                water_selection.append('serialNumber=%i' % wid)
            self.rc("~display :HOH")
            self.rc("display :@/%s" % " or ".join(water_selection))

        # Show all interacting binding site residues
        self.rc("~display #%i & ~:/isHet" % self.model_dict[self.plipname])
        self.rc("display :%s" % ",".join([str(self.atoms[bsid].residue.id) for bsid in self.bs_res_ids]))
        self.rc("color lightblue :HOH")



    def zoom_to_ligand(self):
        """Centers the view on the ligand and its binding site residues."""
        self.rc("center #%i & :%s" % (self.model_dict[self.plipname], self.hetid))

    def refinements(self):
        """Details for the visualization."""
        self.rc("setattr a color gray @CENTROID")
        self.rc("setattr a radius 0.3 @CENTROID")
        self.rc("represent sphere @CENTROID")
        self.rc("setattr a color orange @CHARGE")
        self.rc("setattr a radius 0.4 @CHARGE")
        self.rc("represent sphere @CHARGE")
        self.rc("display :pseudoatoms")
