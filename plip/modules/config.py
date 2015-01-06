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

# Configuration file for Protein-Ligand Interaction Profiler (PLIP)
# Set thresholds for detection of interactions

# Thresholds for detection (global variables)
BS_DIST = 6.0  # Determines maximum distance to include binding site residues
AROMATIC_PLANARITY = 7.5  # Determines allowed deviation from planarity in aromatic rings

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