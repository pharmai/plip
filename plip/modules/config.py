"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
config.py - Store thresholds used by PLIP.
Copyright 2014 Sebastian Salentin

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

# Configuration file for Protein-Ligand Interaction Profiler (PLIP)
# Set thresholds for detection of interactions

# Thresholds for detection (global variables)
BS_DIST = 8.5  # Determines maximum distance to include binding site residues
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