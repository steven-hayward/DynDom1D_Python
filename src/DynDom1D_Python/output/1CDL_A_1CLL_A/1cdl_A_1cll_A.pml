# DynDom Hierarchical Visualization
# Domain-colored structure with hierarchical screw axis arrows
reinitialize
load 1cdl_A_1cll_A.pdb
bg_color grey
color grey

# === HIERARCHICAL DOMAIN STRUCTURE COLORING ===
select domain_0, resi 77-142
color red, domain_0  # Domain 0 (MOVING)

select domain_1, resi 5-67
color blue, domain_1  # Domain 1 (GLOBAL REFERENCE)

# Color bending residues
select bending_residues_1, resi 68-76
color green, bending_residues_1

set dash_gap, 0
set dash_radius, 0.2

# === HIERARCHICAL SCREW AXIS ARROWS ===
load output
load 1cdl_A_1cll_A_arrows.pdb

# Basic protein display
hide everything, output
show cartoon, output
color gray80, output

# Hide arrow atoms initially
hide everything, 1cdl_A_1cll_A_arrows

# Arrow 1: Domain 0 (moving) relative to Domain 1 (reference)
# Shaft color: blue (reference domain), Head color: red (moving domain)
# Rotation: 154.7°

# Select shaft and head atoms by chain and residue
select shaft_1, chain A and resn SHF and resi 100
select head_1, chain A and resn ARH and resi 120

# Display shaft as thick licorice stick (REFERENCE domain color: blue)
show sticks, shaft_1
color blue, shaft_1
set stick_radius, 0.3, shaft_1

# Display arrow head as clean cone (MOVING domain color: red)
show sticks, head_1
color red, head_1
set stick_radius, 0.25, head_1

# Connect atoms ONLY within each section
bond shaft_1, shaft_1
bond head_1, head_1

# Disable automatic bonding between different chains
set auto_bond, 0

# Make arrows more prominent
set stick_transparency, 0.0
set stick_quality, 15
set sphere_quality, 3
set surface_quality, 2

# Final settings
set depth_cue, 0
set ray_shadows, 1
set ray_shadow_decay_factor, 0.1

# Better lighting for 3D arrow heads
set ambient, 0.2
set direct, 0.8
set reflect, 0.5
set shininess, 10

# Clean up selections
delete shaft_*
delete head_*

# === FINAL SETTINGS ===
set stick_transparency, 0.0
set stick_quality, 15
zoom all
orient

# Cleanup selections
delete bending_residues
delete arrow_*

print 'DynDom hierarchical visualization loaded!'
print 'Global reference domain: 1 (blue)'
print 'Analysis pair 1: Domain 0 (red) relative to Domain 1 (blue)'
print '  Rotation: 154.7°'
