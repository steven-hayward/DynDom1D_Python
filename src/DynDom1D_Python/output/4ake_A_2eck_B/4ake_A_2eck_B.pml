# DynDom Hierarchical Visualization
# Domain-colored structure with hierarchical screw axis arrows
reinitialize
load 4ake_A_2eck_B.pdb
bg_color white
color grey

# === HIERARCHICAL DOMAIN STRUCTURE COLORING ===
select domain_0, resi 116-152
color red, domain_0  # Domain 0 (MOVING)

select domain_1, resi 29-58
color yellow, domain_1  # Domain 1 (MOVING)

select domain_2, resi 1-25
select domain_2, domain_2 + resi 63-111
select domain_2, domain_2 + resi 169-210
color blue, domain_2  # Domain 2 (GLOBAL REFERENCE)

# Color bending residues
select bending_residues_1, resi 112-115
color green, bending_residues_1
select bending_residues_2, resi 153-168
color green, bending_residues_2
select bending_residues_3, resi 26-28
color green, bending_residues_3
select bending_residues_4, resi 59-62
color green, bending_residues_4

set dash_gap, 0
set dash_radius, 0.2

# === HIERARCHICAL SCREW AXIS ARROWS ===
load output
load 4ake_A_2eck_B_arrows.pdb

# Basic protein display
hide everything, output
show cartoon, output
color gray80, output

# Hide arrow atoms initially
hide everything, 4ake_A_2eck_B_arrows

# Arrow 1: Domain 0 (moving) relative to Domain 2 (reference)
# Shaft color: blue (reference domain), Head color: red (moving domain)
# Rotation: 53.0°

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

# Arrow 2: Domain 1 (moving) relative to Domain 2 (reference)
# Shaft color: blue (reference domain), Head color: yellow (moving domain)
# Rotation: 46.1°

# Select shaft and head atoms by chain and residue
select shaft_2, chain B and resn SHF and resi 150
select head_2, chain B and resn ARH and resi 170

# Display shaft as thick licorice stick (REFERENCE domain color: blue)
show sticks, shaft_2
color blue, shaft_2
set stick_radius, 0.3, shaft_2

# Display arrow head as clean cone (MOVING domain color: yellow)
show sticks, head_2
color yellow, head_2
set stick_radius, 0.25, head_2

# Connect atoms ONLY within each section
bond shaft_2, shaft_2
bond head_2, head_2

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
print 'Global reference domain: 2 (blue)'
print 'Analysis pair 1: Domain 0 (red) relative to Domain 2 (blue)'
print '  Rotation: 53.0°'
print 'Analysis pair 2: Domain 1 (yellow) relative to Domain 2 (blue)'
print '  Rotation: 46.1°'
