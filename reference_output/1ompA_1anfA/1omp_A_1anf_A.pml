# DynDom Hierarchical Visualization
# Domain-colored structure with hierarchical screw axis arrows
reinitialize
load 1omp_A_1anf_A.pdb
bg_color grey
color grey

# === HIERARCHICAL DOMAIN STRUCTURE COLORING ===
select domain_0, resi 110-252
select domain_0, domain_0 + resi 314-321
select domain_0, domain_0 + resi 331-366
color blue, domain_0  # Domain 0 (GLOBAL REFERENCE)

select domain_1, resi 1-104
select domain_1, domain_1 + resi 259-306
color yellow, domain_1  # Domain 1 (MOVING)

# Color bending residues
select bending_residues_1, resi 105-109
color green, bending_residues_1
select bending_residues_2, resi 253-258
color green, bending_residues_2
select bending_residues_3, resi 307-313
color green, bending_residues_3
select bending_residues_4, resi 322-330
color green, bending_residues_4

set dash_gap, 0
set dash_radius, 0.2

# === HIERARCHICAL SCREW AXIS ARROWS ===
load output
load 1omp_A_1anf_A_arrows.pdb

# Basic protein display
hide everything, output
show cartoon, output
color gray80, output

# Hide arrow atoms initially
hide everything, 1omp_A_1anf_A_arrows

# Arrow 1: Domain 1 (moving) relative to Domain 0 (reference)
# Shaft color: blue (reference domain), Head color: yellow (moving domain)
# Rotation: 36.0°

# Select shaft and head atoms by chain and residue
select shaft_1, chain A and resn SHF and resi 100
select head_1, chain A and resn ARH and resi 120

# Display shaft as thick licorice stick (REFERENCE domain color: blue)
show sticks, shaft_1
color blue, shaft_1
set stick_radius, 0.3, shaft_1

# Display arrow head as clean cone (MOVING domain color: yellow)
show sticks, head_1
color yellow, head_1
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
print 'Global reference domain: 0 (blue)'
print 'Analysis pair 1: Domain 1 (yellow) relative to Domain 0 (blue)'
print '  Rotation: 36.0°'
