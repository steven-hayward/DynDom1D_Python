import traceback
import urllib.request
import numpy as np
import os
from itertools import groupby
from operator import itemgetter
from pathlib import Path

input_command_file_path = "data"
DOMAIN_COLORS = ["blue", "red", "yellow", "pink", "cyan", "purple", "orange", "brown", "black", "white", "magenta", "violet", "indigo", "turquoise", "coral"]

def read_command_file():
    temp_dict = {}
    try:
        fr = open(f"{input_command_file_path}/command.txt", "r")
        lines = fr.readlines()
        for line in lines:
            if not ("#" in line):
                line = line.replace("\n", "")
                line = line.replace(" ", "")
                tokens = line.split("=")
                param_name = tokens[0]
                param_val = tokens[1]
                if "chain" in param_name:
                    param_val = param_val.upper()
                elif "filename" in param_name:
                    param_val = param_val.lower()
                temp_dict[param_name] = param_val
        fr.close()
        check_pdb_exists(temp_dict['input_path'], temp_dict['filename1'])
        check_pdb_exists(temp_dict['input_path'], temp_dict['filename2'])
    except Exception as e:
        print(e)
    return temp_dict


def check_pdb_exists(input_path, file_name):
    """
    Uses HTTPS request to make sure that the pdb file is in the directory
    :param input_path: The directory to the file
    :param file_name:
    :return:
    """
    file_path = f"{input_path}/{file_name}.pdb"
    path = Path(file_path)
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{file_name}.pdb",
                file_path
            )
        except Exception as e:
            traceback.print_exc()
            print(e)


def read_param_file():
    temp_dict = {}
    try:
        fr = open(f"{input_command_file_path}/param.txt", "r")
        lines = fr.readlines()
        for line in lines:
            if not ("#" in line):
                line = line.replace("\n", "")
                line = line.replace(" ", "")
                tokens = line.split("=")
                param_name = tokens[0]
                param_val = tokens[1]
                if param_name == "window":
                    param_val = int(tokens[1])
                elif param_name == "domain":
                    param_val = int(tokens[1])
                elif param_name == "ratio":
                    param_val = float(tokens[1])
                elif param_name == "k_means_n_init":
                    param_val = int(tokens[1])
                elif param_name == "k_means_max_iter":
                    param_val = int(tokens[1])
                temp_dict[param_name] = param_val
        fr.close()
    except Exception as e:
        print(e)
    return temp_dict


def write_rotation_vec_to_pdb(output_path, protein_1, protein_2_name, chain_2, rotation_vectors):
    """
    Writes the rotation vectors of each residue of the slide window into a pdb file
    :param output_path:
    :param protein_1: The Protein object of protein 1
    :param protein_2_name:
    :param chain_2:
    :param rotation_vectors:
    :return:
    """
    try:
        dir_path_str = f"{output_path}/{protein_1.name}_{protein_1.chain_param}_{protein_2_name}_{chain_2}"
        dir_path = Path(dir_path_str)
        if not dir_path.exists():
            dir_path.mkdir(parents=True)
        fw = open(f"{dir_path_str}/{protein_1.name}_rot_vecs.pdb", "w")
        slide_window_indices = protein_1.slide_window_residues_indices
        protein_polymer = protein_1.get_polymer()
        for i in range(rotation_vectors.shape[0]):
            index = protein_1.utilised_residues_indices[i+slide_window_indices[0]]
            residue_name = protein_polymer[index].name
            residue_num = protein_polymer[index].seqid.num
            x = f"{rotation_vectors[i][0]:8.3f}"
            y = f"{rotation_vectors[i][1]:8.3f}"
            z = f"{rotation_vectors[i][2]:8.3f}"
            row = f"ATOM  {i+1:5d}  CA  {residue_name:>3s} A{residue_num:4d}   {x}{y}{z}\n"
            fw.write(row)
        fw.close()
    except Exception as e:
        print(e)
        return False
    return True


def write_final_output_pdb(output_path, protein_1, fitted_protein_2, fitted_protein_2_name, fitted_protein_2_chain, fitted_protein_2_res_ind):
    """
    :param output_path:
    :param protein_1: The Protein object of protein 1
    :param fitted_protein_2: The Chain object of protein 2 fitted to protein 1
    :param fitted_protein_2_name: The name of the protein_2
    :param fitted_protein_2_chain: The chain of protein 2 fitted to protein 1
    :param fitted_protein_2_res_ind: The indices of residues used in the clustering
    :return:
    """
    try:
        folder_name = f"{protein_1.name}_{protein_1.chain_param}_{fitted_protein_2_name}_{fitted_protein_2_chain}"
        dir_path_str = f"{output_path}/{folder_name}"
        dir_path = Path(dir_path_str)
        if not dir_path.exists():
            dir_path.mkdir(parents=True)
        fw = open(f"{dir_path_str}/{folder_name}.pdb", "w")
        
        # Define the protein data for each model
        proteins_data = [
            {
                'model_num': 1,
                'residues': protein_1.get_polymer(),
                'res_indices': protein_1.utilised_residues_indices,
                'chain': protein_1.chain_param
            },
            {
                'model_num': 2,
                'residues': fitted_protein_2.get_polymer(),
                'res_indices': fitted_protein_2_res_ind,
                'chain': fitted_protein_2_chain
            }
        ]
        
        # Write each model using the same logic
        for protein_data in proteins_data:
            fw.write(f"MODEL{protein_data['model_num']:>9}\n")
            atom_count = 1
            
            for i in protein_data['res_indices']:
                r = protein_data['residues'][i]
                res_name = r.name.rjust(3, " ")
                res_num = str(r.seqid.num).rjust(4, " ")
                
                for a in r:
                    atom_num = f"{atom_count:5d}"
                    atom_name = f"{a.name:<4s}"
                    res_name = f"{r.name:>3s}"
                    res_num = f"{r.seqid.num:4d}"
                    x = f"{a.pos.x:8.3f}"
                    y = f"{a.pos.y:8.3f}"
                    z = f"{a.pos.z:8.3f}"
                    row = f"ATOM  {atom_num} {atom_name} {res_name} {protein_data['chain']}{res_num}    {x}{y}{z}\n"
                    fw.write(row)
                    atom_count += 1
            
            fw.write("ENDMDL\n")
        
        fw.close()
        
    except Exception as e:
        traceback.print_exc()
        print(e)

def write_arrows_pdb(output_path, base_name, analysis_pairs, domains, protein_1, pair_specific_results=None):
    """
    Generate arrows PDB file using hierarchical analysis pairs
    FIXED: Use pair-specific screw axis results instead of domain objects
    """
    arrow_pdb_path = f"{output_path}/{base_name}_arrows.pdb"
    
    try:
        if not analysis_pairs:
            print("No analysis pairs found - no arrows to generate")
            return None
            
        # Read protein coordinates for proper arrow scaling
        protein_coordinates = _read_protein_coordinates(f"{output_path}/{base_name}.pdb")
        
        with open(arrow_pdb_path, 'w') as f:
            f.write("REMARK DynDom Hierarchical Arrow Visualization\n")
            f.write("REMARK Shaft atoms = SHF residues (reference domain color)\n")
            f.write("REMARK Head atoms = ARH residues (moving domain color)\n")
            f.write("REMARK Each arrow uses separate chain ID (A, B, C...)\n")
            f.write("REMARK Generated from hierarchical analysis pairs\n")
            
            atom_id = 1000
            arrow_index = 0
            
            for moving_domain_id, reference_domain_id in analysis_pairs:
                moving_domain = domains[moving_domain_id]
                reference_domain = domains[reference_domain_id]
                
                # CRITICAL FIX: Use pair-specific results instead of domain object
                if pair_specific_results:
                    pair_key = (moving_domain_id, reference_domain_id)
                    if pair_key in pair_specific_results:
                        pair_results = pair_specific_results[pair_key]
                        rot_angle = pair_results['rot_angle']
                        screw_axis = pair_results['screw_axis']
                        point_on_axis = pair_results['point_on_axis']
                    else:
                        print(f"Warning: No pair-specific results for ({moving_domain_id}, {reference_domain_id})")
                        rot_angle = moving_domain.rot_angle
                        screw_axis = moving_domain.screw_axis
                        point_on_axis = moving_domain.point_on_axis
                else:
                    # Fallback to domain object
                    rot_angle = moving_domain.rot_angle
                    screw_axis = moving_domain.screw_axis
                    point_on_axis = moving_domain.point_on_axis
                
                f.write(f"REMARK Arrow {arrow_index+1}: Domain {moving_domain_id} moving relative to Domain {reference_domain_id}\n")
                f.write(f"REMARK   Chain {chr(ord('A') + arrow_index)}: Rotation angle: {rot_angle:.1f} degrees\n")
                
                # Generate arrow atoms using pair-specific data
                arrow_atoms = _create_arrow_atoms_hierarchical_with_data(
                    screw_axis, point_on_axis, rot_angle, arrow_index, atom_id, protein_coordinates
                )
                f.write('\n'.join(arrow_atoms))
                f.write('\n')
                
                atom_id += 200  # Large gap between arrows
                arrow_index += 1
                
            f.write("END\n")
            
        print(f"Arrow PDB created: {arrow_pdb_path}")
        return arrow_pdb_path
        
    except Exception as e:
        print(f"Error creating arrow PDB: {e}")
        traceback.print_exc()
        return None
    
    
def _create_arrow_atoms_hierarchical_with_data(screw_axis, point_on_axis, rot_angle, arrow_index, start_atom_id, protein_coordinates):
    """
    Create PDB atom lines for an arrow using provided screw axis data
    """
    # Normalize the screw axis
    screw_axis = screw_axis / np.linalg.norm(screw_axis)
    
    # Calculate appropriate arrow length
    base_shaft_length = _calculate_arrow_length(screw_axis, point_on_axis, protein_coordinates)
    head_length = 3.0
    
    # Extend shaft backwards by 10%
    shaft_extension = base_shaft_length * 0.1
    total_shaft_length = base_shaft_length + shaft_extension
    shaft_spacing = 1.0
    
    pdb_lines = []
    atom_id = start_atom_id
    
    # Use different chain/residue IDs for each arrow
    shaft_res_id = 100 + arrow_index * 50
    head_res_id = shaft_res_id + 20
    chain_id = chr(ord('A') + arrow_index)
    
    # === CREATE SHAFT ===
    n_shaft_atoms = int(total_shaft_length / shaft_spacing)
    for i in range(n_shaft_atoms):
        t = -(base_shaft_length/2 + shaft_extension) + i * shaft_spacing
        pos = point_on_axis + t * screw_axis
        pdb_line = f"ATOM  {atom_id:5d}  CA  SHF {chain_id}{shaft_res_id:4d}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00 30.00           C"
        pdb_lines.append(pdb_line)
        atom_id += 1
        
    pdb_lines.append("TER")
    
    # === CREATE ARROW HEAD ===
    # Create perpendicular vectors for cone base
    if abs(screw_axis[0]) < 0.9:
        temp_vec = np.array([1.0, 0.0, 0.0])
    else:
        temp_vec = np.array([0.0, 1.0, 0.0])
    
    perp1 = np.cross(screw_axis, temp_vec)
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(screw_axis, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)
    
    # Arrow head overlaps with shaft end
    head_overlap = 2.0
    head_start_pos = point_on_axis + (base_shaft_length/2 - head_overlap) * screw_axis
    
    # Create cone structure
    n_layers = 4
    n_points_per_layer = 6
    
    for layer in range(n_layers + 1):
        layer_progress = layer / n_layers
        layer_pos = head_start_pos + layer_progress * head_length * screw_axis
        
        max_radius = 1.2
        layer_radius = max_radius * (1.0 - layer_progress) ** 1.5
        
        if layer == n_layers:
            # Tip of arrow
            pdb_line = f"ATOM  {atom_id:5d}  CA  ARH {chain_id}{head_res_id:4d}    {layer_pos[0]:8.3f}{layer_pos[1]:8.3f}{layer_pos[2]:8.3f}  1.00 50.00           C"
            pdb_lines.append(pdb_line)
            atom_id += 1
        else:
            # Circular layer
            for point in range(n_points_per_layer):
                angle = 2 * np.pi * point / n_points_per_layer
                
                circle_pos = (layer_pos + 
                             layer_radius * np.cos(angle) * perp1 + 
                             layer_radius * np.sin(angle) * perp2)
                
                pdb_line = f"ATOM  {atom_id:5d}  CA  ARH {chain_id}{head_res_id:4d}    {circle_pos[0]:8.3f}{circle_pos[1]:8.3f}{circle_pos[2]:8.3f}  1.00 40.00           C"
                pdb_lines.append(pdb_line)
                atom_id += 1
    
    pdb_lines.append("TER")
    return pdb_lines


def _read_protein_coordinates(pdb_file):
    """Read protein coordinates for arrow scaling"""
    coordinates = []
    try:
        if os.path.exists(pdb_file):
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') and line[12:16].strip() in ['CA', 'N', 'C']:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coordinates.append([x, y, z])
    except Exception as e:
        print(f"Error reading protein coordinates: {e}")
    
    return np.array(coordinates) if coordinates else np.array([])


def _calculate_arrow_length(screw_axis, point_on_axis, protein_coordinates):
    """Calculate appropriate arrow length like axleng.f does"""
    if len(protein_coordinates) == 0:
        return 40.0  # Default length
        
    # Normalize screw axis
    unit_axis = screw_axis / np.linalg.norm(screw_axis)
    
    # Calculate projections of all protein atoms onto the screw axis
    relative_coords = protein_coordinates - point_on_axis
    projections = np.dot(relative_coords, unit_axis)
    
    # Find min and max projections
    tmin = np.min(projections)
    tmax = np.max(projections)
    
    # Add 10% padding on each end
    padding = (tmax - tmin) / 10.0
    tmin_padded = tmin - padding
    tmax_padded = tmax + padding
    
    # Total length along the axis
    total_length = max(tmax_padded - tmin_padded, 30.0)  # Minimum 30Å
    
    return total_length

def _create_arrow_atoms(domain, arrow_index, start_atom_id, protein_coordinates):
    """Create PDB atom lines for an arrow using domain data directly"""
    screw_axis = domain.screw_axis
    point_on_axis = domain.point_on_axis
    
    # Normalize the screw axis
    screw_axis = screw_axis / np.linalg.norm(screw_axis)
    
    # Calculate appropriate arrow length
    base_shaft_length = _calculate_arrow_length(screw_axis, point_on_axis, protein_coordinates)
    head_length = 3.0
    
    # Extend shaft backwards by 10%
    shaft_extension = base_shaft_length * 0.1
    total_shaft_length = base_shaft_length + shaft_extension
    shaft_spacing = 1.0
    
    pdb_lines = []
    atom_id = start_atom_id
    
    # Use different chain/residue IDs for each arrow
    shaft_res_id = 100 + arrow_index * 50
    head_res_id = shaft_res_id + 20
    chain_id = chr(ord('A') + arrow_index)
    
    # === CREATE SHAFT ===
    n_shaft_atoms = int(total_shaft_length / shaft_spacing)
    for i in range(n_shaft_atoms):
        t = -(base_shaft_length/2 + shaft_extension) + i * shaft_spacing
        pos = point_on_axis + t * screw_axis
        
        pdb_line = f"ATOM  {atom_id:5d}  CA  SHF {chain_id}{shaft_res_id:4d}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00 30.00           C"
        pdb_lines.append(pdb_line)
        atom_id += 1
        
    pdb_lines.append("TER")
    
    # === CREATE ARROW HEAD ===
    # Create perpendicular vectors for cone base
    if abs(screw_axis[0]) < 0.9:
        temp_vec = np.array([1.0, 0.0, 0.0])
    else:
        temp_vec = np.array([0.0, 1.0, 0.0])
    
    perp1 = np.cross(screw_axis, temp_vec)
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(screw_axis, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)
    
    # Arrow head overlaps with shaft end
    head_overlap = 2.0
    head_start_pos = point_on_axis + (base_shaft_length/2 - head_overlap) * screw_axis
    
    # Create cone structure
    n_layers = 4
    n_points_per_layer = 6
    
    for layer in range(n_layers + 1):
        layer_progress = layer / n_layers
        layer_pos = head_start_pos + layer_progress * head_length * screw_axis
        
        max_radius = 1.2
        layer_radius = max_radius * (1.0 - layer_progress) ** 1.5
        
        if layer == n_layers:
            # Tip of arrow
            pdb_line = f"ATOM  {atom_id:5d}  CA  ARH {chain_id}{head_res_id:4d}    {layer_pos[0]:8.3f}{layer_pos[1]:8.3f}{layer_pos[2]:8.3f}  1.00 50.00           C"
            pdb_lines.append(pdb_line)
            atom_id += 1
        else:
            # Circular layer
            for point in range(n_points_per_layer):
                angle = 2 * np.pi * point / n_points_per_layer
                
                circle_pos = (layer_pos + 
                             layer_radius * np.cos(angle) * perp1 + 
                             layer_radius * np.sin(angle) * perp2)
                
                pdb_line = f"ATOM  {atom_id:5d}  CA  ARH {chain_id}{head_res_id:4d}    {circle_pos[0]:8.3f}{circle_pos[1]:8.3f}{circle_pos[2]:8.3f}  1.00 40.00           C"
                pdb_lines.append(pdb_line)
                atom_id += 1
    
    pdb_lines.append("TER")
    return pdb_lines

def write_complete_pymol_script(output_path, protein_1, protein_2_name, protein_2_chain, 
                              domains, analysis_pairs, global_reference_id, bending_residues, window_size,
                              pair_specific_results=None):
    """
    Generate single PyMOL script with hierarchical structure coloring and arrows
    Uses hierarchical analysis pairs instead of single fixed domain
    """
    
    folder_name = f"{protein_1.name}_{protein_1.chain_param}_{protein_2_name}_{protein_2_chain}"
    dir_path_str = f"{output_path}/{folder_name}"
    dir_path = Path(dir_path_str)
    if not dir_path.exists():
        dir_path.mkdir(parents=True)
    
    # Generate arrow PDB using hierarchical analysis pairs
    arrow_pdb_path = write_arrows_pdb(dir_path_str, folder_name, analysis_pairs, domains, protein_1, pair_specific_results)
    
    # Generate combined PyMOL script
    pml_file = f"{dir_path_str}/{folder_name}.pml"
    
    try:
        with open(pml_file, 'w') as fw:
            # === HEADER ===
            fw.write("# DynDom Hierarchical Visualization\n")
            fw.write("# Domain-colored structure with hierarchical screw axis arrows\n")
            fw.write("reinitialize\n")
            fw.write(f"load {folder_name}.pdb\n")
            fw.write("bg_color white\n")
            fw.write("color grey\n")
            fw.write("\n")
            
            # === STRUCTURE DOMAIN COLORING ===
            structure_commands = _generate_structure_coloring_commands_hierarchical(
                protein_1, domains, global_reference_id, bending_residues, window_size)
            fw.write('\n'.join(structure_commands))
            fw.write("\n")
            
            # === ARROW VISUALIZATION ===
            if arrow_pdb_path:
                arrow_pdb_name = os.path.basename(arrow_pdb_path)
                arrow_commands = _generate_arrow_commands_hierarchical(
                    analysis_pairs, domains, global_reference_id, arrow_pdb_name, output_path)
                fw.write('\n'.join(arrow_commands))
                fw.write("\n")
            
            # === FINAL SETTINGS ===
            fw.write("# === FINAL SETTINGS ===\n")
            fw.write("set stick_transparency, 0.0\n")
            fw.write("set stick_quality, 15\n")
            fw.write("zoom all\n")
            fw.write("orient\n")
            fw.write("\n")
            
            # Cleanup
            fw.write("# Cleanup selections\n")
            fw.write("delete bending_residues\n")
            fw.write("delete arrow_*\n")
            fw.write("\n")
            
            # Status messages for hierarchical system
            fw.write("print 'DynDom hierarchical visualization loaded!'\n")
            fw.write(f"print 'Global reference domain: {global_reference_id} ({DOMAIN_COLORS[0]})'\n")
            
            # Show analysis pairs
            for i, (moving_id, ref_id) in enumerate(analysis_pairs):
                moving_domain = domains[moving_id]
                ref_color = _get_hierarchical_domain_color_index(ref_id, global_reference_id)
                moving_color = _get_hierarchical_domain_color_index(moving_id, global_reference_id)
                fw.write(f"print 'Analysis pair {i+1}: Domain {moving_id} ({DOMAIN_COLORS[moving_color]}) relative to Domain {ref_id} ({DOMAIN_COLORS[ref_color]})'\n")
                fw.write(f"print '  Rotation: {moving_domain.rot_angle:.1f}°'\n")
            
        print(f"PyMOL script created: {pml_file}")
        return pml_file
        
    except Exception as e:
        traceback.print_exc()
        print(f"Error creating hierarchical PyMOL script: {e}")
        return None

def _generate_structure_coloring_commands_hierarchical(protein_1, domains, global_reference_id, bending_residues, window_size):
    """Generate PyMOL commands for hierarchical structure domain coloring"""
    commands = []
    commands.append("# === HIERARCHICAL DOMAIN STRUCTURE COLORING ===")
    
    mid_point = (window_size - 1) // 2
    util_res = protein_1.utilised_residues_indices
    polymer = protein_1.get_polymer()

    # Collect all bending residue indices
    all_bend_res_indices = []
    for b in bending_residues.values():
        for i in b:
            bb = i  # Segments are already in sliding window indices
            index = polymer[bb].seqid.num
            all_bend_res_indices.append(index)

    # Color each domain with its hierarchical color
    for domain in domains:
        color_index = _get_hierarchical_domain_color_index(domain.domain_id, global_reference_id)
        color = DOMAIN_COLORS[color_index]
        
        domain_res_reg = []
        for s in range(domain.segments.shape[0]):
            reg = []
            for i in range(domain.segments[s][0], domain.segments[s][1]+1):
                j = i  # Segments are already in sliding window indices
                index = util_res[j]
                res_num = polymer[index].seqid.num
                if res_num not in all_bend_res_indices:
                    reg.append(res_num)
            domain_res_reg.extend(group_continuous_regions(reg))
        
        # Generate selection commands for this domain
        for s in range(len(domain_res_reg)):
            if s == 0:
                commands.append(f"select domain_{domain.domain_id}, resi {domain_res_reg[s][0]}-{domain_res_reg[s][1]}")
            else:
                commands.append(f"select domain_{domain.domain_id}, domain_{domain.domain_id} + resi {domain_res_reg[s][0]}-{domain_res_reg[s][1]}")
        
        domain_type = "GLOBAL REFERENCE" if domain.domain_id == global_reference_id else "MOVING"
        commands.append(f"color {color}, domain_{domain.domain_id}  # Domain {domain.domain_id} ({domain_type})")
        commands.append("")

    # Color the bending residues (green)
    bend_res_groups = group_continuous_regions(all_bend_res_indices)
    
    if bend_res_groups:
        commands.append("# Color bending residues")
        for i, g in enumerate(bend_res_groups):
            commands.append(f"select bending_residues_{i+1}, resi {g[0]}-{g[1]}")
            commands.append(f"color green, bending_residues_{i+1}")
        commands.append("")

    # Additional PyMOL settings
    commands.append("set dash_gap, 0")
    commands.append("set dash_radius, 0.2")
    commands.append("")
    
    return commands


def _generate_arrow_commands_hierarchical(analysis_pairs, domains, global_reference_id, arrow_pdb_name, output_path):
    """Generate PyMOL commands for hierarchical arrow display"""
    commands = []
    commands.append("# === HIERARCHICAL SCREW AXIS ARROWS ===")
    commands.append(f'load {output_path}')
    commands.append(f"load {arrow_pdb_name}")
    commands.append("")
    commands.append("# Basic protein display")
    commands.append(f'hide everything, {output_path}')
    commands.append(f'show cartoon, {output_path}')
    commands.append(f'color gray80, {output_path}')
    commands.append("")
    commands.append(f"# Hide arrow atoms initially")
    commands.append("hide everything, " + os.path.splitext(arrow_pdb_name)[0])
    commands.append("")
    
    # Add arrow visualization for each analysis pair
    for i, (moving_domain_id, reference_domain_id) in enumerate(analysis_pairs):
        moving_domain = domains[moving_domain_id]
        
        # Get colors based on hierarchical structure
        reference_color_index = _get_hierarchical_domain_color_index(reference_domain_id, global_reference_id)
        moving_color_index = _get_hierarchical_domain_color_index(moving_domain_id, global_reference_id)
        
        shaft_color = DOMAIN_COLORS[reference_color_index]  # Shaft = reference domain color
        head_color = DOMAIN_COLORS[moving_color_index]      # Head = moving domain color
        
        # Use chain-specific selections to completely separate arrows
        chain_id = chr(ord('A') + i)
        shaft_res_id = 100 + i * 50
        head_res_id = shaft_res_id + 20
        
        commands.extend([
            f"# Arrow {i+1}: Domain {moving_domain_id} (moving) relative to Domain {reference_domain_id} (reference)",
            f"# Shaft color: {shaft_color} (reference domain), Head color: {head_color} (moving domain)",
            f"# Rotation: {moving_domain.rot_angle:.1f}°",
            "",
            f"# Select shaft and head atoms by chain and residue",
            f"select shaft_{i+1}, chain {chain_id} and resn SHF and resi {shaft_res_id}",
            f"select head_{i+1}, chain {chain_id} and resn ARH and resi {head_res_id}",
            "",
            f"# Display shaft as thick licorice stick (REFERENCE domain color: {shaft_color})",
            f"show sticks, shaft_{i+1}",
            f"color {shaft_color}, shaft_{i+1}",
            f"set stick_radius, 0.3, shaft_{i+1}",  # Thicker for licorice style
            "",
            f"# Display arrow head as clean cone (MOVING domain color: {head_color})",
            f"show sticks, head_{i+1}",
            f"color {head_color}, head_{i+1}",
            f"set stick_radius, 0.25, head_{i+1}",
            "",
            f"# Connect atoms ONLY within each section",
            f"bond shaft_{i+1}, shaft_{i+1}",
            f"bond head_{i+1}, head_{i+1}",
            "",
        ])
        
    commands.extend([
        "# Disable automatic bonding between different chains",
        "set auto_bond, 0",
        "",
        "# Make arrows more prominent",
        "set stick_transparency, 0.0",
        "set stick_quality, 15",
        "set sphere_quality, 3",
        "set surface_quality, 2",
        "",
        "# Final settings",
        "set depth_cue, 0",
        "set ray_shadows, 1",
        "set ray_shadow_decay_factor, 0.1",
        "",
        "# Better lighting for 3D arrow heads",
        "set ambient, 0.2",
        "set direct, 0.8",
        "set reflect, 0.5",
        "set shininess, 10",
        "",
        "# Clean up selections", 
        "delete shaft_*",
        "delete head_*",
        "",
    ])
    
    return commands


def _get_hierarchical_domain_color_index(domain_id, global_reference_id):
    """Get color index for domain based on hierarchical structure"""
    if domain_id == global_reference_id:
        return 0  # Global reference always gets blue (index 0)
    
    # For non-reference domains, assign colors sequentially starting from index 1
    # This ensures consistent coloring across visualization components
    return (domain_id + 1) % len(DOMAIN_COLORS)
    
def _generate_structure_coloring_commands(protein_1, domains, fixed_domain_id, bending_residues, window_size):
    """Generate PyMOL commands for structure domain coloring"""
    commands = []
    commands.append("# === DOMAIN STRUCTURE COLORING ===")
    
    mid_point = (window_size - 1) // 2
    fixed_dom_segments = domains[fixed_domain_id].segments
    util_res = protein_1.utilised_residues_indices
    polymer = protein_1.get_polymer()

    # Collect all bending residue indices
    all_bend_res_indices = []
    for b in bending_residues.values():
        for i in b:
            bb = i  # Segments are already in sliding window indices
            index = polymer[bb].seqid.num
            all_bend_res_indices.append(index)
    


    # Color the fixed domain (blue)
    fixed_dom_res_reg = []
    for s in range(fixed_dom_segments.shape[0]):
        reg = []
        for i in range(fixed_dom_segments[s][0], fixed_dom_segments[s][1]+1):
            j = i  # Segments are already in sliding window indices
            index = util_res[j]
            res_num = polymer[index].seqid.num
            if res_num not in all_bend_res_indices:
                reg.append(res_num)
        fixed_dom_res_reg.extend(group_continuous_regions(reg))
    
    # Generate selection commands for fixed domain
    for s in range(len(fixed_dom_res_reg)):
        if s == 0:
            commands.append(f"select fixed_domain, resi {fixed_dom_res_reg[s][0]}-{fixed_dom_res_reg[s][1]}")
        else:
            commands.append(f"select fixed_domain, fixed_domain + resi {fixed_dom_res_reg[s][0]}-{fixed_dom_res_reg[s][1]}")
    
    commands.append(f"color {DOMAIN_COLORS[0]}, fixed_domain")
    commands.append("")

    # Color the dynamic domains
    region_count = 1
    for domain in domains:
        if domain.domain_id == fixed_domain_id:
            continue
            
        dyn_dom_res_reg = []
        segments = domain.segments
        
        for s in range(segments.shape[0]):
            reg = []
            for i in range(segments[s][0], segments[s][1]+1):
                j = i  # Segments are already in sliding window indices
                index = util_res[j]
                res_num = polymer[index].seqid.num
                if res_num not in all_bend_res_indices:
                    reg.append(res_num)
            dyn_dom_res_reg.extend(group_continuous_regions(reg))
        
        # Generate selection commands for this dynamic domain
        for s in range(len(dyn_dom_res_reg)):
            if s == 0:
                commands.append(f"select moving_domain_{domain.domain_id}, resi {dyn_dom_res_reg[s][0]}-{dyn_dom_res_reg[s][1]}")
            else:
                commands.append(f"select moving_domain_{domain.domain_id}, moving_domain_{domain.domain_id} + resi {dyn_dom_res_reg[s][0]}-{dyn_dom_res_reg[s][1]}")
        
        # Use DOMAIN_COLORS array with bounds checking
        color_index = min(region_count, len(DOMAIN_COLORS) - 1)
        commands.append(f"color {DOMAIN_COLORS[color_index]}, moving_domain_{domain.domain_id}")
        commands.append("")
        region_count += 1

    # Color the bending residues (green)
    bend_res_groups = group_continuous_regions(all_bend_res_indices)
    
    if bend_res_groups:
        commands.append("# Color bending residues")
        for i, g in enumerate(bend_res_groups):
            commands.append(f"select bending_residues_{i+1}, resi {g[0]}-{g[1]}")
            commands.append(f"color green, bending_residues_{i+1}")
        commands.append("")

    # Additional PyMOL settings
    commands.append("set dash_gap, 0")
    commands.append("set dash_radius, 0.2")
    commands.append("")
    
    return commands
def _generate_arrow_commands(domains, fixed_domain_id, arrow_pdb_name, output_path):
    """Generate PyMOL commands for arrow display using domain data directly"""
    commands = []
    commands.append("# === SCREW AXIS ARROWS ===")
    commands.append(f'load {output_path}')
    commands.append(f"load {arrow_pdb_name}")
    commands.append("")
    commands.append("# Basic protein display")
    commands.append(f'hide everything, {output_path}')
    commands.append(f'show cartoon, {output_path}')
    commands.append(f'color gray80, {output_path}')
    commands.append("")
    commands.append(f"# Hide arrow atoms initially")
    commands.append("hide everything, " + os.path.splitext(arrow_pdb_name)[0])
    commands.append("")
    
    moving_domains = [d for d in domains if d.domain_id != fixed_domain_id]
    
    # Add arrow visualization for each domain
    for i, domain in enumerate(moving_domains):
        moving_domain_id = domain.domain_id
        
        # Correct assignment: shaft = fixed, head = moving
        shaft_color = DOMAIN_COLORS[0]  # Fixed domain always gets blue (index 0)
        head_color = DOMAIN_COLORS[1 + i]  # Moving domains get red, yellow, etc. (indices 1+)
        
        # Use chain-specific selections to completely separate arrows
        chain_id = chr(ord('A') + i)
        shaft_res_id = 100 + i * 50
        head_res_id = shaft_res_id + 20
        
        commands.extend([
            f"# Arrow {i+1}: Domain {moving_domain_id} (moving) relative to Domain {fixed_domain_id} (fixed)",
            f"# Shaft color: {shaft_color} (fixed domain), Head color: {head_color} (moving domain)",
            f"# Rotation: {domain.rot_angle:.1f}°",
            "",
            f"# Select shaft and head atoms by chain and residue",
            f"select shaft_{i+1}, chain {chain_id} and resn SHF and resi {shaft_res_id}",
            f"select head_{i+1}, chain {chain_id} and resn ARH and resi {head_res_id}",
            "",
            f"# Display shaft as thick licorice stick (FIXED domain color: {shaft_color})",
            f"show sticks, shaft_{i+1}",
            f"color {shaft_color}, shaft_{i+1}",
            f"set stick_radius, 0.3, shaft_{i+1}",  # Thicker for licorice style
            "",
            f"# Display arrow head as clean cone (MOVING domain color: {head_color})",
            f"show sticks, head_{i+1}",
            f"color {head_color}, head_{i+1}",
            f"set stick_radius, 0.25, head_{i+1}",
            "",
            f"# Connect atoms ONLY within each section",
            f"bond shaft_{i+1}, shaft_{i+1}",
            f"bond head_{i+1}, head_{i+1}",
            "",
        ])
        
    commands.extend([
        "# Disable automatic bonding between different chains",
        "set auto_bond, 0",
        "",
        "# Make arrows more prominent",
        "set stick_transparency, 0.0",
        "set stick_quality, 15",
        "set sphere_quality, 3",
        "set surface_quality, 2",
        "",
        "# Final settings",
        "set depth_cue, 0",
        "set ray_shadows, 1",
        "set ray_shadow_decay_factor, 0.1",
        "",
        "# Better lighting for 3D arrow heads",
        "set ambient, 0.2",
        "set direct, 0.8",
        "set reflect, 0.5",
        "set shininess, 10",
        "",
        "# Clean up selections", 
        "delete shaft_*",
        "delete head_*",
        "",
    ])
    
    return commands
def write_w5_info_file(output_path, protein_1_name: str, chain_1, protein_2_name: str, chain_2, window, domain_size, ratio, atoms,
                       domains: list, analysis_pairs: list, global_reference_id: int, protein_1, pair_specific_results=None):
    """
    Write w5_info file with simplified hierarchical domain structure
    Shows fixed domain followed by moving domains relative to it
    """
    try:
        protein_folder = f"{protein_1_name}_{chain_1}_{protein_2_name}_{chain_2}"
        fw = open(f"{output_path}/{protein_folder}/{protein_folder}.w5_info", "w")
        fw.write("DynDom Python Version 1.0\n")
        fw.write(f"{protein_1_name}{chain_1}_{protein_2_name}{chain_2}.w5\n")
        fw.write(f"file name of conformer 1: {protein_1_name}.pdb\n")
        fw.write(f"chain id: {chain_1}\n")
        fw.write(f"file name of conformer 2: {protein_2_name}.pdb\n")
        fw.write(f"chain id: {chain_2}\n")
        fw.write(f"window length: {window}\n")
        fw.write(f"minimum ratio of external to internal motion: {ratio}\n")
        fw.write(f"minimum domain size: {domain_size}\n")
        fw.write(f"atoms to use: {atoms}\n")
        fw.write(f"THERE ARE {len(domains)} DOMAINS\n")
        fw.write("================================================================================\n")
        
        domain_colours = ["blue", "red", "yellow", "pink", "cyan", "purple", "orange", "brown", "black", "white", "magenta", "violet", "indigo", "turquoise", "coral"]
        
        # Group analysis pairs by reference domain
        reference_groups = {}
        for moving_id, reference_id in analysis_pairs:
            if reference_id not in reference_groups:
                reference_groups[reference_id] = []
            reference_groups[reference_id].append(moving_id)
        
        # Write each reference domain followed by its moving domains
        for reference_id in sorted(reference_groups.keys()):
            reference_domain = domains[reference_id]
            moving_domain_ids = reference_groups[reference_id]
            
            # Write fixed domain header
            if reference_id == global_reference_id:
                fw.write("FIXED DOMAIN - Fixed for pymol visualisation\n")
            else:
                fw.write("FIXED DOMAIN\n")
            color_index = _get_domain_color_index(reference_id, global_reference_id)
            fw.write(f"DOMAIN NUMBER: \t {reference_id + 1} (coloured {domain_colours[color_index]} for rasmol)\n")
            
            # Write domain residue ranges
            residue_str = _format_domain_residues(reference_domain, protein_1)
            fw.write(f"RESIDUE NUMBERS: \t{residue_str}\n")
            fw.write(f"SIZE: \t{reference_domain.num_residues} RESIDUES\n")
            fw.write(f"BACKBONE RMSD ON THIS DOMAIN: \t{round(reference_domain.rmsd, 3)}A\n")
            fw.write("------------------------------------------------------------------------------\n")
            
            # Write each moving domain for this reference
            for moving_id in moving_domain_ids:
                moving_domain = domains[moving_id]
                color_index = _get_domain_color_index(moving_id, global_reference_id)
                
                fw.write(f"MOVING DOMAIN (RELATIVE TO Domain {reference_id + 1})\n")
                fw.write(f"DOMAIN NUMBER: \t {moving_id + 1} (coloured {domain_colours[color_index]} for rasmol)\n")
                
                # Write domain residue ranges
                residue_str = _format_domain_residues(moving_domain, protein_1)
                fw.write(f"RESIDUE NUMBERS: \t{residue_str}\n")
                fw.write(f"SIZE: \t{moving_domain.num_residues} RESIDUES\n")
                
                # Use pair-specific results if available
                pair_key = (moving_id, reference_id)
                if pair_specific_results and pair_key in pair_specific_results:
                    pair_data = pair_specific_results[pair_key]
                    rot_angle = pair_data['rot_angle']
                    translation = pair_data['translation']
                    screw_axis = pair_data['screw_axis']
                    point_on_axis = pair_data['point_on_axis']
                    rmsd = pair_data['rmsd']
                    
                    # Calculate ratio (if available in pair data, otherwise use domain object)
                    if 'ratio' in pair_data:
                        ratio_value = pair_data['ratio']
                    else:
                        ratio_value = moving_domain.ratio
                else:
                    # Fallback to domain object if pair-specific data not available
                    rot_angle = moving_domain.rot_angle
                    translation = moving_domain.translation
                    screw_axis = moving_domain.screw_axis
                    point_on_axis = moving_domain.point_on_axis
                    rmsd = moving_domain.rmsd
                    ratio_value = moving_domain.ratio
                
                fw.write(f"BACKBONE RMSD ON THIS DOMAIN: \t{round(rmsd, 3)}A\n")
                
                # Write motion parameters using pair-specific data
                fw.write(f"RATIO OF INTERDOMAIN TO INTRADOMAIN DISPLACEMENT: \t{round(ratio_value, 3)}\n")
                fw.write(f"ANGLE OF ROTATION: \t{round(rot_angle, 3)} DEGREES\n")
                fw.write(f"TRANSLATION ALONG AXIS:\t{round(translation, 3)} A\n")
                fw.write(f"SCREW AXIS DIRECTION (UNIT VECTOR): \t{round(screw_axis[0], 3)} \t{round(screw_axis[1], 3)} \t{round(screw_axis[2], 3)}\n")
                fw.write(f"POINT ON AXIS: \t{round(point_on_axis[0], 3)} \t{round(point_on_axis[1], 3)} \t{round(point_on_axis[2], 3)}\n")
                
                # Write bending residues
                if hasattr(moving_domain, 'bend_res') and moving_domain.bend_res:
                    groups = group_continuous_regions(moving_domain.bend_res)
                    for group in groups:
                        start_pdb_num, end_pdb_num = _convert_bend_res_to_pdb_nums(group, protein_1)
                        fw.write(f"BENDING RESIDUES: \t{start_pdb_num} - {end_pdb_num}\n")
                
                fw.write("------------------------------------------------------------------------------\n")
        
        fw.close()
        
    except Exception as e:
        print(f"Error writing w5_info file: {e}")
        traceback.print_exc()
        return False
    return True


def _format_domain_residues(domain, protein_1):
    """Helper function to format domain residue numbers for output"""
    slide_window_indices = protein_1.slide_window_residues_indices
    util_res = protein_1.utilised_residues_indices
    polymer = protein_1.get_polymer()
    
    residue_str = ""
    for s in range(domain.segments.shape[0]):
        # Convert sequence indices to PDB residue numbers
        start_seq_idx = domain.segments[s][0]
        end_seq_idx = domain.segments[s][1]
        
        start_util_pos = slide_window_indices[0] + start_seq_idx
        end_util_pos = slide_window_indices[0] + end_seq_idx
        
        start_polymer_idx = util_res[start_util_pos]
        end_polymer_idx = util_res[end_util_pos]
        
        start_pdb_num = polymer[start_polymer_idx].seqid.num
        end_pdb_num = polymer[end_polymer_idx].seqid.num
        
        if s == 0:
            residue_str = f"{start_pdb_num} - {end_pdb_num}"
        else:
            residue_str = residue_str + f", {start_pdb_num} - {end_pdb_num}"
    
    return residue_str


def _get_domain_color_index(domain_id, global_reference_id):
    """Get color index for domain based on hierarchical structure"""
    if domain_id == global_reference_id:
        return 0  # Global reference always gets blue (index 0)
    
    # For non-reference domains, assign colors sequentially
    # This is a simplified approach - you might want to make this more sophisticated
    return (domain_id + 1) % 15  # Cycle through available colors


def _convert_bend_res_to_pdb_nums(bend_res_group, protein_1):
    """Convert bending residue sequence indices to PDB residue numbers"""
    slide_window_indices = protein_1.slide_window_residues_indices
    util_res = protein_1.utilised_residues_indices
    polymer = protein_1.get_polymer()
    
    start_seq_idx = bend_res_group[0]
    end_seq_idx = bend_res_group[-1]
    
    start_util_pos = slide_window_indices[0] + start_seq_idx
    end_util_pos = slide_window_indices[0] + end_seq_idx
    
    start_polymer_idx = util_res[start_util_pos]
    end_polymer_idx = util_res[end_util_pos]
    
    start_pdb_num = polymer[start_polymer_idx].seqid.num
    end_pdb_num = polymer[end_polymer_idx].seqid.num
    
    return start_pdb_num, end_pdb_num


def group_continuous_regions(data: list):
    groups = []
    for k, g in groupby(enumerate(data), lambda ix: ix[0] - ix[1]):
        temp = list(map(itemgetter(1), g))
        groups.append([temp[0], temp[-1]])
    return groups


