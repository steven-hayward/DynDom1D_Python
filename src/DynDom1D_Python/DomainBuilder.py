import gemmi
import numpy as np
from sklearn.neighbors import KDTree
from Domain import Domain


def domain_builder(protein_1, protein_2, segments_dict: dict, min_domain: int):
    """
    Builds the domains from each cluster in both protein chains.
    A Binary Connection Matrix is first created in each cluster that determines which segments of the chain are
    connected with another segment in the same cluster.
    :return: List of Domain objects
    """
    domains_1 = []
    domains_2 = []
    break_cluster = False
    protein_1_slide_win_res = protein_1.get_slide_window_residues()
    protein_2_slide_win_res = protein_2.get_slide_window_residues()
    # For each cluster and their segments
    for cluster, cluster_segments in segments_dict.items():
        # https://journals.iucr.org/paper?buy=yes&cnor=jn0027&showscheme=yes&sing=yes
        # Create a binary connection matrix for both proteins from the cluster segments
        bin_mat_1 = create_bin_conn_mat(protein_1_slide_win_res, cluster_segments)
        bin_mat_2 = create_bin_conn_mat(protein_2_slide_win_res, cluster_segments)
        # List of Numpy arrays of the binary connection matrix after reduction
        reduced_mat_1 = row_reduction(bin_mat_1)
        reduced_mat_2 = row_reduction(bin_mat_2)
        domains_1, contains_valid_domains_1 = create_cluster_domains(cluster_segments=cluster_segments,
                                                                     domains=domains_1, reduced_mat=reduced_mat_1,
                                                                     min_domain_size=min_domain, cluster_id=cluster)
        domains_2, contains_valid_domains_2 = create_cluster_domains(cluster_segments=cluster_segments,
                                                                     domains=domains_2, reduced_mat=reduced_mat_2,
                                                                     min_domain_size=min_domain, cluster_id=cluster)
        # If both proteins do not have at least one domain with the minimum number of residues for the cluster,
        # the condition is not met
        if not (contains_valid_domains_1 or contains_valid_domains_2):
            break_cluster = True
            return domains_1, break_cluster

    # For any domain that does not meet the minimum number of residues, take the segments of the tiny domain and put them
    # into a domain that is larger.
    while True:
        small_domain = None
        for domain in domains_1:
            if domain.num_residues < min_domain:
                small_domain = domain
                break
        if small_domain is None:
            break
        domains_1 = join_domains(small_domain, domains_1)

    return domains_1, break_cluster


def create_bin_conn_mat(residues: list, segments: np.array):
    """
    Creates a binary connection matrix that determines which segments are connected with each other.
    :param residues:    List of gemmi.Residue objects from the slide window
    :param segments:    Numpy array of arrays [s, e], where s and e are the indexes where the segment starts and end in the
                        residues list respectively
    :return bin_conn_mat: A matrix of size NxN where N is the length of the segments array
    """
    # Binary connection matrix that stores which segments are connected
    bin_conn_mat = np.zeros((segments.shape[0], segments.shape[0]), dtype=int)
    curr_segment = 0
    # Go through each segment.
    while curr_segment < segments.shape[0]:
        # Get the start and end indices of the residues in the current segment
        curr_segment_start_index = segments[curr_segment][0]
        curr_segment_end_index = segments[curr_segment][1] + 1
        # Check whether current residues segment is connected with any other residues segment
        for s in range(segments.shape[0]):
            # No need to check distances between identical segments.
            if curr_segment == s:
                bin_conn_mat[curr_segment][s] = 1
                continue
            # Get the start and end indices of the segment to be checked
            checked_segment_start_index = segments[s][0]
            checked_segment_end_index = segments[s][1] + 1
            # Get coordinates of all atoms from every residue in the to-be-checked segment
            checked_segment_atoms_coordinates = get_segment_atoms_coordinates(residues, checked_segment_start_index, checked_segment_end_index)
            # Checks the distances between the segments using all the atoms
            for r in range(curr_segment_start_index, curr_segment_end_index):
                # Get the coordinates of all atoms in one residue of the current segment
                residue_atoms_coordinates, has_side_chain = get_residue_atoms_coordinates(residues[r])
                # Check if the atoms are close to any of the atoms of the to-be-checked segment
                is_connected = check_connection_between_atoms(residue_atoms_coordinates, checked_segment_atoms_coordinates,
                                                              has_side_chain)
                # If atoms are close enough, the segments are connected.
                if is_connected:
                    bin_conn_mat[curr_segment][s] = 1
                    bin_conn_mat[s][curr_segment] = 1
                    break
        curr_segment += 1
    return bin_conn_mat


def row_reduction(matrix):
    """
    Performs row reduction on the binary connection matrix to obtain the groups of connected segments. The method uses
    the OR-wise logical operation between the rows containing connections, represented by 1.
    :param matrix: NxN Numpy array
    :return: A list of Numpy arrays containing indices indicating the location of the segments in the cluster.
             If there are 5 segments in the cluster, Segment index = [0, 1, 2, 3, 4], and the segments connected with
             each other are (0, 3), (1), and (2, 4). The list will be returned as [array([0, 3]), array([1]), array([2, 4])]
    """
    rows = matrix.shape[0]
    # This list will store Numpy arrays each containing groups of indices of the segments. These groups correspond to
    # segments connected with each other
    connections_list: list = []
    # Go through each row sequentially starting from the row with the highest index to the lowest index
    for m in range(rows-1, -1, -1):
        # The current row will be used to perform OR-wise operations on other rows
        or_wise_result = matrix[m]
        # Go through each index (column) sequentially starting from the current row number to the lowest index
        for n in range(m, -1, -1):
            # OR-wise operations must only be done between rows where the location (m, n) in the matrix is 1
            if matrix[m][n] == 1:
                # Obtain the array at row n
                list_to_or_wise = matrix[n]
                # Perform OR wise logical operation between the 2 arrays
                or_wise_result |= list_to_or_wise
                # Returns tuple (a, b) where a is the list of indices of or_wise_result where the element is 1.
                # b is empty.
                results = np.where(or_wise_result == 1)
                indices_array = results[0]
                # ind = [i for i in range(len(connections_list)) if any(np.isin(indices_array, connections_list[i]))]
                # This is the multi-lined equivalent. Explanation is shown below.
                ind = []
                # Go through the connections_list. The first run will always be skipped since connections_list is
                # empty.
                for c in range(len(connections_list)):
                    # Checks if any values in indices_array exists in the Numpy array of connections_list at index c.
                    if any(np.isin(indices_array, connections_list[c])):
                        # If yes, append the index to ind
                        ind.append(c)
                # If ind is not empty
                if len(ind) > 0:
                    # Get the Numpy array at index c of connections_list and append the indices_array to it.
                    temp = np.append(connections_list[ind[0]], indices_array)
                    # Get the unique values from temp and set it as the new values for connections_list at index c
                    connections_list[ind[0]] = np.unique(temp, return_counts=False)
                # For the first run (connections_list is empty), just add the first row to the list.
                else:
                    connections_list.append(indices_array)

    return connections_list


def create_cluster_domains(cluster_segments: list, domains: list, reduced_mat: list, min_domain_size: int,
                           cluster_id: int):
    """
    Create domains based on the reduced binary connection matrix
    :param cluster_segments: The segments of the protein chain for the cluster
    :param domains: The list of Domain objects
    :param reduced_mat: The reduced binary connection matrix from connected set algorithm
    :param min_domain_size: Minimum number of residues for a domain
    :param cluster_id:
    :return:
    """
    # Create a copy of the Domain list
    new_domains = domains
    contains_valid_domains = False
    # For each row in the reduction matrix
    for rows in range(len(reduced_mat)):
        # Get the segment from the cluster
        segments = np.array([cluster_segments[r] for r in reduced_mat[rows]])
        # Create a Domain object
        domain = Domain(cluster_id, len(domains), segments)
        # Checks if at least one domain size meets the minimum number of residues
        if domain.num_residues >= min_domain_size:
            contains_valid_domains = True
        new_domains.append(domain)
    return new_domains, contains_valid_domains


def join_domains(tiny_domain: Domain, domains: list):
    """
    Takes the tiny domain's (smaller than minimum domain size) segments and absorbs it into another domain.
    :param tiny_domain: Domain that is smaller than minimum domain size
    :param domains: List of Domain objects
    :return new_domains: Updated List of Domain objects
    """
    tiny_domain_segments = tiny_domain.segments

    # Goes through each segment of the tiny domain
    for ts in range(tiny_domain_segments.shape[0]):
        # First obtain the indices of the previous and next residues connected to the segments
        prev_connecting_index = tiny_domain_segments[ts][0] - 1
        next_connecting_index = tiny_domain_segments[ts][1] + 1
        # If the previous index of the segment is -1, it means the segment is at the tail end of the protein chain.
        if prev_connecting_index == -1:
            for curr_d in domains:
                if next_connecting_index in curr_d.segments:
                    curr_d.add_segment(tiny_domain_segments[ts], add_at_left_side=True)
                    break
            continue
        # Going through the other domains
        for d in range(len(domains)):
            # Ignore any domains that are from the same cluster as the tiny domain. Ignore the identical domain as well.
            if domains[d].cluster_id == tiny_domain.cluster_id or domains[d].domain_id == tiny_domain.domain_id:
                continue

            if prev_connecting_index in domains[d].segments:
                domains[d].add_segment(tiny_domain_segments[ts])
                break

    new_domains = []
    domain_count = 0
    # Remove the tiny domains and reset the domain IDs
    for domain in domains:
        if domain.domain_id == tiny_domain.domain_id:
            continue
        domain.domain_id = domain_count
        new_domains.append(domain)
        domain_count += 1

    return new_domains


def remove_domains(domains: list, domains_to_remove: list):
    """
    Removes domains that smaller than minimum domain size
    :param domains: List of Domains
    :param domains_to_remove: IDs of the to-be-removed domains in the list
    :return: new_domains
    """
    new_domains = []
    domains_added = 0
    for d in range(len(domains)):
        if domains[d].domain_id in domains_to_remove:
            continue
        domains[d].domain_id = domains_added
        new_domains.append(domains[d])
        domains_added += 1

    return new_domains


def get_segment_atoms_coordinates(residues: list, start_index, end_index):
    """
    Gets the residues' atom coordinates
    :param residues: List of gemmi.Residue objects
    :param start_index:
    :param end_index:
    :return:
    """

    pos = []
    for i in range(start_index, end_index):
        for a in residues[i]:
            pos.append(a.pos.tolist())
    atoms_pos = np.asarray(pos)
    return atoms_pos


def get_residue_atoms_coordinates(residue: gemmi.Residue):
    """
    Get coordinates of all atoms of given residue and checks whether residue contains side chains (Residues
    containing more than just N, CA, C, and O).
    :param residue: A gemmi.Residue object
    :return:
    """
    side_chain = False
    # If residue has more than 4 atoms, we can assume that there are side chains
    for atom in residue:
        if atom.name not in ["N", "CA", "C", "O"]:
            side_chain = True
            break
    atoms_pos = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in residue])
    return atoms_pos, side_chain


def check_connection_between_atoms(residue_atoms: np.array, segment_atoms: np.array, side_chain=False):
    """
    Checks if any atoms of a residue are close to a segment
    :param residue_atoms:
    :param segment_atoms:
    :param side_chain:
    :return:
    """
    dist_criteria = 4 if side_chain else 10
    kdt = KDTree(segment_atoms, metric="euclidean")
    count = kdt.query_radius(residue_atoms, r=dist_criteria, count_only=True)
    return True if any(c > 0 for c in count) else False

