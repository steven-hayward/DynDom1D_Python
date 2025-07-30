from datetime import datetime
import math
import numpy as np
import gemmi
import Protein
import DomainBuilder as dom_build
from Domain import Domain
from sklearn.cluster import KMeans
from itertools import groupby
from operator import itemgetter
from clustering_logger import ClusteringLogger
from HierarchySystem import HierarchicalDomainSystem


class Clusterer:
    
    def __init__(self, protein_1: Protein, protein_2: Protein, rotation_vectors: np.array, atoms_to_use, k_means_n_init=1,
             k_means_max_iter=500, window=5, domain=20, ratio=1.0, enable_logging=True, log_dir="clustering_logs"):
        self.protein_1: Protein = protein_1
        self.protein_2: Protein = protein_2
        self.k_means_n_init = k_means_n_init
        self.k_means_max_iter = k_means_max_iter
        self.window = window
        self.domain = domain
        self.ratio = ratio
        # We want a large enough k so that it doesn't get used. If it does manage to get used, good luck. 
        self.max_k = 40
        # Start with 2 clusters to do the clustering
        self.current_k = 2
        # The k value that is valid for domain building (Number of clusters)
        self.valid_k = 0
        # The list of atoms to be used
        self.atoms_to_use = atoms_to_use
        # The number of atoms used in each residue
        self.num_atoms = len(self.atoms_to_use)
        # The rotation vectors of the sliding windows
        self.rotation_vectors: np.array = rotation_vectors
        # A dictionary containing valid clusters to be used for the domain calculation
        self.valid_clusters = {}
        # The clustering result from KMeans
        self.k_means_results = None
        self.segments = {}
        self.domains = []
        # The domain ID number determining which domain will be the fixed point
        self.fixed_domain = None
        self.clusterer_status = 1
        # This condition determines whether a k-value cluster is valid based on number of residues
        self.valid_cluster_found = False
        # This condition determines whether a valid cluster contains domains that have valid domain pair ratio
        self.valid_domains_found = False
        
        # MODIFIED: Only create logger if enabled, otherwise set to None (will be assigned externally)
        if enable_logging:
            protein_name = f"{protein_1.name}_{protein_1.chain_param}_{protein_2.name}_{protein_2.chain_param}"
            self.logger = ClusteringLogger(log_dir, protein_name)
        else:
            self.logger = None  # Will be assigned externally


    def cluster(self):
        fails = 0
        clusters = {0: self.rotation_vectors}
        while self.current_k < self.max_k:
            centroids = self.calc_cluster_centroids(clusters)
            # KMeans the rotation vectors to obtain k number of clusters of rotation vectors
            self.k_means_results, clusters = self.calc_k_means(centroids)
            # Obtain the segments from the KMeans labels.
            # The segments' indices are from the slide windowed residues. So the first and last few indices are not included.
            temp_segments, cluster_residues_small = self.determine_segments()


            # Prepare variables for logging
            temp_domains = None
            cluster_break = False
            ratios = None
            ratio_not_met = False
            decision = None
            reason = None
            print(f'\nWindow size: {self.window}, k: {self.current_k}, fails: {fails} ')
            # If there is a cluster where its total residue is smaller than min domain size:
            if cluster_residues_small and self.valid_domains_found:
                # If there is a previous iteration where valid clusters and valid domain pair ratios are found, clustering
                # can be halted
                print("Found cluster with number of residues smaller than minimum domain size. Clustering halted. ")
                print(f"Final clusterer to use: {self.valid_k}")
                self.clusterer_status = 0
                
                # Log this final attempt
                if self.logger:
                    self.logger.log_k_attempt(
                        k=self.current_k,
                        window_size=self.window,
                        domains=None,
                        cluster_segments=temp_segments,
                        cluster_residues_small=cluster_residues_small,
                        domain_build_failed=False,
                        ratios=None,
                        ratio_decision=None,
                        decision='halt_success',
                        reason=f'Cluster too small but previous valid solution exists (k={self.valid_k})'
                    )

                if self.logger:
                    self.logger.log_final_result(
                        final_k=self.valid_k if hasattr(self, 'valid_k') else None,
                        final_domains=self.domains if hasattr(self, 'domains') else None,
                        final_window=self.window,
                        clusterer_status=self.clusterer_status,
                        total_attempts=None
                    )
                return
                
            if cluster_residues_small and not self.valid_cluster_found:
                # If there are no previous iterations where domains have been found yet, skip to the next k value
                print("Found cluster with number of residues smaller than minimum domain size. Going to next k value")
                self.current_k += 1
                fails += 1
                
                # Log this attempt
                if self.logger:
                    self.logger.log_k_attempt(
                        k=self.current_k - 1,
                        window_size=self.window,
                        domains=None,
                        cluster_segments=temp_segments,
                        cluster_residues_small=cluster_residues_small,
                        domain_build_failed=False,
                        ratios=None,
                        ratio_decision=None,
                        decision='reject_size_early',
                        reason=f'Cluster too small, no previous valid solutions (fails={fails})'
                    )
                
                if fails > 5:
                    print("Too many fails. Increasing window size by 2")
                    self.clusterer_status = -1
                    
                    # Log window size change
                    if self.logger:
                        self.logger.log_window_change(
                            old_window=self.window,
                            new_window=self.window + 2,
                            reason=f'Too many consecutive failures ({fails})'
                        )

                    if self.logger:
                        self.logger.log_final_result(
                            final_k=None,
                            final_domains=None,
                            final_window=self.window,
                            clusterer_status=self.clusterer_status,
                            total_attempts=None
                        )
                    return
                continue
                
            if cluster_residues_small and self.valid_cluster_found:
                # If there was a previous iteration where a cluster had domains but the domain pair ratios are invalid,
                # reset the clustering but this time increase the window size by 2.
                print("Found cluster with number of residues smaller than minimum domain size. Valid cluster previously "
                    "found but no domains. Increasing window size by 2")
                self.clusterer_status = -1
                
                # Log this attempt
                if self.logger:
                    self.logger.log_k_attempt(
                        k=self.current_k,
                        window_size=self.window,
                        domains=None,
                        cluster_segments=temp_segments,
                        cluster_residues_small=cluster_residues_small,
                        domain_build_failed=False,
                        ratios=None,
                        ratio_decision=None,
                        decision='reject_size_restart',
                        reason='Cluster too small, had valid clusters before but no valid domains'
                    )

                if self.logger:
                    self.logger.log_window_change(
                        old_window=self.window,
                        new_window=self.window + 2,
                        reason='Valid clusters found before but no valid domains'
                    )
                
                # Log final result before returning
                if self.logger:
                    self.logger.log_final_result(
                        final_k=None,
                        final_domains=None,
                        final_window=self.window,
                        clusterer_status=self.clusterer_status,
                        total_attempts=None
                    )
                return
                
            self.valid_cluster_found = True
            # Construct the domains from the clusters
            temp_domains, cluster_break = dom_build.domain_builder(self.protein_1, self.protein_2,
                                                                temp_segments, self.domain)

            if cluster_break and self.valid_domains_found:
                print("All domains found are smaller than minimum domain size. Clustering halted.")
                self.clusterer_status = 0
                
                # Log this attempt
                if self.logger:
                    self.logger.log_k_attempt(
                        k=self.current_k,
                        window_size=self.window,
                        domains=None,
                        cluster_segments=temp_segments,
                        cluster_residues_small=cluster_residues_small,
                        domain_build_failed=cluster_break,
                        ratios=None,
                        ratio_decision=None,
                        decision='halt_domain_build',
                        reason='Domain building failed but previous valid solution exists'
                    )

                if self.logger:
                    self.logger.log_final_result(
                        final_k=self.valid_k if hasattr(self, 'valid_k') else None,
                        final_domains=self.domains if hasattr(self, 'domains') else None,
                        final_window=self.window,
                        clusterer_status=self.clusterer_status,
                        total_attempts=None
                    )
                return
                
            elif cluster_break and not self.valid_domains_found:
                print("All domains in the cluster are less than minimum domain size. Increasing k by 1")
                self.current_k += 1
                
                # Log this attempt
                if self.logger:
                    self.logger.log_k_attempt(
                        k=self.current_k - 1,
                        window_size=self.window,
                        domains=None,
                        cluster_segments=temp_segments,
                        cluster_residues_small=cluster_residues_small,
                        domain_build_failed=cluster_break,
                        ratios=None,
                        ratio_decision=None,
                        decision='reject_domain_build',
                        reason='Domain building failed - all domains too small after connectivity analysis'
                    )
                continue
              
            ratio_not_met = not self.check_ratios_hierarchical(temp_domains)

            if ratio_not_met:
                self.current_k += 1
                
                # Log this attempt
                if self.logger:
                    self.logger.log_k_attempt(
                        k=self.current_k - 1,
                        window_size=self.window,
                        domains=temp_domains,
                        cluster_segments=temp_segments,
                        cluster_residues_small=cluster_residues_small,
                        domain_build_failed=cluster_break,
                        ratios=None,  # Could collect ratios from hierarchical analysis if needed
                        ratio_decision='fail',
                        decision='reject_ratios',
                        reason='Hierarchical ratio validation failed'
                    )
                continue
            else:
                # The first time that a k value produces valid domains for all clusters, set it to true
                self.valid_k = self.current_k
                self.domains = temp_domains
                
                # Store the hierarchical fixed domain (primary reference)
                self.fixed_domain = self.get_hierarchical_fixed_domain()
                
                self.segments = temp_segments
                self.valid_domains_found = True
                self.valid_clusters = clusters
                
                # Log this successful attempt
                if self.logger:
                    avg_ratio = 0  # Could calculate from hierarchical analysis if needed
                    self.logger.log_k_attempt(
                        k=self.current_k,
                        window_size=self.window,
                        domains=temp_domains,
                        cluster_segments=temp_segments,
                        cluster_residues_small=cluster_residues_small,
                        domain_build_failed=cluster_break,
                        ratios=None,  # Could collect ratios from hierarchical analysis if needed
                        ratio_decision='pass',
                        decision='accept',
                        reason=f'Hierarchical analysis successful (domains: {len(temp_domains)})'
                    )
                    
            self.current_k += 1

        # If we exit the while loop without returning, log the final result
        if self.logger:
            self.logger.log_final_result(
                final_k=self.valid_k if hasattr(self, 'valid_k') else None,
                final_domains=self.domains if hasattr(self, 'domains') else None,
                final_window=self.window,
                clusterer_status=self.clusterer_status,
                total_attempts=None
            )
            
    def calc_cluster_centroids(self, clusters: dict):
        """
        Find the initial cluster centroids before performing KMeans clustering.
        :param clusters: A dictionary of clusters with the cluster IDs being the key, and array of indices being the value.
        :return:
        """
        # Stores the maximum variance
        max_var = 0
        # The dimension where the variance is largest (x, y, or z) -> (0, 1, or 2)
        dimension = None
        # The cluster ID with the largest variance
        max_var_cluster = None
        # The list of centroids for each cluster
        centroids = []
        # Each cluster
        for cluster in clusters:
            # Get the rotation vectors
            rotation_vectors: np.array = clusters[cluster]
            # Calculate the average of the cluster's rotation vectors. This is the centroid and append it to the list.
            centroids.append(np.average(rotation_vectors, axis=0))
            # Calculate the dimensional variance of the cluster's rotation vectors
            variances = np.var(rotation_vectors, axis=0)
            # If any of the variances is greater than the current max variance,
            if np.max(variances) > max_var:
                # Replace the current max variance with the new one
                max_var = np.max(variances)
                # Assign the cluster ID to the max variance cluster
                max_var_cluster = cluster
                # Get the dimension which has the highest variance
                dimension = np.argmax(variances)

        # Remove the centroid of the cluster with the highest variance. It will be replaced with 2 new centroids
        del centroids[max_var_cluster]

        # The rotation vectors of the highest-variance cluster
        rot_vecs = clusters[max_var_cluster]
        # The average value of the rotation vectors of one dimension only (The dimension with the highest variance).
        avg_dim_val = np.average(rot_vecs[:, dimension])
        # Divide the rotation vectors using the average dimension value.
        below_avg_rot_vecs = rot_vecs[(rot_vecs[:, dimension] < avg_dim_val)]
        above_avg_rot_vecs = rot_vecs[(rot_vecs[:, dimension] >= avg_dim_val)]
        # Calculate the new centroids
        centroids.append(np.average(below_avg_rot_vecs, axis=0))
        centroids.append(np.average(above_avg_rot_vecs, axis=0))
        return centroids

    def calc_k_means(self, centroids="k-means++"):
        """
        Performs KMeans on the rotation vectors.
        :return: KMeans results
        """
        k_means = KMeans(
            n_clusters=self.current_k,
            n_init=self.k_means_n_init,
            max_iter=self.k_means_max_iter,
            init=centroids
        ).fit(self.rotation_vectors)
        clusters = {}

        # Store the rotation vectors corresponding to the clusters
        for i in range(self.current_k):
            clusters[i] = self.rotation_vectors[(k_means.labels_ == i)]

        return k_means, clusters

    def determine_segments(self):
        """
        Gets segments from the array of K-Means labels where the segments contain the same cluster ID labels in a row.
        :return:    A dictionary where the keys are the KMean label cluster IDs and the values are a 2D array of [start, end]
                    where start and end are the indices in the KMeans labels where the segment starts and end.
        """
        # Initialise the label checker with the first label
        current_label_to_check = self.k_means_results.labels_[0]
        # The starting index of the segment
        start_index = 0
        # The end index of the segment
        end_index = 0

        # Initialise the dictionary of arrays to store the segments' indices
        segment_indices = {key: np.array([[]], dtype=int) for key in range(self.current_k)}

        # Iterate through each label
        for i in range(len(self.k_means_results.labels_)):
            # When the label is not equal to the checker, it means the segment ends there and a new segment is to be
            # created.
            if self.k_means_results.labels_[i] != current_label_to_check:
                # The segment's start and end indices are obtained and stored.
                temp = np.append(segment_indices[current_label_to_check], [[start_index, end_index]],
                                 axis=0 if segment_indices[current_label_to_check].shape[1] > 0 else 1)
                segment_indices[current_label_to_check] = temp
                # Set the label checker with the new label found
                current_label_to_check = self.k_means_results.labels_[i]
                # Reset the segment's start index
                start_index = i
            end_index = i
            # This is used if the final label is reached
            if i == len(self.k_means_results.labels_) - 1:
                temp = np.append(segment_indices[current_label_to_check], [[start_index, end_index]],
                                 axis=0 if segment_indices[current_label_to_check].shape[1] > 0 else 1)
                segment_indices[current_label_to_check] = temp
                # segment_indices[current_element_to_check].append((start_index, i))

        cluster_residues_too_little = False
        # Checks whether each cluster's total residues have minimum domain size. If a cluster is less than minimum
        # domain size, it means that the clustering will stop.
        for indices in segment_indices.values():
            num_residues = indices[:, 1] + 1 - indices[:, 0]
            if sum(num_residues) < 20:
                cluster_residues_too_little = True
                break

        return segment_indices, cluster_residues_too_little

    def find_fixed_domain(self, domains):
        """
        Takes the domains of one of the protein conformation and determines which domain has the most number of domains
        connected to it.
        :return chosen_domain:  The index/id of the domain with the most connected domains. Domain ID is the same as the
                                index of the domain
        """
        # If there are only 2 domains, the fixed domain is the domain with the largest number of residues
        if len(domains) <= 2:
            domain_sizes = [sum(d.segments[:, 1] + 1 - d.segments[:, 0]) for d in domains]
            chosen_domain = np.argmax(domain_sizes)
            return chosen_domain
        # Get the segments from each domain
        segments = [d.segments for d in domains]
        # Prepare a dictionary that stores the domain ID as key and the value is an array that will hold other domain
        # IDs that the domain is connected to
        connectivity = {d: np.array([]) for d in range(len(domains))}
        # For each domain
        for curr_d in domains:
            # Get the indices that are before and after the indices of the domain's segments
            prev_indices = segments[curr_d.domain_id][:, 0] - 1
            next_indices = segments[curr_d.domain_id][:, 1] + 1
            # For each of the other domains
            for d in domains:
                if curr_d.domain_id == d.domain_id:
                    continue
                # Check if any before and after indices are in the segment list of the domain
                prev_hits = np.in1d(prev_indices, segments[d.domain_id])
                next_hits = np.in1d(next_indices, segments[d.domain_id])
                if np.any([prev_hits, next_hits]):
                    temp = np.append(connectivity[curr_d.domain_id], d.domain_id)
                    connectivity[curr_d.domain_id] = temp
                    temp = np.append(connectivity[d.domain_id], curr_d.domain_id)
                    connectivity[d.domain_id] = temp
        # Can there be multiple max connectivities?
        # https://datagy.io/python-get-dictionary-key-with-max-value/
        connectivity = {key: np.unique(value).size for key, value in connectivity.items()}
        chosen_domain = max(connectivity, key=connectivity.get)
        return chosen_domain

    def mass_weighted_fit(self, dynamic_domain: Domain, fixed_domain: Domain):
        """
        Performs a mass-weighted protein best fit on a domain pair to get a new set of coordinates of a domain pair.
        :param dynamic_domain: The domain connected to the fixed domain
        :param fixed_domain: The fixed domain
        # :return r_1:    A SupResult object containing RMSD, Center 1, Center 2, Rotation Matrix, and Translation Vector
        #                 of Protein 1
        :return slide_window_1: The slide window residue chain of Protein 1 after fitting (transformation)
                                to Protein 2's position
        """
        slide_window_1 = self.protein_1.get_slide_window_residues()
        slide_window_2 = self.protein_2.get_slide_window_residues()
        # Lists to store the backbone atom coordinates of dynamic and fixed domain residues together in Protein 1 and 2.
        coords_1 = []
        coords_2 = []
        # The weights for each atom
        weights = []

        # Get number of residues in the dynamic and fixed domains and calculate the weights for the atoms.
        dynamic_domain_num_atoms = dynamic_domain.num_residues * len(self.atoms_to_use)
        fixed_domain_num_atoms = fixed_domain.num_residues * len(self.atoms_to_use)
        dynamic_weight = 1 / dynamic_domain_num_atoms
        fixed_weight = 1 / fixed_domain_num_atoms

        # In the dynamic domain, get the coordinates of the atoms in the segments and their weights and store into the lists
        for s in dynamic_domain.segments:
            for i in range(s[0], s[1] + 1):
                for a in self.atoms_to_use:
                    coords_1.append(slide_window_1[i][a][0].pos)
                    coords_2.append(slide_window_2[i][a][0].pos)
                    weights.append(dynamic_weight)

        # Same with the fixed domain
        for s in fixed_domain.segments:
            for i in range(s[0], s[1] + 1):
                for a in self.atoms_to_use:
                    coords_1.append(slide_window_1[i][a][0].pos)
                    coords_2.append(slide_window_2[i][a][0].pos)
                    weights.append(fixed_weight)

        # The superposition result of Protein 2 fitting onto Protein 1
        r: gemmi.SupResult = gemmi.superpose_positions(coords_1, coords_2, weights)
        slide_window_2.transform_pos_and_adp(r.transform)
        # return r_1, r_2, slide_window_1, slide_window_2
        return r, slide_window_2

    def fit_fixed_dom(self, dynamic_domain: Domain, fixed_domain: Domain):
        print(f"Fitting fixed domain {fixed_domain.domain_id} onto dynamic domain {dynamic_domain.domain_id}...")
        """
        Performs a fixed-domain-only protein best fit following Fortran DynDom approach.
        Fits conformation 2's fixed domain onto conformation 1's fixed domain, moving entire conformation 2.
        :param dynamic_domain: The domain connected to the fixed domain (not used in fitting)
        :param fixed_domain: The fixed domain (used as sole basis for fitting)
        :return r: A SupResult object containing RMSD, Center 1, Center 2, Rotation Matrix, and Translation Vector
        :return slide_window_2: The slide window residue chain of Protein 2 after fitting (transformation)
                                to Protein 1's fixed domain position
        """
        slide_window_1 = self.protein_1.get_slide_window_residues()
        slide_window_2 = self.protein_2.get_slide_window_residues()
        # Lists to store the backbone atom coordinates of ONLY the fixed domain residues in Protein 1 and 2.
        coords_1 = []
        coords_2 = []
        
        # Extract coordinates from ONLY the fixed domain segments (dynamic domain is ignored for fitting)
        for s in fixed_domain.segments:
            for i in range(s[0], s[1] + 1):
                for a in self.atoms_to_use:
                    coords_1.append(slide_window_1[i][a][0].pos)
                    coords_2.append(slide_window_2[i][a][0].pos)
        
        # The superposition result of Protein 2's fixed domain fitting onto Protein 1's fixed domain
        # This transformation will be applied to the ENTIRE Protein 2 slide window
        r: gemmi.SupResult = gemmi.superpose_positions(coords_1, coords_2)
        slide_window_2.transform_pos_and_adp(r.transform)
        
        return r, slide_window_2

    def check_ratios(self, transformed_protein, dynamic_domain: Domain, fixed_domain: Domain):
        """
        Checks the ratio of external to internal domain movements of each fixed-connected domain pairs.
        :param transformed_protein: The mass-weight-fitted Protein slide window chain .
        :param dynamic_domain: The domain connected to the fixed domain
        :param fixed_domain: The fixed domain
        :return:
        """
        # Superimpose the fixed domain of the original protein onto the mass-weighted fitted protein
        fixed_domain_r: gemmi.SupResult = self.superimpose_domain(transformed_protein, fixed_domain)
        fixed_domain_num_atoms = fixed_domain.num_residues * len(self.atoms_to_use)

        # Calculate the internal and external motions
        fixed_domain_int_msf = self.calc_domain_int_msf(fixed_domain_num_atoms, fixed_domain_r)
        fixed_domain_ext_msf = self.calc_domain_ext_msf(fixed_domain, fixed_domain_r)

        # Same for the dynamic domain
        dynamic_domain_r: gemmi.SupResult = self.superimpose_domain(transformed_protein, dynamic_domain)
        domain_num_atoms = dynamic_domain.num_residues * len(self.atoms_to_use)
        connected_domain_int_msf = self.calc_domain_int_msf(domain_num_atoms, dynamic_domain_r)
        connected_domain_ext_msf = self.calc_domain_ext_msf(dynamic_domain, dynamic_domain_r)

        ratio = self.calc_ext_int_ratio(fixed_domain_ext_msf, fixed_domain_int_msf, fixed_domain_num_atoms,
                                        connected_domain_ext_msf, connected_domain_int_msf, domain_num_atoms)
        if ratio < self.ratio:
            print("Ratio below minimum criteria. Break.")
            return False
        else:
            dynamic_domain.ratio = ratio
            return True

    def superimpose_domain(self, residue_span, domain: Domain, protein=1):
        """
        Superimposes each Protein domain's original residue atoms onto the transformed residue atoms
        :param residue_span: A ResidueSpan object
        :param domain: Domain object
        :param protein: The ID of the Protein to be used to superimpose onto residue_span
        :return results: A list of gemmi.SupResult objects
        """
        # Get the original chain from the Protein that will be used to superimpose onto the target
        protein: Protein = self.protein_1 if protein == 1 else self.protein_2
        fitting_chain: gemmi.ResidueSpan = protein.get_slide_window_residues()
        # Original residue atoms as the fitting
        fitting_coords = []
        # The transformed residue atoms as the target
        target_coords = []
        for s in domain.segments:
            for si in range(s[0], s[1] + 1):
                for a in self.atoms_to_use:
                    fitting_coords.append(fitting_chain[si][a][0].pos)
                    target_coords.append(residue_span[si][a][0].pos)
        # Superpose the domain residue atoms onto the target residue atoms
        r = gemmi.superpose_positions(target_coords, fitting_coords)
        return r

    def calc_domain_int_msf(self, domain_atoms: int, r: gemmi.SupResult):
        """
        Calculates the domain's internal Mean Square Fluctuation
        :param domain_atoms: The specified domain's number of atoms
        :param r: The gemmi.SupResult associated with the domain
        :return:
        """
        return (r.rmsd ** 2) * domain_atoms

    def calc_domain_ext_msf(self, domain, r, protein_id=1):
        """
        Calculates the domain's external Mean Square Fluctuation. The function first transforms the domain chain using
        the given r, then calculates the displacement vectors between atoms of the original domain chain and the
        transformed domain chain.
        :param domain: A Domain object
        :param r: The gemmi.SupResult associated with the domain
        :param protein_id: The ID of the Protein used
        :return:
        """
        protein: Protein = self.protein_1 if protein_id == 1 else self.protein_2
        slide_residues = protein.get_slide_window_residues()
        original_chain: gemmi.Chain = gemmi.Chain(protein.chain_param)
        transformed_chain: gemmi.Chain = gemmi.Chain(protein.chain_param)
        for s in domain.segments:
            for i in range(s[0], s[1]+1):
                original_chain.add_residue(slide_residues[i])
                transformed_chain.add_residue(slide_residues[i])
        transformed_polymer: gemmi.ResidueSpan = transformed_chain.get_polymer()
        transformed_polymer.transform_pos_and_adp(r.transform)
        ext_msf = 0
        for r in range(len(original_chain)):
            for a in self.atoms_to_use:
                atom_coords = np.asarray(original_chain[r][a][0].pos.tolist())
                transformed_atom_coords = np.asarray(transformed_chain[r][a][0].pos.tolist())
                disp_vec = atom_coords - transformed_atom_coords
                sum_disp = np.sum(disp_vec ** 2)
                ext_msf += sum_disp

        return ext_msf

    def calc_ext_int_ratio(self, fixed_ext, fixed_int, fixed_num_residues,
                           dynamic_ext, dynamic_int, dynamic_num_residues):
        """
        Calculate the ratio of external to internal domain motions
        :param fixed_ext: External motion of the fixed domain
        :param fixed_int: Internal motion of the fixed domain
        :param fixed_num_residues: Number of residues in the fixed domain
        :param dynamic_ext: External motion of the dynamic domain
        :param dynamic_int: Internal motion of the dynamic domain
        :param dynamic_num_residues: Number of residues in the dynamic domain
        :return:
        """
        sum_exts = (fixed_ext/fixed_num_residues) + (dynamic_ext / dynamic_num_residues)
        sum_ints = (fixed_int/fixed_num_residues) + (dynamic_int / dynamic_num_residues)
        ratio = math.sqrt(sum_exts/sum_ints)
        return ratio

    def check_ratios_hierarchical(self, temp_domains):
        """
        Check ratios using hierarchical domain reference system
        Replace the existing single-fixed-domain ratio checking
        """
        # Create hierarchical system for this set of domains
        hierarchy_system = HierarchicalDomainSystem(temp_domains)
        analysis_pairs = hierarchy_system.get_analysis_pairs()
        
        if not analysis_pairs:
            print("No analysis pairs found - single domain or no connections")
            return True  # Single domain case
        
        all_ratios_met = True
        ratios = []
        
        for moving_domain_id, reference_domain_id in analysis_pairs:
            moving_domain = temp_domains[moving_domain_id]
            reference_domain = temp_domains[reference_domain_id]
            
            print(f"Analyzing Domain {moving_domain_id} relative to Domain {reference_domain_id}")
            
            # Perform mass-weighted fit between specific domain pair
            r, transformed_slide_window = self.mass_weighted_fit(moving_domain, reference_domain)
            
            # Check ratios for this specific pair
            ratio_met = self.check_ratios(transformed_slide_window, moving_domain, reference_domain)
            
            if not ratio_met:
                print(f"Ratio not met for Domain {moving_domain_id} relative to Domain {reference_domain_id}")
                all_ratios_met = False
                break
            else:
                print(f"Ratio satisfied: {moving_domain.ratio:.3f}")
                ratios.append(moving_domain.ratio)
        
        # Store the hierarchy system for later use
        if all_ratios_met:
            self.hierarchy_system = hierarchy_system
            self.analysis_pairs = analysis_pairs
        
        return all_ratios_met

    def get_hierarchical_fixed_domain(self):
        """
        Get the global reference domain from hierarchy (first in hierarchy)
        This replaces the single fixed domain for global operations
        """
        if hasattr(self, 'hierarchy_system') and self.hierarchy_system:
            return self.hierarchy_system.domain_hierarchy[0]
        else:
            # Fallback to original method
            return self.find_fixed_domain(self.domains)

