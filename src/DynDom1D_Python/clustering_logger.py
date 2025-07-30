import os
from datetime import datetime


class ClusteringLogger:
    """
    Simple logger that writes clustering progress to a text file.
    """
    
    def __init__(self, output_dir, protein_name, run_id=None):
        """Initialize simple text logger."""
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.protein_name = protein_name
        self.run_id = run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Simple text file
        self.log_file = os.path.join(output_dir, f"{protein_name}_{self.run_id}_clustering.log")
        
        # Initialize log file
        with open(self.log_file, 'w', encoding='utf-8') as f:
            f.write(f"DynDom Clustering Log - {protein_name}\n")
            f.write(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 60 + "\n\n")
        
        print(f"Clustering log: {self.log_file}")
    
    def log_k_attempt(self, k, window_size, domains, cluster_segments, 
                     cluster_residues_small, domain_build_failed, 
                     ratios=None, ratio_decision=None, decision=None, reason=None):
        """Log a k-value attempt with simple text output."""
        
        with open(self.log_file, 'a', encoding='utf-8') as f:
            f.write(f"Window {window_size}, K={k} - {datetime.now().strftime('%H:%M:%S')}\n")
            
            # Basic cluster info
            if cluster_segments:
                num_clusters = len(cluster_segments)
                cluster_sizes = []
                for segments in cluster_segments.values():
                    size = sum(seg[1] + 1 - seg[0] for seg in segments)
                    cluster_sizes.append(size)
                
                f.write(f"  Clusters: {num_clusters} (sizes: {cluster_sizes})\n")
                f.write(f"  Min cluster size: {min(cluster_sizes) if cluster_sizes else 0}\n")
            
            # Check results
            if cluster_residues_small:
                f.write(f"  FAIL: Cluster too small (< 20 residues)\n")
            elif domain_build_failed:
                f.write(f"  FAIL: Domain building failed\n")
            elif ratios and ratio_decision == 'fail':
                avg_ratio = sum(ratios) / len(ratios)
                f.write(f"  FAIL: Ratio check failed (avg: {avg_ratio:.3f})\n")
            elif decision == 'accept':
                if ratios:
                    avg_ratio = sum(ratios) / len(ratios)
                    f.write(f"  PASS: All checks passed (avg ratio: {avg_ratio:.3f})\n")
                else:
                    f.write(f"  PASS: All checks passed\n")
                
                if domains:
                    domain_sizes = [d.num_residues for d in domains]
                    f.write(f"  Domains: {len(domains)} (sizes: {domain_sizes})\n")
            
            if reason:
                f.write(f"  Reason: {reason}\n")
            
            f.write("\n")
    
    def log_window_change(self, old_window, new_window, reason):
        """Log window size change."""
        with open(self.log_file, 'a', encoding='utf-8') as f:
            f.write(f"WINDOW SIZE CHANGE: {old_window} -> {new_window}\n")
            f.write(f"  Reason: {reason}\n")
            f.write(f"  Time: {datetime.now().strftime('%H:%M:%S')}\n")
            f.write("-" * 40 + "\n\n")
    
    def log_final_result(self, final_k, final_domains, final_window, 
                        clusterer_status, total_attempts):
        """Log final result."""
        with open(self.log_file, 'a', encoding='utf-8') as f:
            f.write("=" * 60 + "\n")
            f.write("FINAL RESULT\n")
            f.write(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Final K: {final_k}\n")
            f.write(f"Final window: {final_window}\n")
            f.write(f"Status: {clusterer_status}\n")
            f.write(f"Number of domains: {len(final_domains) if final_domains else 0}\n")
            
            if final_domains:
                domain_sizes = [d.num_residues for d in final_domains]
                f.write(f"Domain sizes: {domain_sizes}\n")
            
            f.write("=" * 60 + "\n")
        
        # print(f"Clustering complete - see {self.log_file}")