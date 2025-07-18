# DynDom-Py Methodology

This document provides detailed algorithmic descriptions and mathematical foundations for the DynDom-Py implementation.

## Table of Contents

1. [Theoretical Foundation](#theoretical-foundation)
2. [Algorithm Details](#algorithm-details)
3. [Mathematical Formulations](#mathematical-formulations)
4. [Validation Criteria](#validation-criteria)
5. [Implementation Notes](#implementation-notes)

## Theoretical Foundation

### Chasles' Theorem and Screw Motion

The mathematical foundation of DynDom is based on **Chasles' theorem**, which states that the most general displacement of a rigid body can be described as a screw motion - a combination of rotation about an axis and translation along that same axis.

For any rigid body transformation, the motion can be decomposed into:
- **Rotation angle** (θ): Amount of rotation about the screw axis
- **Translation distance** (d): Movement along the screw axis
- **Screw axis**: A line in 3D space defined by direction and location

### Domain Motion Paradigm

DynDom operates under the assumption that protein conformational changes can be meaningfully described as the motion of quasi-rigid domains connected by flexible hinges. This paradigm is valid when:

1. **Interdomain deformation ≥ Intradomain deformation**
2. **Domains maintain internal rigidity** during the conformational change
3. **Connecting regions** allow sufficient flexibility for domain movement

## Algorithm Details

### Step 1: Global Structure Alignment

**Purpose**: Establish a common reference frame for comparing the two conformational states.

**Method**: 
- Performs least-squares superposition using all backbone atoms (N, Cα, C)
- Uses GEMMI library's `calculate_superposition()` with MainChain selection
- Applies transformation to structure 2, moving it into structure 1's coordinate system

**Output**: 
- Global RMSD value
- Transformation matrix for structure 2

### Step 2: Local Motion Detection via Sliding Window

**Purpose**: Detect local conformational differences along the protein backbone.

**Parameters**:
- Window size (w): Default 5 residues
- Overlap: Windows overlap by (w-1) residues
- Atoms used: Backbone atoms (N, Cα, C) or Cα only

**Process**:
For each window position i (from residue w/2 to n-w/2):
1. Extract backbone atoms from residues [i-w/2, i+w/2] in both structures
2. Perform local least-squares superposition
3. Extract rotation matrix R_i and translation vector t_i
4. Assign motion parameters to central residue i

**Mathematical Details**:
The rotation matrix R_i describes the optimal rotation to superimpose window i from structure 1 onto structure 2. This is obtained by solving:

```
min Σ ||R_i * x_1j + t_i - x_2j||²
```

where x_1j and x_2j are corresponding atom positions in the window.

### Step 3: Rotation Vector Conversion

**Purpose**: Convert rotation matrices to rotation vectors for clustering analysis.

**Method**: 
- Uses Rodrigues' rotation formula via scipy's `Rotation.from_matrix()`
- Converts to rotation vectors in degrees: **r_i** = θ **û**
  - θ: rotation angle
  - **û**: unit vector along rotation axis

**Properties**:
- Rotation vector magnitude = rotation angle
- Rotation vector direction = rotation axis
- Enables Euclidean clustering in 3D rotation space

### Step 4: Iterative K-means Clustering with Adaptive Window Sizing

**Purpose**: Group residues with similar local rotational behavior through systematic exploration of cluster numbers and automatic parameter optimization.

**Overall Strategy**: 
- For current window size, start with K=2 clusters and incrementally increase K
- For each K value, attempt clustering and validate results
- Continue until valid domains found or maximum K reached
- If no valid result found across all K values, increase window size by 2 and restart

**K-means Algorithm**: 
- Initialize centroids using custom variance-based splitting
- Apply standard K-means with configurable iterations and restarts
- Evaluate cluster quality and domain validity for each K

**Centroid Initialization**:
1. Calculate variance along each dimension (x, y, z) for all rotation vectors
2. Identify cluster with maximum variance
3. Split this cluster along its highest-variance dimension
4. Replace original centroid with two new centroids above/below the mean

**K-means Convergence**: For each individual K-means run at a given K value:
- Maximum iterations reached (default: 500), OR
- Centroids change by less than tolerance, OR
- No improvement in within-cluster sum of squares

**Domain Validation**: After each K-means run, check if resulting domains meet criteria:
- All domains ≥ 20 residues, AND  
- All domain pairs have motion ratio ≥ 1.0

**Failure Counting**: If validation fails, increment failure counter and try K+1. After 5 consecutive validation failures, increase window size.

**Adaptive Window Size Optimization**:

**When Applied**: When 5 consecutive clustering attempts fail validation criteria, regardless of whether K_max has been reached.

**Failure Tracking**: Monitor failures across the K-means iteration process:
1. **Size Failure**: Any domain contains fewer than 20 residues
2. **Ratio Failure**: Any domain pair has interdomain/intradomain motion ratio < 1.0

**Adaptive Response**:
- Track consecutive failures across all K values at current window size
- After **5 consecutive failed attempts** (across K iterations), increase window by **2**
- Reset to K=2 and restart the entire clustering process with new window size
- Continue until valid domains found or maximum window size reached

**Window Size Selection**:
- **Starting Point**: Window = 5 residues
  - Provides good balance between local resolution and noise reduction
  - Captures local backbone conformational changes effectively
  - Small enough to detect fine-grained domain boundaries
- **Typical Range**: Most proteins succeed with windows 5-7
  - Window 7 often resolves cases where window 5 creates overly fragmented domains
  - Window sizes >7 are rare and typically indicate very noisy conformational changes
- **Upper Practical Limit**: Rarely exceeds window = 9-11
  - Beyond this, risk of merging distinct functional domains
  - May lose biologically relevant domain boundaries

**Mechanistic Effects of Larger Windows**:
- **Noise Smoothing**: Average out local backbone fluctuations and crystal packing effects
- **Coherent Motion Detection**: Better identification of truly rigid regions
- **Reduced Fragmentation**: Fewer artificial domain breaks due to local noise
- **Improved Connectivity**: Spatially adjacent regions more likely to cluster together
- **Enhanced Signal-to-Noise**: Interdomain motion signal preserved while intradomain noise averaged out

**Process Flow**:
```
Window 5: Try K=2,3,4... → Valid result found → SUCCESS
Window 5: Try K=2,3,4... → All K failed → Count failure
[Repeat until 5 failures]
Window 7: Try K=2,3,4... → Valid result found → SUCCESS
Window 7: Try K=2,3,4... → All K failed → Count failure
[Continue as needed]
```

**Important**: If a valid solution is found at any K value, the algorithm terminates successfully. Window size adaptation only occurs when **no valid result exists** after trying all reasonable K values.

### Step 5: Dynamic Domain Construction

**Purpose**: Build spatially connected domains from clustered rotation vectors.

#### 5.1 Binary Connectivity Matrix

For each cluster, create an N×N binary matrix where N = number of segments in cluster:
- Matrix[i,j] = 1 if segments i and j are spatially connected
- Matrix[i,j] = 0 otherwise

**Connectivity Criterion**: Two segments are connected if any atom in segment i is within distance threshold (4.0 Å for side chains, 10.0 Å for backbone-only) of any atom in segment j.

#### 5.2 Connected Component Analysis

Apply row reduction algorithm to identify connected components:
1. Start with connectivity matrix
2. For each row i, perform OR-wise operations with all connected rows
3. Merge overlapping components
4. Result: List of maximally connected segment groups

**Output**: Each connected component becomes a potential domain.

### Step 6: Domain Validation

#### 6.1 Size Filtering

**Minimum Domain Size**: Default 20 residues
- Domains smaller than threshold are merged with larger neighboring domains
- If all domains in a cluster are too small, this constitutes a **failure condition**

#### 6.2 Motion Ratio Validation

For each pair of connected domains, calculate the ratio:

```
R = √[(Σ E_ext/N_ext) / (Σ E_int/N_int)]
```

Where:
- **E_ext**: External (interdomain) mean square displacement
- **E_int**: Internal (intradomain) mean square displacement  
- **N_ext, N_int**: Number of atoms in external/internal calculations

**Calculation Process**:
1. Perform mass-weighted superposition using both domains
2. Superimpose each domain individually onto the mass-weighted result
3. Calculate displacement vectors for domain atoms
4. Decompose displacements into internal (domain superposition error) and external (relative domain motion) components

**Acceptance Criterion**: R ≥ R_min (default 1.0)
- If any domain pair has R < 1.0, this constitutes a **failure condition**

**Important**: If a valid solution has already been found at the current window size, the algorithm terminates successfully. Window size adaptation only occurs when **no valid result exists yet** after trying all K values.

### Step 7: Screw Axis Calculation

**Purpose**: Determine the screw axis describing motion between domain pairs.

#### Mathematical Framework

For a domain moving from position **P₁** to **P₂**, the screw motion is characterized by:

1. **Rotation Matrix R**: Obtained from domain superposition
2. **Displacement Vector d**: Centroid motion vector
3. **Rotation Axis û**: Eigenvector of R with eigenvalue 1
4. **Rotation Angle θ**: From trace(R) = 1 + 2cos(θ)

#### Screw Axis Location

The screw axis passes through point **p₀** such that:

```
p₀ = p_cent + 0.5 * d_rot - (|d_rot| / (2 * tan(θ/2))) * (d_rot × û) / |d_rot × û|
```

Where:
- **p_cent**: Domain centroid
- **d_rot**: Rotational component of displacement (d - d_parallel)
- **d_parallel**: Component of d parallel to û

#### Translation Component

Translation along the screw axis:
```
T = d · û
```

### Step 8: Hinge Residue Identification

**Purpose**: Locate residues involved in interdomain bending.

#### Statistical Method

1. **Domain Distributions**: Model rotation vectors within each domain as 3D normal distributions
2. **Covariance Matrices**: Calculate covariance matrix Σ for each domain's rotation vectors
3. **Outlier Detection**: Identify residues at domain boundaries with "unusual" rotations

#### Mahalanobis Distance Test

For residue i at a domain boundary:
```
Q_i = (r_i - μ)ᵀ Σ⁻¹ (r_i - μ)
```

Where:
- **r_i**: Rotation vector for residue i  
- **μ**: Mean rotation vector for the domain
- **Σ**: Covariance matrix for the domain

**Classification**: Residue i is a "bending residue" if Q_i > χ²₃(P) where P = 0.2 (80% confidence threshold).

#### Directional Search

Starting from domain boundaries, extend search in both directions until:
- Mahalanobis distance falls below threshold, OR
- Domain boundary reached, OR
- Maximum search distance exceeded

## Validation Criteria

### Geometric Constraints

1. **Minimum Domain Size**: ≥ 20 residues (configurable)
   - **Failure condition** if violated
   - Larger windows help by creating more coherent domains
2. **Spatial Connectivity**: 4.0 Å threshold for atom-atom contacts
3. **Backbone Continuity**: Domains must form connected backbone segments

### Physical Constraints  

1. **Motion Ratio**: R ≥ 1.0 (interdomain ≥ intradomain motion)
   - **Failure condition** if violated  
   - Larger windows improve ratios by averaging out intradomain noise
2. **Rotation Validity**: Rotation matrices must be proper (det(R) = 1)
3. **Displacement Magnitude**: Reasonable displacement vectors (< 100 Å)

### Statistical Validation

1. **Cluster Quality**: Within-cluster sum of squares minimization
2. **Domain Homogeneity**: Low variance within domain rotation vectors
3. **Boundary Significance**: Clear separation between domain rotation distributions

### Adaptive Optimization Strategy

The window size adaptation mechanism addresses common failure modes:

**Why Larger Windows Help**:
- **Noise Reduction**: Smooth out local backbone fluctuations that create artificial domain boundaries
- **Coherent Domains**: Neighboring residues with similar motion get grouped together more effectively  
- **Size Compliance**: Reduced fragmentation leads to naturally larger domains
- **Improved Ratios**: Signal-to-noise improvement where interdomain motion signal is preserved while intradomain noise is averaged out
- **Better Connectivity**: Spatially adjacent regions more likely to have similar motion vectors

**Trade-offs**:
- **Resolution Loss**: Larger windows may miss fine-grained hinge locations
- **Boundary Blurring**: Very large windows might merge distinct but small domains
- **Computational Cost**: More iterations required for convergence

## Implementation Notes

### Numerical Stability

1. **Matrix Decomposition**: Uses SVD for robust superposition calculations
2. **Rotation Conversion**: Handles numerical precision in rotation matrix → vector conversion
3. **Clustering Convergence**: Multiple random initializations to avoid local minima

### Performance Optimizations

1. **Sparse Connectivity**: Only calculate distances for spatially close segments
2. **Incremental Clustering**: Stop when motion ratio criteria cannot be satisfied
3. **Memory Management**: Process large structures in chunks when necessary

### Error Handling

1. **Insufficient Data**: Graceful handling of proteins too small for analysis
2. **Degenerate Cases**: Detection and handling of very small rotations
3. **Convergence Failure**: Fallback strategies for clustering non-convergence

### Coordinate System Considerations

1. **Reference Frame**: Structure 1 defines the reference coordinate system
2. **Transformation Chain**: Careful tracking of coordinate transformations
3. **Axis Orientation**: Consistent orientation of screw axes using right-hand rule

## References

1. Hayward, S. and Berendsen, H.J.C. (1998) Systematic analysis of domain motions in proteins from conformational change. *Proteins* **30**, 144-154.

2. Wriggers, W. and Schulten, K. (1997) Protein domain movements: detection of rigid domains and visualization of hinges in comparisons of atomic coordinates. *Proteins* **29**, 1-14.

3. Gerstein, M. and Krebs, W. (1998) A database of macromolecular motions. *Nucleic Acids Res* **26**, 4280-4290.

4. Flores, S., Echols, N., Milburn, D., Hespenheide, B., Keating, K., Lu, J., Wells, S., Yu, E.Z., Thorpe, M. and Gerstein, M. (2006) The Database of Macromolecular Motions: new features added at the decade mark. *Nucleic Acids Res* **34**, D296-301.
