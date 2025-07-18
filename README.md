# DynDom-Py

A Python implementation of the DynDom algorithm for analyzing protein domain movements between conformational states.

## Overview

DynDom-Py is a Python implementation of the DynDom algorithm for analyzing protein conformational changes in terms of rigid-body domain movements. The software determines protein domains, hinge axes, and amino acid residues involved in hinge bending through fully automated analysis.

Given two conformational states of the same protein (e.g., from X-ray crystallography, NMR, or molecular dynamics simulations), DynDom-Py identifies how the protein moves by treating different regions as quasi-rigid bodies. This approach transforms complex conformational changes into an easily understood view of domain movements connected by flexible hinges.

The analysis reveals:
- **Dynamic domains**: Regions that move as rigid bodies
- **Hinge regions**: Flexible areas that allow domain movement  
- **Screw axes**: Mathematical description of how domains rotate and translate relative to each other

## Features

- **Automated Domain Detection**: Identifies rigid domains using rotation vector clustering
- **Hierarchical Analysis**: Uses advanced hierarchical domain reference system for complex multi-domain proteins
- **Hinge Analysis**: Locates flexible regions and mechanical hinges between domains
- **Screw Axis Calculation**: Determines interdomain rotation axes, angles, and translations
- **Motion Validation**: Applies geometric constraints and interdomain/intradomain motion ratios
- **Comprehensive Visualization**: Generates PyMOL scripts with domain coloring and 3D motion arrows
- **Flexible Input**: Supports any two conformational states in PDB format
- **Detailed Output**: Provides motion statistics, hinge residue identification, and structural files
- **Quality Control**: Validates domain connectivity and motion significance
- **Clustering Logs**: Detailed logging of the clustering process for debugging and analysis

## Installation

### Requirements

- Python 3.7 or higher
- Required packages listed in `requirements.txt`

### Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/dyndom-py.git
   cd dyndom-py
   ```

2. **Create a virtual environment** (recommended):
   ```bash
   python -m venv dyndom-env
   source dyndom-env/bin/activate  # On Windows: dyndom-env\Scripts\activate
   ```

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

   Or install manually:
   ```bash
   pip install numpy scipy scikit-learn gemmi
   ```

4. **Verify installation**:
   ```bash
   python -c "import numpy, scipy, sklearn, gemmi; print('All dependencies installed successfully')"
   ```

5. **Set up directories**:
   ```bash
   mkdir -p data output
   ```

## Usage

### Basic Usage

1. **Prepare input files** in the `data/` directory:

   **data/command.txt**:
   ```
   input_path=data
   output_path=output
   filename1=1hng
   chain1id=A
   filename2=2hng
   chain2id=A
   ```

   **data/param.txt**:
   ```
   window=5
   domain=20
   ratio=1.0
   k_means_n_init=1
   k_means_max_iter=500
   atoms=backbone
   ```

2. **Prepare PDB files**:
   - Place your PDB files in the `data/` directory, or
   - Use PDB codes (e.g., `1hng`, `2hng`) and the software will automatically download them from RCSB

3. **Run the analysis**:
   ```bash
   python main.py
   ```

### Parameters

#### Command Parameters (`data/command.txt`)
- `input_path`: Directory containing PDB files
- `output_path`: Directory for output files
- `filename1`, `filename2`: PDB IDs or filenames (without .pdb extension)
- `chain1id`, `chain2id`: Chain identifiers

#### Analysis Parameters (`data/param.txt`)
- `window`: Sliding window size for local motion analysis (default: 5)
- `domain`: Minimum domain size in residues (default: 20)
- `ratio`: Minimum ratio of interdomain to intradomain motion (default: 1.0)
- `atoms`: Atoms to use for analysis (`backbone` for N,CA,C or `ca` for CA only)
- `k_means_n_init`: Number of K-means initializations (default: 1)
- `k_means_max_iter`: Maximum K-means iterations (default: 500)

## Output Files

The software generates several output files in the specified output directory:

### Structure Files
- **`{protein1}_{chain1}_{protein2}_{chain2}.pdb`**: Superimposed structures (two models)
- **`{protein1}_{chain1}_{protein2}_{chain2}_arrows.pdb`**: 3D motion arrows for PyMOL visualization
- **`{protein1}_rot_vecs.pdb`**: Rotation vectors for each residue

### Visualization Files
- **`{protein1}_{chain1}_{protein2}_{chain2}.pml`**: PyMOL visualization script with:
  - Domain coloring (different colors for each domain)
  - Hinge regions highlighted in green
  - 3D motion arrows showing screw axes
  - Automatic camera positioning and lighting

### Analysis Reports
- **`{protein1}_{chain1}_{protein2}_{chain2}.w5_info`**: Detailed analysis results including:
  - Domain definitions and residue ranges
  - Rotation angles and screw axis parameters
  - Motion statistics and RMSD values
  - Hierarchical domain relationships

### Debug Files
- **`clustering_logs/`**: Directory containing detailed clustering logs
- **`domain_*_comparison.pdb`**: Individual domain motion comparison files (debug mode)

### Visualization

1. **Load in PyMOL**:
   ```bash
   pymol output/{protein1}_{chain1}_{protein2}_{chain2}.pml
   ```

2. **Visualization features**:
   - **Domain coloring**: Different domains shown in different colors
   - **Global reference domain**: Shown in blue
   - **Moving domains**: Shown in red, yellow, pink, etc.
   - **Hinge regions**: Flexible regions highlighted in green
   - **Motion arrows**: 3D arrows showing screw axes and rotation directions
   - **Arrow coloring**: Shaft shows reference domain color, tip shows moving domain color

## Algorithm Overview

DynDom-Py implements a 7-step automated workflow with hierarchical domain analysis:

1. **Global Structure Alignment**: Performs whole-protein best-fit superposition
2. **Local Motion Detection**: Analyzes rotation vectors using sliding windows
3. **Rotation Vector Clustering**: Uses adaptive K-means clustering with validation
4. **Hierarchical Domain Construction**: Builds domains with connectivity analysis
5. **Domain Validation**: Applies size and motion ratio criteria
6. **Hierarchical Screw Axis Calculation**: Determines motion parameters for each domain pair
7. **Hinge Residue Identification**: Locates flexible regions using statistical analysis

### Key Concepts

- **Hierarchical Analysis**: Domains are analyzed relative to their most appropriate reference domain
- **Global Reference Domain**: The most connected domain serves as the global reference
- **Analysis Pairs**: Each domain is analyzed relative to its optimal reference domain
- **Screw Axis**: Mathematical description of rigid body motion (rotation + translation)
- **Bending Residues**: Residues with rotations outside the main domain distribution

## Troubleshooting

### Common Issues

1. **"Sequence Identity less than 40%"**:
   - Check that both PDB files contain the same protein
   - Verify chain IDs are correct

2. **"Too many fails. Increasing window size"**:
   - Normal behavior - algorithm automatically adjusts parameters

3. **No domains found**:
   - Reduce `domain` parameter (try 15 or 10)
   - Reduce `ratio` parameter (try 0.8 or 0.6)
   - Increase `window` parameter


Original DynDom algorithm:
```
Hayward, S. and Berendsen, H.J.C. (1998) Systematic analysis of domain motions in proteins from conformational change: new results on citrate synthase and T4 lysozyme. Proteins 30, 144-154.
```
```
S. Hayward, R. A. Lee
"Improvements in the analysis of domain motions in proteins from conformational change: DynDom version 1.50" J Mol Graph Model, Dec, 21(3), 181-3, 2002. 
```

## License

BSD 3-Clause License 
Refer to LICENSE file

