# Topological Polymer Analysis - Group 6 Alpha Graph

## Project Overview

This project analyzes the **alpha graph topological polymer** using both molecular dynamics simulations (LAMMPS) and analytical graph theory calculations, based on the methodology presented in Cantarella et al. (2022).

## Objectives

1. **MD Simulation**: Use LAMMPS to simulate the alpha graph polymer and compute radius of gyration
2. **Analytical Calculation**: Apply graph theory to calculate theoretical radius of gyration using Kirchhoff index
3. **Comparison**: Compare simulation results with analytical predictions
4. **Validation**: Verify results against Figure 4 values from the research paper

## Alpha Graph Structure

The alpha graph (Group 6) is one of six topological polymer architectures studied in the paper. It represents a specific graph topology with:
- Multiple vertices (junction points)
- Interconnecting edges (polymer chains)
- Unique topological constraints

## Methodology

### 1. Molecular Dynamics Simulation

**LAMMPS Setup:**
- **Force Field**: Kremer-Grest model with FENE bonds
- **Temperature**: 300 K (NPT ensemble)
- **Duration**: 1,000,000 time steps
- **Analysis**: Radius of gyration computed using `compute gyration` command

**Key Parameters:**
```
pair_style lj/cut 14.0
bond_style fene
temperature 300.0 K
simulation_time 1ns
```

### 2. Graph Theory Analysis

**Kirchhoff Index Calculation:**
- Compute graph Laplacian matrix L
- Calculate Moore-Penrose pseudoinverse L⁺
- Determine resistance distances rᵢⱼ
- Sum over all vertex pairs: Kf(G) = Σᵢ<ⱼ rᵢⱼ

**Expected Radius of Gyration:**
Using Theorem 9 from the paper:
```
E(R²ᵍ; G) = (d/v²) × Kf(G)
```
where:
- d = spatial dimension (3)
- v = number of vertices
- Kf(G) = Kirchhoff index

**Contraction Factor:**
Using Theorem 5 for the asymptotic limit:
```
g(G∞) = (3/e²) × [Tr L⁺(G) + (1/3)Loops(G) - 1/6]
```

## File Structure

```
/group6_project/
├── lammps_files/
│   ├── alpha_polymer.data          # Polymer structure
│   ├── alpha_polymer.in            # LAMMPS input script
│   ├── run_simulation.sh           # Cluster submission script
│   └── analysis_scripts/
│       ├── compute_gyration.py     # MD analysis
│       └── plot_results.py         # Visualization
├── analytical_calculation/
│   ├── kirchhoff_calculation.py    # Graph theory calculations
│   └── comparison_analysis.py      # Compare MD vs theory
├── documentation/
│   ├── README.md                   # This file
│   ├── methodology.md              # Detailed methods
│   └── results_summary.md          # Final results
└── presentation/
    ├── slides.pdf                  # Presentation materials
    └── figures/                    # Generated plots
```

## Expected Results

Based on the research paper, we expect:

1. **Radius of Gyration**: Specific value from Figure 4 table (column 2)
2. **Contraction Factor**: Theoretical g-factor (column 3)
3. **Good Agreement**: MD simulation should match analytical predictions within ~5%

## Running the Analysis

### 1. LAMMPS Simulation
```bash
# On cluster
sbatch run_simulation.sh

# Local execution
lmp_serial -in alpha_polymer.in
```

### 2. Analysis Scripts
```bash
# Process MD results
python3 compute_gyration.py

# Analytical calculations
python3 kirchhoff_calculation.py

# Compare results
python3 comparison_analysis.py
```

## Key Equations

**Radius of Gyration:**
```
Rg² = (1/2v²) Σᵢ,ⱼ ||xᵢ - xⱼ||²
```

**Kirchhoff Index:**
```
Kf(G) = Σᵢ<ⱼ rᵢⱼ = (v/2) Tr(L⁺)
```

**Resistance Distance:**
```
rᵢⱼ = L⁺ᵢᵢ + L⁺ⱼⱼ - L⁺ᵢⱼ - L⁺ⱼᵢ
```

## References

1. Cantarella, J., Deguchi, T., Shonkwiler, C., & Uehara, E. (2022). Radius of gyration, contraction factors, and subdivisions of topological polymers. *Journal of Physics A: Mathematical and Theoretical*, 55(47), 475202.

2. Plimpton, S. (1995). Fast parallel algorithms for short-range molecular dynamics. *Journal of Computational Physics*, 117(1), 1-19.

## Project Timeline

- **Setup Phase**: Create LAMMPS input files and structure
- **Simulation Phase**: Run MD simulations (24-48 hours)
- **Analysis Phase**: Process results and perform calculations
- **Comparison Phase**: Validate against theoretical predictions
- **Documentation Phase**: Prepare final report and presentation

## Contact

Group 6 - Topological Polymer Analysis Project  
Chemistry/Engineering Department  
Based on instructions from Professor  
Submission Date: November 11, 2025