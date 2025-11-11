#!/usr/bin/env python3
"""
Tetracyclic vs. Tree G-Factor Analysis
Calculates a custom ratio = <Rg_alpha>^2 / <Rg_tree>^2

Usage: python3 compute_gyration.py <gyration_file_alpha> <gyration_file_tree>
       (This script is called automatically by run_simulation.sh)
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def get_avg_rg_from_file(filename):
    """
    Helper function to read a LAMMPS gyration file, calculate the
    equilibrated average Rg and Rg^2, and return all data for plotting.
    """
    if not os.path.exists(filename):
        print(f"❌ ERROR: {filename} not found!")
        return None

    try:
        # Read gyration data, skip LAMMPS header lines (lines starting with #)
        data = np.loadtxt(filename, comments='#')
        
        # Check if file is empty or has wrong format
        if data.ndim == 0 or data.shape[1] < 2:
             print(f"❌ ERROR: {filename} contains no valid data.")
             return None

        timesteps = data[:, 0]
        rg_values = data[:, 1]

        # Calculate equilibrated average (last 25% of simulation)
        equilibrated_start = (len(rg_values) * 3)// 4
        
        # Handle case with very few data points
        if equilibrated_start == 0 and len(rg_values) > 1:
            equilibrated_start = 1
        elif len(rg_values) <= 1:
            print(f"❌ ERROR: Not enough data in {filename} to analyze.")
            return None

        equilibrated_rg = rg_values[equilibrated_start:]
        
        avg_rg = np.mean(equilibrated_rg)
        std_rg = np.std(equilibrated_rg)
        
        # --- CRITICAL FIX ---
        # The g-factor uses the mean-square Rg, <Rg^2>, 
        # NOT the square-of-the-mean, <Rg>^2.
        avg_rg_squared = np.mean(equilibrated_rg ** 2)
        # --- END FIX ---
        
        print(f"--- Analysis for: {filename} ---")
        print(f"  Avg. Rg: {avg_rg:.3f} +/- {std_rg:.3f}")
        print(f"  Avg. <Rg^2>: {avg_rg_squared:.3f}")
        print(f"  Equilibrated data points: {len(equilibrated_rg)}")
        
        # Return all components for plotting and calculation
        return (avg_rg_squared, avg_rg, timesteps, rg_values, 
                equilibrated_start, equilibrated_rg)

    except Exception as e:
        print(f"❌ ERROR: An error occurred processing {filename}: {str(e)}")
        return None

def main(file_alpha, file_tree):
    """
    Main analysis function to compare two simulation outputs.
    """
    print("=" * 60)
    print("GROUP 6: ALPHA vs. TREE CUSTOM RATIO ANALYSIS")
    print(f"  Alpha (Tetracyclic) file: {file_alpha}")
    print(f"  Tree (Dendrimer) file:    {file_tree}")
    print("=" * 60)
    
    # --- Get data for both files ---
    alpha_data = get_avg_rg_from_file(file_alpha)
    tree_data = get_avg_rg_from_file(file_tree)
    
    if alpha_data is None or tree_data is None:
        print("\n❌ ERROR: One or both files failed to load. Exiting analysis.")
        sys.exit(1) # Exit with error code

    # Unpack the data
    (rg2_alpha, rg_alpha, ts_alpha, 
     val_alpha, start_alpha, hist_alpha) = alpha_data
    (rg2_tree,  rg_tree,  ts_tree,  
     val_tree,  start_tree,  hist_tree) = tree_data

    # --- Calculate the Custom Ratio ---
    # Ratio = <Rg^2>(alpha) / <Rg^2>(tree)
    custom_ratio = rg2_alpha / rg2_tree

    print("\n" + "-" * 60)
    print("CUSTOM RATIO CALCULATION:")
    print(f"  <Rg^2>_alpha: {rg2_alpha:.3f} (Tetracyclic)")
    print(f"  <Rg^2>_tree:  {rg2_tree:.3f} (Dendrimer)")
    print(f"  Custom Ratio = <Rg^2>_alpha / <Rg^2>_tree = {custom_ratio:.4f}")
    print("-" * 60)
    
    # --- Plotting and Saving ---
    fig, axs = plt.subplots(2, 2, figsize=(20, 14), dpi=100) # 2x2 grid
    fig.suptitle(f"Custom Ratio Analysis: Ratio = {custom_ratio:.4f}  (<Rg^2>_alpha / <Rg^2>_tree)", fontsize=20)
    
    # --- Alpha Polymer Plots ---
    axs[0, 0].plot(ts_alpha, val_alpha, 'b-', linewidth=1, alpha=0.7, label='Rg (Alpha)')
    axs[0, 0].axhline(y=rg_alpha, color='r', linestyle='--', linewidth=2, label=f'Avg Rg = {rg_alpha:.3f}')
    axs[0, 0].axvline(x=ts_alpha[start_alpha], color='g', linestyle=':', label='Equilibration')
    axs[0, 0].set_xlabel('Timestep')
    axs[0, 0].set_ylabel('Radius of Gyration (Å)')
    axs[0, 0].set_title(f'Alpha Polymer (Tetracyclic)')
    axs[0, 0].legend()
    axs[0, 0].grid(True, alpha=0.3)
    
    axs[1, 0].hist(hist_alpha, bins=30, alpha=0.7, color='skyblue', edgecolor='black', density=True)
    axs[1, 0].axvline(x=rg_alpha, color='r', linestyle='--', linewidth=2, label=f'Mean = {rg_alpha:.3f}')
    axs[1, 0].set_xlabel('Radius of Gyration (Å)')
    axs[1, 0].set_ylabel('Probability Density')
    axs[1, 0].set_title('Alpha Polymer Rg Distribution')
    axs[1, 0].legend()
    axs[1, 0].grid(True, alpha=0.3)

    # --- Tree Polymer Plots ---
    axs[0, 1].plot(ts_tree, val_tree, 'g-', linewidth=1, alpha=0.7, label='Rg (Tree)')
    axs[0, 1].axhline(y=rg_tree, color='r', linestyle='--', linewidth=2, label=f'Avg Rg = {rg_tree:.3f}')
    axs[0, 1].axvline(x=ts_tree[start_tree], color='g', linestyle=':', label='Equilibration')
    axs[0, 1].set_xlabel('Timestep')
    axs[0, 1].set_ylabel('Radius of Gyration (Å)')
    axs[0, 1].set_title(f'Tree Polymer (Dendrimer)')
    axs[0, 1].legend()
    axs[0, 1].grid(True, alpha=0.3)
    
    axs[1, 1].hist(hist_tree, bins=30, alpha=0.7, color='lime', edgecolor='black', density=True)
    axs[1, 1].axvline(x=rg_tree, color='r', linestyle='--', linewidth=2, label=f'Mean = {rg_tree:.3f}')
    axs[1, 1].set_xlabel('Radius of Gyration (Å)')
    axs[1, 1].set_ylabel('Probability Density')
    axs[1, 1].set_title('Tree Polymer Rg Distribution')
    axs[1, 1].legend()
    axs[1, 1].grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust for suptitle
    plt.savefig('g_factor_analysis_results.png', dpi=300, bbox_inches='tight')
    print(f"Plot saved as: g_factor_analysis_results.png")

    # Save summary results
    with open('g_factor_results.txt', 'w') as f:
        f.write(f"ALPHA vs. TREE CUSTOM RATIO ANALYSIS\n")
        f.write("=" * 50 + "\n")
        f.write(f"Alpha File: {file_alpha}\n")
        f.write(f"Tree File:  {file_tree}\n")
        f.write("-" * 50 + "\n")
        f.write(f"Avg. <Rg^2>_alpha: {rg2_alpha:.6f}\n")
        f.write(f"Avg. <Rg^2>_tree:  {rg2_tree:.6f}\n")
        f.write(f"Custom Ratio:       {custom_ratio:.6f}\n")

    print(f"Results saved as: g_factor_results.txt")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("❌ ERROR: Invalid arguments.")
        print("This script is designed to be run by 'run_simulation.sh'.")
        print("Usage: python3 compute_gyration.py <gyration_file_alpha> <gyration_file_tree>")
        sys.exit(1) # Exit with an error code
    
    file_alpha = sys.argv[1]
    file_tree = sys.argv[2]
    main(file_alpha, file_tree)
    
    print("\n" + "=" * 60)
    print("CUSTOM RATIO ANALYSIS COMPLETE! ✅")
    print("=" * 60)
