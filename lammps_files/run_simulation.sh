#!/bin/bash
#SBATCH --job-name=tetracyclic_md
#SBATCH --output=tetracyclic_%j.out
#SBATCH --error=tetracyclic_%j.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=01:00:00 # Increased time for 2 runs
#SBATCH --partition=phd_student
#SBATCH --qos=phd_student

module load openmpi-4.1.5
module load lammps-openmpi
module load anaconda3

export OMP_NUM_THREADS=1

echo "Starting Group 6 MD Simulation Loop"
echo "Date: $(date)"
echo "=========================================="

# --- MODIFIED SECTION ---

# 1. ADD YOUR TWO DATA FILES TO THIS LIST
#    IMPORTANT: Put the alpha file FIRST and the tree file SECOND.
DATA_FILES=(
    "alpha_polymer.data"
    "tree_polymer.data"
)

# 2. Define the *exact* output filenames that will be generated
#    These MUST match the data file names from the list above.
GYR_FILE_ALPHA="gyration.alpha_polymer.txt"
GYR_FILE_TREE="gyration.tree_polymer.txt"

# 3. Loop over each file and run LAMMPS
for F in "${DATA_FILES[@]}"; do
    
    # Create a unique base name for all output files
    BASENAME=$(basename "$F" .data)

    echo "--- Starting LAMMPS for: $F ---"
    
    # Define unique output filenames for this run
    LOG_FILE="log.${BASENAME}.lammps"
    GYR_FILE="gyration.${BASENAME}.txt"
    DUMP_RELAX="traj.${BASENAME}.relax.lammpstrj"
    DUMP_PROD="traj.${BASENAME}.prod.lammpstrj"

    # Pass all filenames as variables to the .in script
    mpirun -np 1 lmp -in alpha_polymer.in \
        -var DATAFILE "$F" \
        -var LOGFILE "$LOG_FILE" \
        -var GYRFILE "$GYR_FILE" \
        -var DUMP_RELAX "$DUMP_RELAX" \
        -var DUMP_PROD "$DUMP_PROD"

    if [ $? -eq 0 ]; then
        echo "SUCCESS: LAMMPS simulation completed for $F!"
        echo "=========================================="
    else
        echo "ERROR: LAMMPS simulation FAILED for $F!"
        exit 1
    fi
done

# 4. RUN ANALYSIS *AFTER* BOTH SIMULATIONS ARE DONE
echo "--- All simulations complete. Running g-factor analysis... ---"

# Pass BOTH gyration files to the new python script
python3 compute_gyration.py "$GYR_FILE_ALPHA" "$GYR_FILE_TREE"

if [ $? -eq 0 ]; then
    echo "Analysis completed successfully!"
    echo "Check g_factor_results.txt and g_factor_analysis_results.png."
else
    echo "ERROR: Python analysis script failed!"
    exit 1
fi

echo "All jobs completed at: $(date)"
