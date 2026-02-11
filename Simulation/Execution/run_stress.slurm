#!/bin/bash
# Stress Tensor Post-Processing Job (SLURM)
# Author: Noelia Olivera Rodríguez (2026)
# This script performs batch post-processing of ESyS-Particle interaction
# output files (partForce.*.dat) to compute coarse-grained stress fields.
# For each saved timestep, the script:
#   1. Reads particle interaction forces.
#   2. Computes the stress tensor on a regular 3D grid.
#   3. Exports the result in VTK format (.vti).
# The workflow is consistent with the temporal discretization used in the
# simulation (dt, t_snap, total simulation time).

# Slurm resource configuration

#SBATCH --job-name=stress
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4098
#SBATCH --time=4:00:00
#SBATCH --tmp=9G
#SBATCH --partition=normal
#SBATCH --qos=normal

# Temporal parameters (must match simulation)
DT=0.0000003        # DEM time step (s)
TE=1.8              # Total simulation time (s)
TSNAP=0.00375       # Snapshot interval (s)


# Spatial grid parameters for stress coarse-graining
XMIN=0
XMAX=0.4
YMIN=0
YMAX=0.3970
ZMIN=0
ZMAX=0.4

NX=18
NY=18
NZ=18
# Number of grid cells in each spatial direction

# Environment initialization
source /clusteruy/home/gonzalo/ENV/bin/activate
# Activate Python virtual environment containing required packages


# Step size calculation

# Compute number of time steps between snapshots
DSTEP=$(echo "$TSNAP/$DT" | bc)
DSTEP=$(printf "%d" $DSTEP)

# Explicitly enforce expected snapshot stride (safety override)
DSTEP=6250

STEP=0
ts=$DT   # Initialize physical time


# Main loop over saved timesteps
while true; do

    # Check if current simulation time is within total duration
    if [ $(echo "$TE >= $ts" | bc) -eq 1 ]; then

        INPUT_FILE="/clusteruy/home/nolivera/ESyS/CompAndDescCyc/3CiclosBajo/partForce.${STEP}.dat"
        OUTPUT_FILE="/clusteruy/home/nolivera/ESyS/CompAndDescCyc/stress/stress3CiclosBajo/stress.${STEP}"

        # Verify existence of input file
        if [ -f "$INPUT_FILE" ]; then

            echo "Processing $INPUT_FILE"

            python /clusteruy/home/nolivera/ESyS/CompAndDescCyc/stress/stress3vti.py \
                "$INPUT_FILE" "$OUTPUT_FILE" \
                $XMIN $XMAX $YMIN $YMAX $ZMIN $ZMAX \
                $NX $NY $NZ

        else
            echo "File $INPUT_FILE not found. Stopping."
            break
        fi

        # Advance physical time and timestep index
        ts=$(echo "$ts + $DT" | bc)
        STEP=$(echo "$STEP + $DSTEP" | bc)

        echo "Snapshot stride: $DSTEP"
        echo "Current STEP: $STEP"

    else
        echo "Final time exceeded (ts > TE). Exiting loop."
        break
    fi

done

echo "DONE"

deactivate
