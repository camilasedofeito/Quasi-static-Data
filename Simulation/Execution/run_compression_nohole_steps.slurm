# SLURM submission script for parallel bonded block compression simulations
# using ESyS-Particle (v2.3.5).
#
# This script:
#   1. Requests computational resources via SLURM.
#   2. Loads the required MPI and ESyS-Particle modules.
#   3. Automatically determines the 3D domain decomposition from the
#      simulation parameter file.
#   4. Computes the total number of MPI processes (master + slaves).
#   5. Generates the wall-status parameter file dynamically to ensure
#      consistency with the simulation geometry.
#   6. Launches the parallel compression simulation.
#
# The workflow guarantees consistency between geometric configuration,
# domain decomposition, and MPI execution.
# Original author: Gonzalo Tancredi
# Modified and adapted by: Noelia Olivera Rodríguez


# SLURM RESOURCE CONFIGURATION
#SBATCH --job-name=run_compression          # Job identifier in SLURM queue
#SBATCH --nodes=1                          # Number of compute nodes
#SBATCH --ntasks-per-node=19               # MPI tasks per node
#SBATCH --mem-per-cpu=4096                 # Memory per CPU (MB)
#SBATCH --time=120:00:00                   # Maximum wall-clock time (HH:MM:SS)
#SBATCH --tmp=9G                           # Temporary disk space
#SBATCH --partition=normal                 # Cluster partition
#SBATCH --qos=normal                       # Quality of service
#SBATCH --exclude=node06,node40,node41,node42,node43,node44,node45,node46,node47
# Exclude specific nodes (e.g., unstable or reserved nodes)

source /etc/profile.d/modules.sh            # Initialize module system
module load mpi/mpich-3.2-x86_64            # Load MPI implementation
module load esys-particle/2.3.5             # Load ESyS-Particle


# Domain decomposition (MPI)
# The number of subdomains in each spatial direction (X, Y, Z)
# is read automatically from the simulation parameter file.
# This ensures that the MPI configuration is consistent with
# the internal ESyS-Particle domain decomposition.

nSlavesCompressionX=$(grep numSubdomainsX scripts/Parameters_blockCompression.py | awk -F\= '{print $2}' | awk '{print $1}')
nSlavesCompressionY=$(grep numSubdomainsY scripts/Parameters_blockCompression.py | awk -F\= '{print $2}' | awk '{print $1}')
nSlavesCompressionZ=$(grep numSubdomainsZ scripts/Parameters_blockCompression.py | awk -F\= '{print $2}' | awk '{print $1}')

# Total number of slave processes (3D Cartesian decomposition)
nSlavesCompression=$((nSlavesCompressionX * nSlavesCompressionY * nSlavesCompressionZ))

# Total number of MPI processes:
nProcsCompression=$((nSlavesCompressionX * nSlavesCompressionY * nSlavesCompressionZ + 1))

# The vertical extent of the simulation domain is extracted from the
# general parameter file to dynamically define wall positions.
# This prevents inconsistencies between geometry and boundary conditions.
heightRoof=$(grep "blockSizeY =" scripts/Parameters_general.py | awk -F= '{print $2}')
floorPosition=0


# Generate (overwrite) wall-status parameter file
echo "floorPosition = $floorPosition" > scripts/Parameters_statusWall.py
echo "heightRoof = $heightRoof" >> scripts/Parameters_statusWall.py
echo Compression
		
# Launch bonded block compression simulation in parallel.
# Input geometry corresponds to a pre-equilibrated tagged configuration.  
mpirun -np $nProcsCompression esysparticle scripts/bondedblockCompression.py data/block_equilibrium_tagged_001.geo
