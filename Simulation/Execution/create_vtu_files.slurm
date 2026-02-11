#!/bin/bash
#SBATCH --job-name=run_compression
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=19
#SBATCH --mem-per-cpu=1024
#SBATCH --time=120:00:00
#SBATCH --tmp=9G
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --exclude=node06,node40,node41,node42,node43,node44,node45,node46,node47

source /etc/profile.d/modules.sh
module load mpi/mpich-3.2-x86_64
module load esys-particle/2.3.5
#dump2vtk -i cptPartParedes2/cpt -o cptPartParedes2/frame -t 0 1800 37040
source /clusteruy/home/gonzalo/ENV/bin/activate
#python scripts/analyze_geo.py data/block_con_caras.geo -1.0
#python scripts/tag_geo_py31.py data/block_equilibrium_0.geo data/block_equilibrium_tagged_0004.geo
source /clusteruy/home/gonzalo/ENV/bin/activate
#python scripts/generate_particles_no_base.py
#python scripts/merge_geo.py data/block_equilibrium_0.geo walls_no_base.geo data/block_con_caras.geo
#python scripts/cpt2geo.py
#python scripts/geo2vtu.py
python scripts/PDFforce_components.py
#python scripts/FULLanalysis.py
