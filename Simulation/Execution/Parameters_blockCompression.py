# Original script developed by Gonzalo Tancredi
# Modified and adapted by Noelia Olivera Rodríguez (2026)
#
# This file defines the temporal parameters and loading protocol
# for cyclic compression–decompression simulations of bonded
# granular blocks using ESyS-Particle.

from Parameters_general import *
# Imports global geometric and material parameters


# Cyclic loading parameters
floordisp = 0.012        # Vertical displacement imposed per cycle (m)
num_cycles = 4          # Number of compression–decompression cycles
num_phases = num_cycles * 2
# Each cycle consists of:
#   1 compression phase + 1 decompression phase

# Temporal control

te = 1.8                # Total simulation time (s)
trelax = 0              # Initial relaxation time (s)

# Duration of each loading phase
time_per_phase = (te - trelax) / num_phases


# Required wall velocity to achieve floordisp within one phase
velCompWall = floordisp / time_per_phase


# Time discretization

nsteps = 180                        # Number of increments per phase
tup = time_per_phase / nsteps       # Duration of each increment

startinc = tup * 2.0                # Ramp-up time increment
stopinc = tup / 2.0                 # Ramp-down time increment

# Snapshot interval (output frequency)
t_snap = tup


# Explicit time step (DEM stability condition scaling)
dt = 1.e-1 * 2. * 5.e-5 * minParticleRadius
# Time step proportional to particle size

# Domain decompositio
# Number of subdomains in each spatial direction
# (must match SLURM/MPI configuration)

numSubdomainsX = 3
numSubdomainsY = 2
numSubdomainsZ = 3


print floordisp, velCompWall, nsteps, startinc, stopinc, tup, dt, te, t_snap
# Prints key parameters for runtime verification
