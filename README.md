This repository contains simulation data and post-processing scripts to reproduce all figures from the paper studying stress transmission, force anisotropy, and hysteresis in quasi-statically loaded granular
   assemblies.


Overview
---------------
  The dataset comes from a 3D Discrete Element Method (DEM) simulation of a granular block undergoing cyclic compression-decompression. The simulation was run with ESyS-Particle on an MPI cluster.
  Post-processing scripts extract stress tensors, compute force statistics, and generate publication-ready figures.

  Physical setup:
  - Box: 0.4 × 0.397 × 0.4 m filled with polydisperse spheres (radius 3–7 mm, density 3000 kg/m³)
  - Loading: 4 cyclic compression-decompression cycles, 12 mm displacement amplitude
  - Interaction: Hertzian viscoelastic contact law (E = 46 GPa, ν = 0.245, μ = 0.16)
  - Three spatial sampling regions: Center, Medium, Corner — each split into top, middle, and bottom sub-zones

  ---
  Repository Structure

  Quasi-static-Data-main/
  
  └── Simulation/
  
      ├── Data/
      
      │   ├── Forces/
      
      │   │   ├── ForcesFilt/          # filtered_partForce.*.dat  (961 files)
      
      │   │   ├── xyz_Forces_Center.csv
      
      │   │   ├── xyz_Forces_Medium.csv
      
      │   │   └── xyz_Forces_Corner.csv
      
      │   ├── Stress/
      
      │   │   └── data.*.txt           # stress tensor fields (955 files)
      
      │   └── Walls/
      
      │       ├── floorPosition.dat
      
      │       ├── floorForce.dat
      
      │       ├── roofForce.dat
      
      │       ├── x+WallForce.dat / x-WallForce.dat
      
      │       └── z+WallForce.dat / z-WallForce.dat
      
      ├── Data-Processing/             # MATLAB analysis scripts
      
      │   ├── ForceComponentsAnalysis.m
      
      │   ├── StressAnalysis.m
      
      │   ├── WallForcesAnalysis.m
      
      │   └── errormc.m               # helper: polynomial fit with errors
      
      ├── Execution/                   # Simulation and data-extraction scripts
      │   ├── blockCompression.py      # main DEM simulation
      │   ├── extract_data.py          # VTK → plain-text converter
      │   ├── PDFforce_components.py   # spatial force statistics
      │   ├── stress3vti.py            # stress tensor on regular grid
      │   ├── Parameters_blockCompression.py
      │   └── Parameters_general.py
      └── Figures/                     # output figures (PNG, 300 DPI)
          ├── FloorForcevsStrain.png
          ├── Fx/Fy/Fzmax.png, Fx/Fy/Fzmean.png
          ├── InternalSressvsConfPressure_lateral.png
          ├── Center/
          ├── Medium/
          └── Corner/

Requirements

  MATLAB (figure generation)

  - MATLAB R2019b or later
  - No additional toolboxes required (uses built-in polyfit, movmean, smooth)

  Python (simulation and data extraction)

  - Python 3.7+
  - https://launchpad.net/esys-particle (for running blockCompression.py)
  - numpy, vtk (for extract_data.py and stress3vti.py)

  pip install numpy vtk

  ---
  Reproducing the Figures

  All figure-generating scripts are in Simulation/Data-Processing/. Run them from MATLAB with the working directory set to Simulation/. The scripts read data from Data/ and write figures to Figures/.

  Figure: Macroscopic stress–strain curve

  Script: WallForcesAnalysis.m
  Reads: Data/Walls/floorPosition.dat, floorForce.dat, roofForce.dat, x±WallForce.dat, z±WallForce.dat
  Output: Figures/FloorForcevsStrain.png

  Computes macroscopic strain ε = (L₀ − floor displacement)/L₀ and wall stress σ = Force / 0.16 m². Plots compression and decompression branches separately with cycle-averaged curves.

Data File Formats

  filtered_partForce.*.dat

  ESyS-Particle RAW_WITH_POS_ID format. Each row is one contact:
  particle_id_1  particle_id_2  x  y  z  Fx  Fy  Fz

  data.*.txt

  Plain ASCII, one grid cell per row:
  x  y  z  sigma_xx  sigma_xy  sigma_xz  sigma_yy  sigma_yz  sigma_zz

  xyz_Forces_*.csv

  timestep, zone, N, Fy_mean, Fy_std, Fy_max, frac_strong_y,
  Fx_mean, Fx_std, Fx_max, frac_strong_x,
  Fz_mean, Fz_std, Fz_max, frac_strong_z, anis_xz
  - frac_strong: fraction of contacts with F > ⟨F⟩ (strong network)
  - anis_xz: horizontal anisotropy index (⟨Fx²⟩ − ⟨Fz²⟩)/(⟨Fx²⟩ + ⟨Fz²⟩)

  Wall files

  Two-column ASCII: simulation time and scalar value.


Re-running the Simulation

 Note: Requires a cluster with ESyS-Particle. Not needed to reproduce figures from pre-computed data in Data/.

  cd Simulation/Execution
  mpirun -np 18 python blockCompression.py <geometry_file>

  After simulation, convert to analysis-ready formats:
  python stress3vti.py partForce.N.dat output.vti Xmin Xmax Ymin Ymax Zmin Zmax Nx Ny Nz
  python extract_data.py       # converts .vti → data.*.txt
  python PDFforce_components.py  # computes xyz_Forces_*.csv


Contact

  For questions about the data or scripts, please open an issue in this repository.
