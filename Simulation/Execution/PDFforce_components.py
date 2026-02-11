#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Force Component Analysis inside a Cylindrical Subdomain
# Author: Noelia Olivera Rodríguez (2026)
#
# This script performs a spatially resolved analysis of interparticle forces
# obtained from ESyS-Particle simulations (partForce.*.dat files).
# For each timestep:
#   • A cylindrical subdomain (defined in the x–z plane) is selected.
#   • Forces are decomposed into:
#         - Vertical component (Fy)
#         - Radial component (projected in the x–z plane)
#         - Tangential component
#   • The system is divided into three vertical regions:
#         bottom, middle, and top.
#   • For each region, the script computes:
#         - Mean force
#         - Standard deviation
#         - Maximum force
#         - Fraction of strong forces (F > ⟨F⟩)
#         - In-plane anisotropy index
#
# Results are exported to a CSV file for further statistical analysis.
#
# This post-processing step is used to characterize force transmission
# heterogeneity and anisotropy during cyclic compression of a bonded
# granular specimen.

import numpy as np
import matplotlib
matplotlib.use('Agg') # Non-interactive backend for HPC environments
import matplotlib.pyplot as plt
import os
import re
import csv
import sys

# Configuration parameters

# Directory containing partForce.*.dat files
folder = "/clusteruy/home/nolivera/ESyS/CompAndDescCyc/3CiclosBajo"

# Radial selection (cylindrical probe in x–z plane)
R_cyl = 0.04       # Cylinder radius [m]

# Cylinder center (user-defined)
x0 = 0.01           # x-coordinate of cylinder center [m]
z0 = 0.01           # z-coordinate of cylinder center [m]

# Vertical limits of the domain
ymin, ymax = 0.0, 0.3970

# Auxiliary functions
def timestep_number(filename):
     """
    Extract timestep number from filenames of the form:
    partForce.<STEP>.dat
    """
    m = re.search(r'partForce\.(\d+)\.dat', filename)
    return int(m.group(1)) if m else -1

#File collection
files = [f for f in os.listdir(folder)
         if f.startswith("partForce") and f.endswith(".dat")]
files = sorted(files, key=timestep_number)
print("Files found:", len(files))
sys.stdout.flush()

rows = []  # Storage for statistical results

# Main processing loop
for file in files:

    print("\nProcessing:", file)
    sys.stdout.flush()

    Fx_all, Fy_all, Fz_all = [], [], []
    x_all, z_all, y_all = [], [], []

    path = os.path.join(folder, file)
    
    # Read contact forces
    with open(path, "r") as f:
        for line in f:
            if len(line.split()) < 14:
                continue

            v = list(map(float, line.split()))
            # Contact position
            x, y, z = v[8], v[9], v[10]
            # Contact force components
            Fx, Fy, Fz = v[11], v[12], v[13]

           
            # Radial filtering in x–z plane
            dx = x - x0
            dz = z - z0
            r = np.sqrt(dx*dx + dz*dz)

            if r > R_cyl:
                continue

            Fx_all.append(Fx)
            Fy_all.append(Fy)
            Fz_all.append(Fz)
            x_all.append(x)
            z_all.append(z)
            y_all.append(y)

    if len(Fy_all) == 0:
        print("No valid data")
        continue

    # Convert to numpy arrays
    Fx_all = np.array(Fx_all)
    Fy_all = np.array(Fy_all)
    Fz_all = np.array(Fz_all)
    x_all  = np.array(x_all)
    z_all  = np.array(z_all)
    y_all  = np.array(y_all)

    # Radial projection
    # Decomposition of horizontal forces into radial direction
    dx_all = x_all - x0
    dz_all = z_all - z0

    r_mag = np.sqrt(dx_all**2 + dz_all**2)
    r_mag[r_mag == 0] = 1e-12 

    Fx_radial = (dx_all / r_mag) * Fx_all + (dz_all / r_mag) * Fz_all
    Fz_radial = Fy_all

    # Vertical segmentation
    # Domain divided into three equal vertical regions
    H = ymax - ymin
    zones = {
        "bottom": y_all < ymin + H/3,
        "middle": (y_all >= ymin + H/3) & (y_all < ymin + 2*H/3),
        "top": y_all >= ymin + 2*H/3
    }

    print("   Centro cilindro: x0 =", x0, " z0 =", z0)
    print("   y min =", np.min(y_all), " y max =", np.max(y_all))
    print("   bottom:", np.sum(zones["bottom"]),
          "middle:", np.sum(zones["middle"]),
          "top:", np.sum(zones["top"]))

    # Statistical analysis
    for zone, mask in zones.items():
        if np.sum(mask) == 0:
            continue

        Fy = Fz_radial[mask]   # vertical
        Fx = Fx_radial[mask]   # radial
        Fz = Fz_all[mask]      # tangencial
        # First-order statistics
        Fy_mean = np.mean(Fy)
        Fy_std  = np.std(Fy)
        Fy_max  = np.max(Fy)
        frac_strong = np.sum(Fy > Fy_mean) / float(len(Fy))

        Fx_mean = np.mean(Fx)
        Fx_std  = np.std(Fx)
        Fx_max  = np.max(Fx)
        frac_strong_x = np.sum(Fx > Fx_mean) / float(len(Fx))

        Fz_mean = np.mean(Fz)
        Fz_std  = np.std(Fz)
        Fz_max  = np.max(Fz)
        frac_strong_z = np.sum(Fz > Fz_mean) / float(len(Fz))

        N = len(Fy)
        
        # Anisotropy indicator in horizontal plane
        anis_xz = (np.mean(Fx**2) - np.mean(Fz**2)) / \
                  (np.mean(Fx**2) + np.mean(Fz**2))

        rows.append([
            timestep_number(file),
            zone,
            N,
            Fy_mean,
            Fy_std,
            Fy_max,
            frac_strong,
            Fx_mean,
            Fx_std,
            Fx_max,
            frac_strong_x,
            Fz_mean,
            Fz_std,
            Fz_max,
            frac_strong_z,
            anis_xz
        ])

        print(" ", zone, "<Fy> =", round(Fy_mean, 4))

# Output CSV
csv_path = os.path.join(folder, "xyz_Forces.csv")

with open(csv_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerow([
        "timestep",
        "zona",
        "N",
        "Fy_mean",
        "Fy_std",
        "Fy_max",
        "frac_strong",
        "Fx_mean",
        "Fx_std",
        "Fx_max",
        "frac_strong_x",
        "Fz_mean",
        "Fz_std",
        "Fz_max",
        "frac_strong_z",
        "anis_xz"
    ])
    for row in rows:
        writer.writerow(row)

print("\nANALYSIS COMPLETED")
print("CSV saved in:", csv_path)

