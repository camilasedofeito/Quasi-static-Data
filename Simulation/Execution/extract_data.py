# Original workflow concept, conversion and adaptation to VTK-based stress extraction:
# Noelia Olivera Rodríguez (2026)
# This script converts ESyS-Particle VTK image data files (.vti)
# containing cell-based stress tensor components into plain-text (.txt)
# files. Each output file contains:
#   x  y  z  sigma_xx  sigma_xy  sigma_xz  sigma_yy  sigma_yz  sigma_zz
# where (x, y, z) are cell-center coordinates and the remaining columns
# correspond to the six independent components of the Cauchy stress tensor.
# The script processes all valid stress*.vti files in a given directory.


import vtk
import numpy as np
import os

# FUNCTION: Read VTI file and export stress tensor to TXT
def read_vti_and_save(filename, output_file):
    """
    Reads a .vti file containing cell-based stress components and
    exports a text file with cell-center coordinates and stress tensor values.
    """

    # Load VTK image data
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()

    # Extract stress tensor components
    cell_data = data.GetCellData()

    components = [
        "sigma_xx", "sigma_xy", "sigma_xz",
        "sigma_yy", "sigma_yz", "sigma_zz"
    ]

    stress_data = []

    for component in components:
        array = cell_data.GetArray(component)

        # Extract scalar value for each cell
        stress_data.append(
            np.array([array.GetTuple1(i)
                      for i in range(array.GetNumberOfTuples())])
        )

    # Shape: (n_cells, 6)
    stress_data = np.array(stress_data).T


    # Compute cell-center coordinates
    coords = np.array([
        data.GetCell(i).GetBounds()
        for i in range(data.GetNumberOfCells())
    ])

    # Bounds format: [x0, x1, y0, y1, z0, z1]
    cell_centers = np.array([
        [(b[0] + b[1]) / 2,
         (b[2] + b[3]) / 2,
         (b[4] + b[5]) / 2]
        for b in coords
    ])


    # Combine spatial and stress data
    combined_data = np.hstack((cell_centers, stress_data))


    # Save to text file
    header = "x y z sigma_xx sigma_xy sigma_xz sigma_yy sigma_yz sigma_zz"

    np.savetxt(output_file, combined_data,
               header=header, fmt="%.6f")

    print(f"Data saved in {output_file}")

# Batch processing of stress files
# The script scans the specified directory and processes all files:
#   - Starting with 'stress'
#   - Ending with '.vti'
#   - Excluding transient files ending in 't.vti'

dir = 'C:\\Users\\Noeli\\Desktop\\stress\\stress3CiclosBajo\\'

contents = os.listdir(dir)

for file in contents:

    if (os.path.isfile(os.path.join(dir, file)) and
        file.endswith('.vti') and
        file.startswith('stress') and
        not file.endswith('t.vti')):

        print(f"Processing: {file}")

        filename = os.path.join(dir, file)

        outname = file.split('stress.')[1].split('.vti')[0]

        output_file = (
            f"C:\\Users\\Noeli\\Desktop\\stress\\data3CiclosBajo\\"
            f"data.{outname}.txt"
        )

        read_vti_and_save(filename, output_file)
