import numpy as np
from math import inf, pi
import meshio
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
import copy
import gmsh
import math
import sys

def write_to_vtk_with_faces(mshfile):
    # Temporarily store mesh_data in copy:
    Mesh = meshio.read(mshfile)
    mesh = copy.copy(Mesh)

    available_geometries = ['hexahedron', 'wedge', 'tetra', 'quad', 'triangle']

    cell_property = ['CellEntityIds']
    props_num = len(cell_property)

    # Matrix
    geom_id = 0
    Mesh.cells = {}
    cell_data = {}
    for ith_geometry in mesh.cells_dict.keys():
        if ith_geometry in available_geometries:
            Mesh.cells[ith_geometry] = mesh.cells_dict[ith_geometry]
            # Add matrix data to dictionary:
            for i in range(props_num):
                if cell_property[i] not in cell_data: cell_data[cell_property[i]] = []
                cell_data[cell_property[i]].append(np.abs(np.array(mesh.cell_data_dict['gmsh:physical'][ith_geometry], dtype=np.int64), dtype=np.int64))
        geom_id += 1

    #calc_tetra_volumes(mesh.cells[1].data[mesh.cell_data_dict['gmsh:physical']['tetra'] == 97], mesh.points)

    # Store solution for each time-step:
    mesh = meshio.Mesh(
        Mesh.points,
        Mesh.cells,
        cell_data=cell_data)
    meshio.write(mshfile.split('.')[0] + '.vtu', mesh)

    return 0

write_to_vtk_with_faces(mshfile='new_setup.msh')
