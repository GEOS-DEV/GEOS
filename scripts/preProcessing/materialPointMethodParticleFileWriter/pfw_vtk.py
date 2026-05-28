# -*- coding: utf-8 -*-
import numpy as np                   
import os                                  
import importlib
import re
import math
import argparse
import xml.etree.ElementTree as ET
import sys

# Example usage:
# python3 pfw_vtk.py "/g/g19/homel1/geosxRuns/compactionAndTension_r10" -t 8.0 -r -f "particleCenter" "particleMass" "particleVolume" "particleStress" "particleDamage"

def read_vtk(jobDir, vtk_fields, timesteps=[], read_grid_metadata = True):
    # Read first time step values
    # vtkOutput.pvd
    # vktOutput/
    # |
    # --- 000000/
    # --- 000000.vtm
    # --- ...

    jobName = os.path.basename(os.path.normpath(jobDir))
    # print("Job Directory", jobDir)
    # print("Job Name", jobName)
    vtk_filepath = os.path.join(jobDir, "vtkOutput.pvd")
    vtk_dirpath = os.path.join(jobDir, "vtkOutput")

    data = {} # Store particle data as a dictionary, each field has key with name
    assert vtk_fields is not None, "Must specify which fields to read from vtk"

    # Split particle and grid fields, they have different vtk subdirectories
    particle_fields = []
    grid_fields = []
    for f in vtk_fields:
        if f[:min(len(f),8)] == "particle":
            particle_fields.append(f)
        if f[:min(len(f),4)] == "grid":
            grid_fields.append(f)

    readAll = False
    cycle_times_to_read = []
    print("Reading timesteps = {" + ", ".join(timesteps) + "}")
    if len(timesteps) == 0:
        readAll = True
    else:
        cycle_times_to_read = [ float(i) for i in timesteps]

    vtk_tree = ET.parse(vtk_filepath)
    vtk_root = vtk_tree.getroot()
    vtk_cycles = vtk_root.findall(".//Collection//")

    for vtk_cycle in vtk_cycles:

        cycle_file = vtk_cycle.attrib["file"] 
        cycle_time = float(vtk_cycle.attrib["timestep"])

        if np.all(np.absolute(np.array(cycle_times_to_read) - cycle_time) > 1e-12 ) and not readAll:
            continue

        timestep_data = {}
        for pf in particle_fields:
            timestep_data[pf] = []
        # vtm_string = "{:d}".format(cycle).zfill(6)

        print("Cycle time=", cycle_time)

        vtm_filepath = os.path.join(jobDir, cycle_file) 
        print(vtm_filepath)
        vtm_tree = ET.parse(vtm_filepath)
        vtm_root = vtm_tree.getroot()
        
        # numRanks = len(ranks)/numParticleRegions

        # Count number of particles first and resize 
        # Read particle fields
        if particle_fields:
            particle_regions = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='particles']//Block//Block//Block")
            numParticleRegions = len(particle_regions)
            particle_ranks = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='particles']//Block//Block//Block//")
            for rank in particle_ranks:
                vtu_filepath = os.path.join(jobDir, "vtkOutput", rank.attrib['file'])
                print("\t" + vtu_filepath)
                vtu_tree = ET.parse(vtu_filepath)
                vtu_root = vtu_tree.getroot()

                # Loop over particle fields and assign data
                for i, pf in enumerate(particle_fields):
                    particle_field_nodes = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='" + pf + "']")
                    print(particle_field_nodes)
                    if particle_field_nodes:
                        particle_field_node = particle_field_nodes[0]
                        fieldType = particle_field_node.attrib["type"]
                        if "NumberOfComponents" in particle_field_node.attrib:
                            numComponents = int(particle_field_node.attrib["NumberOfComponents"])
                        else:
                            numComponents = 1

                        if fieldType == "Float64":
                            elements = [float(d) for d in particle_field_node.text.split()]
                        elif fieldType == "Int64":
                            elements = [int(d) for d in particle_field_node.text.split()]
                        else:
                            raise Exception("Unknown vtk field type reading field:", pf)
                        
                        ## We use the first field of every rank to compute the number of particles per rank
                        ## We track this for every timestep for output and to check for particle deletion
                        # if i == 0:
                        #     timestep_data["numParticles"] += int(len(elements)/numComponents)

                        for i in range(0,len(elements),numComponents):
                            timestep_data[pf].append(elements[i:i+numComponents])

        # These are required to read meta data from grid (e.g. element centers, nodal positions, and connectivity)
        grid_regions = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='backgroundGrid']//Block//Block//Block")
        numGridRegions = len(grid_regions)
        grid_ranks = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='backgroundGrid']//Block//Block//Block//")
        if grid_fields:
            for gf in grid_fields:
                timestep_data[gf] = []

            for rank in grid_ranks:
                vtu_filepath = os.path.join(jobDir, "vtkOutput", rank.attrib['file'])
                print("\t" + vtu_filepath)
                vtu_tree = ET.parse(vtu_filepath)
                vtu_root = vtu_tree.getroot()

                # Loop over particle fields and assign data
                # Need to check if field is point or cell data
                for gf in grid_fields:
                    grid_field_node = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='" + gf + "']")
                    if len(grid_field_node) == 0:
                        grid_field_node = vtu_root.findall(".//UnstructuredGrid//Piece//PointData//DataArray[@Name='" + gf + "']")
                        assert len(grid_field_node) > 0, "Could not find backgroundGrid field " + gf + " in vtk output"

                    grid_field_node = grid_field_node[0]

                    fieldType = grid_field_node.attrib["type"]
                    if "NumberOfComponents" in grid_field_node.attrib:
                        numComponents = int(grid_field_node.attrib["NumberOfComponents"])
                    else:
                        numComponents = 1

                    if fieldType == "Float64":
                        elements = [float(d) for d in grid_field_node.text.split()]
                    elif fieldType == "Int64":
                        elements = [int(d) for d in grid_field_node.text.split()]
                    else:
                        raise Exception("Unknown vtk field type reading field:", gf)
                    
                    for i in range(0,len(elements),numComponents):
                        timestep_data[gf].append(elements[i:i+numComponents])

        # Import cell mesh separately from background grid
        # connectivity
        # element center
        # node positions
        if read_grid_metadata:
            timestep_data["cellCenters"] = []
            timestep_data["nodalPositions"] = []
            timestep_data["connectivity"] = []
            connectivityOffset = 0
            for rank in grid_ranks:
                vtu_filepath = os.path.join(jobDir, "vtkOutput", rank.attrib['file'])
                print("\t" + vtu_filepath)
                vtu_tree = ET.parse(vtu_filepath)
                vtu_root = vtu_tree.getroot()

                # Read elementCenters
                vtk_node = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='elementCenter']")
                numComponents = 3

                cell_centers = [float(d) for d in vtk_node[0].text.split()]
                for i in range(0, len(cell_centers), numComponents):
                    timestep_data["cellCenters"].append(cell_centers[i:i+numComponents])

                # Read nodal positions
                vtk_node = vtu_root.findall(".//UnstructuredGrid//Piece//Points//DataArray[@Name='Points']")
                numComponents = 3

                nodal_positions = [float(d) for d in vtk_node[0].text.split()]
                for i in range(0, len(nodal_positions), numComponents):
                    timestep_data["nodalPositions"].append(nodal_positions[i:i+numComponents])

                # Read connectivity of element nodal connectivity
                vtk_node = vtu_root.findall(".//UnstructuredGrid//Piece//Cells//DataArray[@Name='connectivity']")
                numComponents = 8  # Assumes Hexahedral element

                connectivity = [int(d) + connectivityOffset for d in vtk_node[0].text.split()]
                for i in range(0, len(connectivity), numComponents):
                    timestep_data["connectivity"].append(connectivity[i:i+numComponents])
                connectivityOffset += int(len(connectivity)/numComponents)

        data[cycle_time] = timestep_data

    return data


def output_data(data, output="output"):
        # OUTPUT data to files

        # Post process particle and grid fields
        # Output data for all timesteps read
        times = list(data.keys())
        times.sort()

        # Use fields from first timestep read
        fields = list(data[times[0]].keys())       

        p_fields = [x for x in fields if x in particle_fields]
        file_name = output + "_particles"
        num_fields = len(p_fields)
        if num_fields > 0:
            num_particles = len(data[times[0]][p_fields[0]])
            with open(file_name, 'w') as file:
                for t in times:
                    file.write("Timestep: " + str(t) + "\n")
                    for k, f in enumerate(p_fields):
                        file.write(f + "\n")
                        for p in range(num_particles):
                            pdat = data[t][f][p]
                            # if hasattr(p, '__iter__'):
                            file.write(", ".join([str(i) for i in pdat]) + "\n")
            print("Output particle data to", file_name)

        g_fields = [x for x in fields if x in grid_fields]
        file_name = output + "_grid"
        num_fields = len(g_fields)
        if num_fields > 0:
            num_nodes = len(data[times[0]][g_fields[0]])
            with open(file_name, 'w') as file:
                for t in times:
                    file.write("Timestep: " + str(t) + "\n")
                    for k, f in enumerate(g_fields):
                        file.write(f + "\n")
                        for p in range(num_nodes):
                            pdat = data[t][f][p]
                            file.write(", ".join([str(i) for i in pdat]) + "\n")
            print("Output grid data to", file_name)

        m_fields = ["cellCenters", "nodalPositions", "connectivity"]
        file_name = output + "_grid_metadata"
        num_cells = len(data[times[0]][m_fields[0]]) # number of elements
        num_nodes = len(data[times[0]][m_fields[1]]) # grid nodes
        with open(file_name, 'w') as file:
            for t in times:
                file.write("Timestep: " + str(t) + "\n")
                file.write("cellCenters\n")
                for p in range(num_cells):
                    file.write(", ".join([str(i) for i in data[t]["cellCenters"][p]]) + "\n")
                file.write("nodalPositions\n")
                for p in range(num_nodes):
                    file.write(", ".join([str(i) for i in data[t]["nodalPositions"][p]]) + "\n")
                file.write("connectivity\n")
                for p in range(num_cells):
                    file.write(", ".join([str(i) for i in data[t]["connectivity"][p]]) + "\n")

        print("Output grid metadata to", file_name)


def read_data(jobDir, particle_fields, grid_fields, timesteps):
    data = {} # Store particle data as a dictionary, each field has key with name

    readAll = False
    # cycle_times_to_read = []
    # print("Reading timesteps = {" + ", ".join(timesteps) + "}")
    # if len(timesteps) == 0:
    #     readAll = True
    # else:
    #     cycle_times_to_read = [ float(i) for i in timesteps]

    cycle_times_to_read = timesteps

    vtk_filepath = os.path.join(jobDir, "vtkOutput.pvd")
    vtk_tree = ET.parse(vtk_filepath)
    vtk_root = vtk_tree.getroot()
    vtk_cycles = vtk_root.findall(".//Collection//")

    for vtk_cycle in vtk_cycles:

        cycle_file = vtk_cycle.attrib["file"] 
        cycle_time = float(vtk_cycle.attrib["timestep"])

        if np.all(np.absolute(np.array(cycle_times_to_read) - cycle_time) > 1e-12 ) and not readAll:
            continue

        timestep_data = {}
        for pf in particle_fields:
            timestep_data[pf] = []
        # vtm_string = "{:d}".format(cycle).zfill(6)

        # print("Cycle time=", cycle_time)

        vtm_filepath = os.path.join(jobDir, cycle_file) 
        # print(vtm_filepath)
        vtm_tree = ET.parse(vtm_filepath)
        vtm_root = vtm_tree.getroot()
        
        # numRanks = len(ranks)/numParticleRegions

        # Count number of particles first and resize 
        # Read particle fields
        if particle_fields:
            particle_regions = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='particles']//Block//Block//Block")
            numParticleRegions = len(particle_regions)
            particle_ranks = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='particles']//Block//Block//Block//")
            for rank in particle_ranks:
                vtu_filepath = os.path.join(jobDir, "vtkOutput", rank.attrib['file'])
                # print("\t" + vtu_filepath)
                vtu_tree = ET.parse(vtu_filepath)
                vtu_root = vtu_tree.getroot()

                # Loop over particle fields and assign data
                for i, pf in enumerate(particle_fields):
                    particle_field_node = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='" + pf + "']")[0]
                    fieldType = particle_field_node.attrib["type"]
                    if "NumberOfComponents" in particle_field_node.attrib:
                        numComponents = int(particle_field_node.attrib["NumberOfComponents"])
                    else:
                        numComponents = 1

                    if fieldType == "Float64":
                        elements = [float(d) for d in particle_field_node.text.split()]
                    elif fieldType == "Int64" or fieldType == "Int32": 
                        elements = [int(d) for d in particle_field_node.text.split()]
                    else:
                        raise Exception("Unknown vtk field type reading field:", pf)
                    
                    ## We use the first field of every rank to compute the number of particles per rank
                    ## We track this for every timestep for output and to check for particle deletion
                    # if i == 0:
                    #     timestep_data["numParticles"] += int(len(elements)/numComponents)

                    for i in range(0,len(elements),numComponents):
                        timestep_data[pf].append(elements[i:i+numComponents])

        # These are required to read meta data from grid (e.g. element centers, nodal positions, and connectivity)
        grid_regions = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='backgroundGrid']//Block//Block//Block")
        numGridRegions = len(grid_regions)
        grid_ranks = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='backgroundGrid']//Block//Block//Block//")
        if grid_fields:
            for gf in grid_fields:
                timestep_data[gf] = []

            for rank in grid_ranks:
                vtu_filepath = os.path.join(jobDir, "vtkOutput", rank.attrib['file'])
                # print("\t" + vtu_filepath)
                vtu_tree = ET.parse(vtu_filepath)
                vtu_root = vtu_tree.getroot()

                # Loop over particle fields and assign data
                # Need to check if field is point or cell data
                for gf in grid_fields:
                    # grid_field_node = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='" + gf + "']")
                    grid_field_node = vtu_root.findall(".//UnstructuredGrid//Piece//PointData//DataArray[@Name='" + gf + "']")
                    if len(grid_field_node) == 0:
                        grid_field_node = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='" + gf + "']")
                        # grid_field_node = vtu_root.findall(".//UnstructuredGrid//Piece//PointData//DataArray[@Name='" + gf + "']")
                        assert len(grid_field_node) > 0, "Could not find backgroundGrid field " + gf + " in vtk output"

                    grid_field_node = grid_field_node[0]

                    fieldType = grid_field_node.attrib["type"]
                    if "NumberOfComponents" in grid_field_node.attrib:
                        numComponents = int(grid_field_node.attrib["NumberOfComponents"])
                    else:
                        numComponents = 1

                    if fieldType == "Float64":
                        elements = [float(d) for d in grid_field_node.text.split()]
                    elif fieldType == "Int64" or fieldType == "Int32":
                        elements = [int(d) for d in grid_field_node.text.split()]
                    else:
                        raise Exception("Unknown vtk field type reading field:", gf)
                    
                    for i in range(0,len(elements),numComponents):
                        timestep_data[gf].append(elements[i:i+numComponents])

        # Import cell mesh separately from background grid
        # connectivity
        # element center
        # node positions
        read_grid_metadata = True
        if read_grid_metadata:
            timestep_data["cellCenters"] = []
            timestep_data["nodalPositions"] = []
            timestep_data["connectivity"] = []
            connectivityOffset = 0
            for rank in grid_ranks:
                vtu_filepath = os.path.join(jobDir, "vtkOutput", rank.attrib['file'])
                # print("\t" + vtu_filepath)
                vtu_tree = ET.parse(vtu_filepath)
                vtu_root = vtu_tree.getroot()

                # Read elementCenters
                vtk_node = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='elementCenter']")
                numComponents = 3

                cell_centers = [float(d) for d in vtk_node[0].text.split()]
                for i in range(0, len(cell_centers), numComponents):
                    timestep_data["cellCenters"].append(cell_centers[i:i+numComponents])

                # Read nodal positions
                vtk_node = vtu_root.findall(".//UnstructuredGrid//Piece//Points//DataArray[@Name='Points']")
                numComponents = 3

                nodal_positions = [float(d) for d in vtk_node[0].text.split()]
                for i in range(0, len(nodal_positions), numComponents):
                    timestep_data["nodalPositions"].append(nodal_positions[i:i+numComponents])

                # Read connectivity of element nodal connectivity
                vtk_node = vtu_root.findall(".//UnstructuredGrid//Piece//Cells//DataArray[@Name='connectivity']")
                numComponents = 8  # Assumes Hexahedral element

                connectivity = [int(d) + connectivityOffset for d in vtk_node[0].text.split()]
                for i in range(0, len(connectivity), numComponents):
                    timestep_data["connectivity"].append(connectivity[i:i+numComponents])
                connectivityOffset += int(len(connectivity)/numComponents)

        data[cycle_time] = timestep_data
    return data


def read_fields(jobDir):
    particle_fields = []
    grid_fields = []

    cycle_times_to_read = [0.0]

    vtk_filepath = os.path.join(jobDir, "vtkOutput.pvd")
    vtk_tree = ET.parse(vtk_filepath)
    vtk_root = vtk_tree.getroot()
    vtk_cycles = vtk_root.findall(".//Collection//")

    for vtk_cycle in vtk_cycles:

        cycle_file = vtk_cycle.attrib["file"] 
        cycle_time = float(vtk_cycle.attrib["timestep"])

        # Only need to read first timestep and first rank (all ranks and timesteps should have the same fields)
        if np.all(np.absolute(np.array(cycle_times_to_read) - cycle_time) > 1e-12 ):
            continue
        
        print("Cycle", vtk_cycle)

        vtm_filepath = os.path.join(jobDir, cycle_file) 
        vtm_tree = ET.parse(vtm_filepath)
        vtm_root = vtm_tree.getroot()

        # Read particle fields
        particle_regions = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='particles']//Block//Block//Block")
        numParticleRegions = len(particle_regions)
        particle_ranks = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='particles']//Block//Block//Block//")
        
        print("VTK Particle Fields:")
        for i, rank in enumerate(particle_ranks):
            # Dirty method to exit loop after first iteration
            if i > 0:
                break

            vtu_filepath = os.path.join(jobDir, "vtkOutput", rank.attrib['file'])
            vtu_tree = ET.parse(vtu_filepath)
            vtu_root = vtu_tree.getroot()

            pnodes = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray")
            for f in pnodes:
                print("\t" + f.attrib["Name"])
                
        # These are required to read meta data from grid (e.g. element centers, nodal positions, and connectivity)
        grid_regions = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='backgroundGrid']//Block//Block//Block")
        numGridRegions = len(grid_regions)
        grid_ranks = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='backgroundGrid']//Block//Block//Block//")

        print("VTK Grid Fields:")
        for i, rank in enumerate(grid_ranks):
            # Dirty method to exit loop after first iteration
            if i > 0:
                break
                
            vtu_filepath = os.path.join(jobDir, "vtkOutput", rank.attrib['file'])
            vtu_tree = ET.parse(vtu_filepath)
            vtu_root = vtu_tree.getroot()

            gnodes = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray")
            for f in gnodes:
                print("\t" + f.attrib["Name"])

            gnodes = vtu_root.findall(".//UnstructuredGrid//Piece//PointData//DataArray")
            for f in gnodes:
                print("\t" + f.attrib["Name"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Manipulate VTK output')
    parser.add_argument('jobdir', help="directory with vtk output")
    parser.add_argument('-r', '--read', action='store_true', help="display plot interactively")
    parser.add_argument('-t', '--timesteps', default=[], nargs='+', help="specify which timesteps to read")
    parser.add_argument('-l', '--list', action='store_true', default=False, help="List all fields in vtkOutput")
    parser.add_argument('-f', '--fields', nargs='+', help="specify which fields to read")
    parser.add_argument('-c', '--clean', action='store_true', default=False, help="directory with vtk to clean")
    parser.add_argument('-b', '--rebuild', action='store_true', default=False, help="rebuild vtkOutput.pvd file from slurm output and vtkOutput subdirectories")
    parser.add_argument('-d', '--delete', action='store_true', default=False, help="delete fields from vtkOutput")
    parser.add_argument('-o', '--output', default="output", help="Name of output file to save particle data")
    args = parser.parse_args()
    print("Command line args",args)

    # Removes extra plots output before job times out and is restarted
    # Run from command line by
    # python3 pfw_clean_vtk.py <job directory>
    jobDir = args.jobdir

    jobName = os.path.basename(os.path.normpath(jobDir))
    print("Job Directory", jobDir)
    print("Job Name", jobName)
    vtk_filepath = os.path.join(jobDir, "vtkOutput.pvd")
    vtk_dirpath = os.path.join(jobDir, "vtkOutput")

    assert os.path.isdir(jobDir), "Input must be a job directory"
    assert os.path.isfile(vtk_filepath), "Job directory has no vtkOutput.pvd file!"
    assert os.path.isdir(vtk_dirpath), "Job directory has no vtkOutput directory!"

    work_dir = os.path.dirname(os.path.realpath(__file__))
    print("Working directory =", work_dir)

    # Rebuild vtkOutput.pvd
    if args.rebuild:
        # Read number of vtk cycle outputs from vtkOutput directory
        vtk_dir_items = sorted([ int(f) for f in os.listdir(vtk_dirpath) if ".vtm" not in f ])
        output_cycles = dict.fromkeys(vtk_dir_items, None)
        # print(vtk_dir_items, output_cycles)

        # Search slurm files for time code for cycle number
        slurm_files = [ f for f in os.listdir(jobDir) if "slurm-" in f ]
        # Sort slurm_files with smallest job id first, ensures that last cycle entry corresponds to latest run, incase of overwritting of timestep data
        slurm_job_ids = [int(s.replace("slurm-","").replace('.out',"")) for s in slurm_files ]
        # print(slurm_job_ids)
        slurm_job_id_sorted_indices = np.argsort(slurm_job_ids)
        # print(slurm_job_id_sorted_indices)
        slurm_files = [ slurm_files[s] for s in slurm_job_id_sorted_indices ]
        # print(slurm_files)

        # Read slurm_files and search for cycles numbers
        for sf in slurm_files:
            with open(os.path.join(jobDir, sf), 'r') as sfile:
                sfile.readline()
                for line in sfile:
                    if 'Time:' in line and 'Cycle:' in line:
                        regmatch = re.search('Time: (.*) s, dt: (.*) s, Cycle: (.*)', line)
                        cycle_num = int(regmatch.group(3))
                        if cycle_num in output_cycles:
                            cycle_time = float(regmatch.group(1))
                            output_cycles[cycle_num] = cycle_time
                            # print(cycle_num, '=', cycle_time)
        # print(output_cycles)

        # Regenerate vtk text
        vtk_text = """<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1">
\t<Collection>
"""

        for cycle_num, cycle_time in output_cycles.items():
            print(cycle_num, cycle_time)
            if cycle_time is not None:
                vtk_text += """\t\t\t<DataSet timestep=""" + '"' + "{:.17f}".format(cycle_time) + '"' + """ file="vtkOutput/""" + str(cycle_num).zfill(6) + """.vtm" />\n"""
            else:
                print("Warning: No cycle time not found for cycle number", cycle_num, "in slurm outputs")
		
        vtk_text += """\t</Collection>
</VTKFile>"""

        # print(vtk_text)
        with open(vtk_filepath, 'w') as f:
            f.write(vtk_text)


    # Clean vtk (remove extra writes from restarts from vtkOutput.pvd file)
    if args.clean:
        # Get plotting interval from input
        inputFilePath = os.path.join(jobDir, "pfw_input_" + jobName + ".py")
        print("Py input file", inputFilePath)
        sys.path.insert(0,jobDir)
        inputFile = __import__("pfw_input_" + jobName)
        plotInterval = inputFile.pfw["plotInterval"]

        # Read vtkOutput.pvd
        print("Reading vtk output file")
        new_text = ""
        unique_times = set()
        with open(vtk_filepath, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if "DataSet" in line:
                    time = float(re.search('<DataSet timestep="(.*)" file=".*', line).group(1))
                    if abs(round(time/plotInterval % 1, 0) - (time/plotInterval % 1)) > 1e-2:
                        print(time, time/plotInterval % 1)
                        continue

                    if time in unique_times:
                        print(time, time/plotInterval % 1)
                        continue
                    
                    unique_times.add(time)

                new_text += line

        with open(vtk_filepath, 'w') as f:
            f.write(new_text)

    if args.delete:
        vtk_tree = ET.parse(vtk_filepath)
        vtk_root = vtk_tree.getroot()
        vtk_cycles = vtk_root.findall(".//Collection//")

        num_fields = len(args.fields)
        assert num_fields > 0, "Must specify one or more fields to delete from vtk output"
        found_fields = [False for f in range(num_fields)]

        for vtk_cycle in vtk_cycles:
            cycle_file = vtk_cycle.attrib["file"] 

            vtm_filepath = os.path.join(jobDir, cycle_file) 
            print(vtm_filepath)
            vtm_tree = ET.parse(vtm_filepath)
            vtm_root = vtm_tree.getroot()

            # First search the particle fields
            particle_ranks = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='particles']//Block//Block//Block//")
            for rank in particle_ranks:
                vtu_filepath = os.path.join(jobDir, "vtkOutput", rank.attrib['file'])
                print("\t" + vtu_filepath)
                vtu_tree = ET.parse(vtu_filepath)
                vtu_root = vtu_tree.getroot()

                # Loop over particle fields and assign data
                for ff, field in enumerate(args.fields):
                    parent_node = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='" + field + "']...") # Returns list
                    particle_field_node = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='" + field + "']") # Returns list
                    if particle_field_node:
                        found_fields[ff] = True
                        parent_node[0].remove(particle_field_node[0])
                vtu_tree.write(vtu_filepath)
            
            # Then check the grid fields
            grid_ranks = vtm_root.findall(".//vtkMultiBlockDataSet//Block[@name='backgroundGrid']//Block//Block//Block//")
            for rank in grid_ranks:
                vtu_filepath = os.path.join(jobDir, "vtkOutput", rank.attrib['file'])
                print("\t" + vtu_filepath)
                vtu_tree = ET.parse(vtu_filepath)
                vtu_root = vtu_tree.getroot()

                for ff, field in enumerate(args.fields):
                    parent_node = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='" + field + "']...") # Returns list
                    grid_field_node = vtu_root.findall(".//UnstructuredGrid//Piece//CellData//DataArray[@Name='" + field + "']") # Returns list
                    if grid_field_node:
                        found_fields[ff] = True
                        parent_node[0].remove(grid_field_node[0])
                vtu_tree.write(vtu_filepath)

        for ff, field in enumerate(args.fields):
            if not found_fields[ff]:
                print("Could not find field", field, "to delete")

    if args.list:
        # Lists particle and grid fields available in the vtk file
        read_fields(jobDir)


    # Read particle data fields from file and generate csv
    if args.read:
        # Read first time step values
        # vtkOutput.pvd
        # vktOutput/
        # |
        # --- 000000/
        # --- 000000.vtm
        # --- ...

        data = {} # Store particle data as a dictionary, each field has key with name
        assert args.fields is not None, "Must specify which fields to read from vtk"

        # Split particle and grid fields, they have different vtk subdirectories
        particle_fields = []
        grid_fields = []
        for f in args.fields:
            if f[:min(len(f),8)] == "particle":
                particle_fields.append(f)
            if f[:min(len(f),4)] == "grid":
                grid_fields.append(f)

        data = read_data(jobDir,  particle_fields, grid_fields, args.timesteps)

        # OUTPUT data to files
        # Post process particle and grid fields
        # Output data for all timesteps read
        times = list(data.keys())
        times.sort()

        # Use fields from first timestep read
        fields = list(data[times[0]].keys())       

        p_fields = [x for x in fields if x in particle_fields]
        file_name = args.output + "_particles"
        num_fields = len(p_fields)
        if num_fields > 0:
            num_particles = len(data[times[0]][p_fields[0]])
            with open(file_name, 'w') as file:
                for t in times:
                    file.write("Timestep: " + str(t) + "\n")
                    for k, f in enumerate(p_fields):
                        file.write(f + "\n")
                        for p in range(num_particles):
                            pdat = data[t][f][p]
                            # if hasattr(p, '__iter__'):
                            file.write(", ".join([str(i) for i in pdat]) + "\n")
            print("Output particle data to", file_name)

        g_fields = [x for x in fields if x in grid_fields]
        file_name = args.output + "_grid"
        num_fields = len(g_fields)
        if num_fields > 0:
            num_nodes = len(data[times[0]][g_fields[0]])
            with open(file_name, 'w') as file:
                for t in times:
                    file.write("Timestep: " + str(t) + "\n")
                    for k, f in enumerate(g_fields):
                        file.write(f + "\n")
                        for p in range(num_nodes):
                            pdat = data[t][f][p]
                            file.write(", ".join([str(i) for i in pdat]) + "\n")
            print("Output grid data to", file_name)

        m_fields = ["cellCenters", "nodalPositions", "connectivity"]
        file_name = args.output + "_grid_metadata"
        num_cells = len(data[times[0]][m_fields[0]]) # number of elements
        num_nodes = len(data[times[0]][m_fields[1]]) # grid nodes
        with open(file_name, 'w') as file:
            for t in times:
                file.write("Timestep: " + str(t) + "\n")
                file.write("cellCenters\n")
                for p in range(num_cells):
                    file.write(", ".join([str(i) for i in data[t]["cellCenters"][p]]) + "\n")
                file.write("nodalPositions\n")
                for p in range(num_nodes):
                    file.write(", ".join([str(i) for i in data[t]["nodalPositions"][p]]) + "\n")
                file.write("connectivity\n")
                for p in range(num_cells):
                    file.write(", ".join([str(i) for i in data[t]["connectivity"][p]]) + "\n")

        print("Output grid metadata to", file_name)
