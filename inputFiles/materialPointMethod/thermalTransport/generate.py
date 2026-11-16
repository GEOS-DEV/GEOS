#!/usr/bin/env python3
"""Generate small, stationary GEOS MPM thermal verification jobs (standard library only)."""
import argparse
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def generate(destination, case="groups", cells=16, resistance=0.2, end_time=0.2):
    if cells < 4 or cells % 2 or resistance < 0 or end_time <= 0:
        raise ValueError("Use an even cells >= 4, nonnegative resistance, and positive end time")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    h = 1.0 / cells
    rows = max(2, cells // 4)
    height = rows * h
    contact = case in ("groups", "dfg", "isothermal", "insulated", "materials")
    two_materials = case == "materials"
    conductivity = 0.0 if case == "source" else 0.05
    columns = ["ID", "PositionX", "PositionY", "PositionZ", "ParticleType", "MaterialType",
               "ContactGroup", "SurfaceFlag", "Damage", "Temperature", "RVectorXX",
               "SurfaceNormalX", "SurfaceNormalY", "SurfaceNormalZ",
               "SurfacePositionX", "SurfacePositionY", "SurfacePositionZ"]
    with (destination / "particles.tsv").open("w") as stream:
        stream.write("\t".join(columns) + "\n")
        pid = 0
        for j in range(2 * rows):
            for i in range(2 * cells):
                pid += 1
                x, y, z = (i + 0.5) * h / 2, (j + 0.5) * h / 2, h / 2
                left = x < 0.5
                temperature = (400.0 if left else 200.0) if contact else 350.0
                if case == "bulk":
                    temperature = 350.0 + 50.0 * math.cos(math.pi * x)
                group = int(not left) if contact and case != "dfg" else 0
                row = [pid, x, y, z, 0, int(not left) if two_materials else 0,
                       group, 2 if contact else 0, 1.0 if case == "dfg" else 0.0,
                       temperature, (h / 2) ** 2 * h, (1.0 if left else -1.0) if contact else 0.0,
                       0, 0, (0.5 - x) if contact else 0.0, 0, 0]
                stream.write("\t".join(str(v) for v in row) + "\n")
    root = ET.Element("Problem")
    mesh = ET.SubElement(root, "Mesh")
    ET.SubElement(mesh, "InternalMesh", name="backgroundGrid", elementTypes="{ C3D8 }",
                  xCoords=f"{{ {-h}, {1+h} }}", yCoords=f"{{ {-h}, {height+h} }}",
                  zCoords=f"{{ {-h}, {2*h} }}", nx=f"{{ {cells+2} }}", ny=f"{{ {rows+2} }}",
                  nz="{ 3 }", periodic="{ 0, 0, 0 }", cellBlockNames="{ cb1 }")
    names = ["solidA", "solidB"] if two_materials else ["solidA"]
    blocks = ["pb0", "pb1"] if two_materials else ["pb0"]
    ET.SubElement(mesh, "ParticleMesh", name="particles", particleFile="particles.tsv",
                  particleBlockNames="{ " + ", ".join(blocks) + " }",
                  particleTypes="{ " + ", ".join("SinglePoint" for _ in blocks) + " }")
    regions = ET.SubElement(root, "ElementRegions")
    ET.SubElement(regions, "CellElementRegion", name="gridRegion", meshBody="backgroundGrid",
                  cellBlocks="{ cb1 }", materialList="{ null }")
    particle_regions = ET.SubElement(root, "ParticleRegions")
    for index, material in enumerate(names):
        ET.SubElement(particle_regions, "ParticleRegion", name=f"body{index}", meshBody="particles",
                      particleBlocks=f"{{ {blocks[index]} }}", materialList=f"{{ {material} }}")
    solvers = ET.SubElement(root, "Solvers", gravityVector="{ 0, 0, 0 }")
    target = "{ backgroundGrid/gridRegion, " + ", ".join(f"particles/body{i}" for i in range(len(names))) + " }"
    tracer = "{ " + ", ".join(f"{{ {x}, {height/2}, {h/2} }}" for x in (h/4, 0.25+h/4, 0.5-h/4, 0.5+h/4, 0.75-h/4, 1-h/4)) + " }"
    solver = ET.SubElement(solvers, "SolidMechanics_MPM", name="mpm", discretization="FE1",
                           targetRegions=target, timeIntegrationOption="ExplicitDynamic", updateMethod="FLIP",
                           planeStrain="1", cflFactor="0.25", initialDt="0.001", enableContact="0",
                           damageFieldPartitioning="1" if case == "dfg" else "0",
                           needsNeighborList="1" if case == "dfg" else "0", neighborRadius=str(2*h),
                           maxSingleFieldStateFractionForSeparability="-1", thinFeatureDFGThreshold="10",
                           boundaryConditionTypes="{ 1, 1, 1, 1, 1, 1 }", subdivideParticles="0",
                           tracerHistory="1", tracerCycleInterval="1", tracerCoordinates=tracer,
                           tracerVariables="{ particleID, mass, temperature }", tracerOutputPrefix="thermalTracer")
    ET.SubElement(solver, "MPMEvents")
    ET.SubElement(solver, "CohesiveZoneRegions")
    thermal = ET.SubElement(solver, "Thermal", enabled="1", defaultContactMode="Resistance" if contact else "Insulated",
                           defaultContactResistance=str(resistance), areaQuadrature="40", diagnostics="1")
    if case == "materials":
        ET.SubElement(thermal, "ContactResistance", name="materialPair", materials="{ solidA, solidB }", resistance=str(resistance))
    if case == "insulated":
        ET.SubElement(thermal, "ContactResistance", name="excluded", contactGroups="{ 0, 1 }", mode="Insulated")
    if case == "source":
        ET.SubElement(thermal, "HeatSource", name="uniformSource", specificPower="4")
    if case == "isothermal":
        ET.SubElement(thermal, "Isothermal", name="hotBath", contactGroups="{ 0 }", temperature="400")
    constitutive = ET.SubElement(root, "Constitutive")
    ET.SubElement(constitutive, "ElasticIsotropic", name="null", defaultDensity="1", defaultBulkModulus="1", defaultShearModulus="1")
    for name in names:
        ET.SubElement(constitutive, "ThermalElasticIsotropic", name=name, defaultDensity="1",
                      defaultBulkModulus="1", defaultShearModulus="1", thermalSpecificHeat="1",
                      thermalConductivity=str(conductivity), thermalReferenceTemperature="300")
    events = ET.SubElement(root, "Events", maxTime=str(end_time))
    ET.SubElement(events, "PeriodicEvent", name="advance", target="/Solvers/mpm")
    numerical = ET.SubElement(root, "NumericalMethods")
    finite = ET.SubElement(numerical, "FiniteElements")
    ET.SubElement(finite, "FiniteElementSpace", name="FE1", order="1")
    ET.indent(root)
    ET.ElementTree(root).write(destination / "thermal.xml", encoding="utf-8", xml_declaration=True)
    (destination / "case.json").write_text(json.dumps({"case": case, "cells": cells, "end_time": end_time,
                                                     "diffusivity": conductivity, "resistance": resistance}, indent=2) + "\n")
    return destination / "thermal.xml"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--case", choices=("source", "bulk", "groups", "dfg", "materials", "isothermal", "insulated"), default="groups")
    parser.add_argument("--cells", type=int, default=16)
    parser.add_argument("--resistance", type=float, default=0.2)
    parser.add_argument("--end-time", type=float, default=0.2)
    args = parser.parse_args()
    print(generate(args.output, args.case, args.cells, args.resistance, args.end_time))
