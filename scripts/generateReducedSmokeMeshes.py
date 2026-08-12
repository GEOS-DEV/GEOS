#!/usr/bin/env python3
"""Generate the reduced external meshes used by the long-running smoke tests."""

from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk


ROOT = Path(__file__).resolve().parents[1]


def read_grid(path: Path) -> vtk.vtkUnstructuredGrid:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = vtk.vtkUnstructuredGrid()
    grid.DeepCopy(reader.GetOutput())
    return grid


def write_grid(grid: vtk.vtkUnstructuredGrid, path: Path) -> None:
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(grid)
    writer.SetCompressorTypeToZLib()
    writer.SetDataModeToBinary()
    if not writer.Write():
        raise RuntimeError(f"Failed to write {path}")


def add_array(attributes, name: str, values, vtk_type=None) -> None:
    values = np.asarray(values)
    array = numpy_to_vtk(values, deep=True, array_type=vtk_type)
    array.SetName(name)
    attributes.AddArray(array)


def make_hybrid_property_mesh() -> None:
    source = ROOT / "inputFiles/solidMechanics/hybridHexPrismMesh.vtk"
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(str(source))
    reader.Update()
    grid = vtk.vtkUnstructuredGrid()
    grid.DeepCopy(reader.GetOutput())
    attribute = grid.GetCellData().GetArray("attribute")
    for cell_id in range(grid.GetNumberOfCells()):
        attribute.SetTuple1(cell_id, 0)

    centers_filter = vtk.vtkCellCenters()
    centers_filter.SetInputData(grid)
    centers_filter.Update()
    centers = np.array(
        [centers_filter.GetOutput().GetPoint(i) for i in range(grid.GetNumberOfCells())]
    )
    x, y, z = centers.T
    add_array(grid.GetCellData(), "PORO", 0.075 + 0.025 * (x + y + z) / 3.0)
    add_array(grid.GetCellData(), "NTG", np.full(len(centers), 0.9))
    add_array(grid.GetCellData(), "DENSITY", 2600.0 + 100.0 * z)
    add_array(grid.GetCellData(), "BULKMOD", 4.5e9 + 1.0e9 * x)

    output = ROOT / "inputFiles/singlePhaseFlow/pebi3d_with_properties_smoke.vtu"
    write_grid(grid, output)


def make_fault_mesh() -> None:
    mesh_dir = ROOT / "inputFiles/poromechanicsFractures/MESH"
    x_values = [-6000.0, -3600.0, -2400.0, -1200.0, 0.0, 1200.0, 2400.0, 3600.0, 6000.0]
    y_values = [-6000.0, -2400.0, -1200.0, 0.0, 1200.0, 6000.0]
    z_values = [
        -4000.0, -2256.0, -2038.0, -1820.0, -1810.0, -1805.0,
        -1800.0, -1780.0, -1650.0, -1575.0, -1350.0, -900.0, -450.0, 0.0,
    ]
    fault_x = x_values.index(1200.0)
    fault_y = range(y_values.index(-2400.0), y_values.index(1200.0) + 1)
    fault_z = range(z_values.index(-1810.0), z_values.index(-1575.0) + 1)

    points = vtk.vtkPoints()
    base_ids = {}
    for iz, z in enumerate(z_values):
        for iy, y in enumerate(y_values):
            for ix, x in enumerate(x_values):
                base_ids[ix, iy, iz] = points.InsertNextPoint(x, y, z)

    duplicate_ids = {}
    for iz in list(fault_z)[1:-1]:
        for iy in list(fault_y)[1:-1]:
            duplicate_ids[iy, iz] = points.InsertNextPoint(1200.0, y_values[iy], z_values[iz])

    domain = vtk.vtkUnstructuredGrid()
    domain.SetPoints(points)
    attributes = []
    young_modulus = []
    poisson_ratio = []

    def domain_point(ix, iy, iz, cell_ix):
        if ix == fault_x and cell_ix == fault_x and (iy, iz) in duplicate_ids:
            return duplicate_ids[iy, iz]
        return base_ids[ix, iy, iz]

    for iz in range(len(z_values) - 1):
        for iy in range(len(y_values) - 1):
            for ix in range(len(x_values) - 1):
                ids = [
                    domain_point(ix, iy, iz, ix),
                    domain_point(ix + 1, iy, iz, ix),
                    domain_point(ix + 1, iy + 1, iz, ix),
                    domain_point(ix, iy + 1, iz, ix),
                    domain_point(ix, iy, iz + 1, ix),
                    domain_point(ix + 1, iy, iz + 1, ix),
                    domain_point(ix + 1, iy + 1, iz + 1, ix),
                    domain_point(ix, iy + 1, iz + 1, ix),
                ]
                domain.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, ids)
                center_x = 0.5 * (x_values[ix] + x_values[ix + 1])
                center_y = 0.5 * (y_values[iy] + y_values[iy + 1])
                center_z = 0.5 * (z_values[iz] + z_values[iz + 1])
                in_fault_zone = (
                    -2400.0 <= center_y <= 1200.0
                    and -1810.0 <= center_z <= -1575.0
                )
                if in_fault_zone and -2400.0 <= center_x < 1200.0:
                    attributes.append(0)
                elif in_fault_zone and 1200.0 < center_x <= 3600.0:
                    attributes.append(1)
                else:
                    attributes.append(2)
                young_modulus.append(1.495e10)
                poisson_ratio.append(0.0)

    add_array(domain.GetPointData(), "GLOBAL_IDS_POINTS", np.arange(domain.GetNumberOfPoints()), vtk.VTK_LONG_LONG)
    add_array(domain.GetCellData(), "attribute", attributes, vtk.VTK_INT)
    add_array(domain.GetCellData(), "YoungModulus", young_modulus)
    add_array(domain.GetCellData(), "PoissonRatio", poisson_ratio)
    add_array(domain.GetCellData(), "GLOBAL_IDS_CELLS", np.arange(domain.GetNumberOfCells()), vtk.VTK_LONG_LONG)
    domain.GetPointData().SetGlobalIds(domain.GetPointData().GetArray("GLOBAL_IDS_POINTS"))
    domain.GetCellData().SetGlobalIds(domain.GetCellData().GetArray("GLOBAL_IDS_CELLS"))

    fracture_points = vtk.vtkPoints()
    collocated = []
    fracture_point_ids = {}
    for iz in fault_z:
        for iy in fault_y:
            fracture_point_ids[iy, iz] = fracture_points.InsertNextPoint(1200.0, y_values[iy], z_values[iz])
            collocated.append((base_ids[fault_x, iy, iz], duplicate_ids.get((iy, iz), -1)))

    fracture = vtk.vtkUnstructuredGrid()
    fracture.SetPoints(fracture_points)
    for iz in list(fault_z)[:-1]:
        for iy in list(fault_y)[:-1]:
            ids = [
                fracture_point_ids[iy, iz],
                fracture_point_ids[iy + 1, iz],
                fracture_point_ids[iy + 1, iz + 1],
                fracture_point_ids[iy, iz + 1],
            ]
            fracture.InsertNextCell(vtk.VTK_POLYGON, 4, ids)
    add_array(fracture.GetPointData(), "collocated_nodes", collocated, vtk.VTK_LONG_LONG)
    add_array(fracture.GetPointData(), "GLOBAL_IDS_POINTS", np.arange(fracture.GetNumberOfPoints()), vtk.VTK_LONG_LONG)
    add_array(fracture.GetCellData(), "GLOBAL_IDS_CELLS", np.arange(fracture.GetNumberOfCells()), vtk.VTK_LONG_LONG)
    fracture.GetPointData().SetGlobalIds(fracture.GetPointData().GetArray("GLOBAL_IDS_POINTS"))
    fracture.GetCellData().SetGlobalIds(fracture.GetCellData().GetArray("GLOBAL_IDS_CELLS"))

    write_grid(domain, mesh_dir / "Domain_Mesh_FaultModel_smoke.vtu")
    write_grid(fracture, mesh_dir / "Fault_Mesh_FaultModel_smoke.vtu")


def make_curved_fracture_mesh() -> None:
    mesh_dir = ROOT / "inputFiles/poromechanicsFractures/MESH"
    original = read_grid(mesh_dir / "curvedFrac_domain.vtu")
    y_samples = list(range(0, 31, 3))
    z_samples = list(range(0, 81, 4))

    def original_base_id(ix, iy, iz):
        return (iz * 31 + iy) * 6 + ix

    points = vtk.vtkPoints()
    base_ids = {}
    duplicate_ids = {}
    for new_z, old_z in enumerate(z_samples):
        for new_y, old_y in enumerate(y_samples):
            for ix in range(6):
                base_ids[ix, new_y, new_z] = points.InsertNextPoint(original.GetPoint(original_base_id(ix, old_y, old_z)))
            for ix in (2, 3):
                duplicate_ids[ix, new_y, new_z] = points.InsertNextPoint(original.GetPoint(original_base_id(ix, old_y, old_z)))

    domain = vtk.vtkUnstructuredGrid()
    domain.SetPoints(points)
    attributes = []
    reduced_centers = []
    for new_z in range(len(z_samples) - 1):
        old_z_cell = (z_samples[new_z] + z_samples[new_z + 1] - 1) // 2
        for new_y in range(len(y_samples) - 1):
            old_y_cell = (y_samples[new_y] + y_samples[new_y + 1] - 1) // 2
            for ix in range(5):
                left = duplicate_ids[2, new_y, new_z] if ix == 2 else base_ids[ix, new_y, new_z]
                right = duplicate_ids[3, new_y, new_z] if ix == 2 else base_ids[ix + 1, new_y, new_z]
                left_y = duplicate_ids[2, new_y + 1, new_z] if ix == 2 else base_ids[ix, new_y + 1, new_z]
                right_y = duplicate_ids[3, new_y + 1, new_z] if ix == 2 else base_ids[ix + 1, new_y + 1, new_z]
                left_z = duplicate_ids[2, new_y, new_z + 1] if ix == 2 else base_ids[ix, new_y, new_z + 1]
                right_z = duplicate_ids[3, new_y, new_z + 1] if ix == 2 else base_ids[ix + 1, new_y, new_z + 1]
                left_yz = duplicate_ids[2, new_y + 1, new_z + 1] if ix == 2 else base_ids[ix, new_y + 1, new_z + 1]
                right_yz = duplicate_ids[3, new_y + 1, new_z + 1] if ix == 2 else base_ids[ix + 1, new_y + 1, new_z + 1]
                ids = [left, right, right_y, left_y, left_z, right_z, right_yz, left_yz]
                domain.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, ids)
                reduced_centers.append(np.mean([points.GetPoint(point_id) for point_id in ids], axis=0))
                old_cell = (old_z_cell * 30 + old_y_cell) * 5 + ix
                attributes.append(int(original.GetCellData().GetArray("attribute").GetTuple1(old_cell)))

    # Thin material bands can fall between the coarse sampling planes. Preserve
    # every original 3-D region by assigning a nearby coarse cell to any missing
    # attribute, which keeps the smoke deck's Region/Reservoir selectors valid.
    original_centers_filter = vtk.vtkCellCenters()
    original_centers_filter.SetInputData(original)
    original_centers_filter.Update()
    original_centers = original_centers_filter.GetOutput()
    original_attributes = original.GetCellData().GetArray("attribute")
    reduced_centers = np.asarray(reduced_centers)
    counts = {attribute: attributes.count(attribute) for attribute in range(1, 16)}
    reserved = set()
    for required_attribute in range(1, 16):
        if counts[required_attribute]:
            continue
        source_points = [
            original_centers.GetPoint(cell_id)
            for cell_id in range(original.GetNumberOfCells())
            if int(original_attributes.GetTuple1(cell_id)) == required_attribute
        ]
        target = np.mean(source_points, axis=0)
        distances = np.linalg.norm(reduced_centers - target, axis=1)
        for cell_id in np.argsort(distances):
            old_attribute = attributes[cell_id]
            if cell_id not in reserved and counts[old_attribute] > 1:
                counts[old_attribute] -= 1
                attributes[cell_id] = required_attribute
                counts[required_attribute] += 1
                reserved.add(cell_id)
                break

    for new_z in range(len(z_samples) - 1):
        for new_y in range(len(y_samples) - 1):
            for ix, attribute in ((2, 21), (3, 22)):
                ids = [
                    base_ids[ix, new_y, new_z],
                    base_ids[ix, new_y + 1, new_z],
                    base_ids[ix, new_y + 1, new_z + 1],
                    base_ids[ix, new_y, new_z + 1],
                ]
                domain.InsertNextCell(vtk.VTK_QUAD, 4, ids)
                attributes.append(attribute)

    add_array(domain.GetPointData(), "GLOBAL_IDS_POINTS", np.arange(domain.GetNumberOfPoints()), vtk.VTK_LONG_LONG)
    add_array(domain.GetCellData(), "attribute", attributes, vtk.VTK_INT)
    add_array(domain.GetCellData(), "GLOBAL_IDS_CELLS", np.arange(domain.GetNumberOfCells()), vtk.VTK_LONG_LONG)
    domain.GetPointData().SetGlobalIds(domain.GetPointData().GetArray("GLOBAL_IDS_POINTS"))
    domain.GetCellData().SetGlobalIds(domain.GetCellData().GetArray("GLOBAL_IDS_CELLS"))

    fracture_points = vtk.vtkPoints()
    fracture_ids = {}
    collocated = []
    for new_z in range(len(z_samples)):
        for new_y in range(len(y_samples)):
            for ix in (2, 3):
                fracture_ids[ix, new_y, new_z] = fracture_points.InsertNextPoint(
                    original.GetPoint(original_base_id(ix, y_samples[new_y], z_samples[new_z]))
                )
                if ix == 2:
                    collocated.append((duplicate_ids[ix, new_y, new_z], base_ids[ix, new_y, new_z]))
                else:
                    collocated.append((base_ids[ix, new_y, new_z], duplicate_ids[ix, new_y, new_z]))

    fracture = vtk.vtkUnstructuredGrid()
    fracture.SetPoints(fracture_points)
    fracture_attributes = []
    for new_z in range(len(z_samples) - 1):
        for new_y in range(len(y_samples) - 1):
            for ix, attribute in ((2, 21), (3, 22)):
                ids = [
                    fracture_ids[ix, new_y, new_z],
                    fracture_ids[ix, new_y + 1, new_z],
                    fracture_ids[ix, new_y + 1, new_z + 1],
                    fracture_ids[ix, new_y, new_z + 1],
                ]
                fracture.InsertNextCell(vtk.VTK_POLYGON, 4, ids)
                fracture_attributes.append(attribute)
    add_array(fracture.GetPointData(), "collocated_nodes", collocated, vtk.VTK_LONG_LONG)
    add_array(fracture.GetPointData(), "GLOBAL_IDS_POINTS", np.arange(fracture.GetNumberOfPoints()), vtk.VTK_LONG_LONG)
    add_array(fracture.GetCellData(), "attribute", fracture_attributes, vtk.VTK_INT)
    add_array(fracture.GetCellData(), "GLOBAL_IDS_CELLS", np.arange(fracture.GetNumberOfCells()), vtk.VTK_LONG_LONG)
    fracture.GetPointData().SetGlobalIds(fracture.GetPointData().GetArray("GLOBAL_IDS_POINTS"))
    fracture.GetCellData().SetGlobalIds(fracture.GetCellData().GetArray("GLOBAL_IDS_CELLS"))

    write_grid(domain, mesh_dir / "curvedFrac_domain_smoke.vtu")
    write_grid(fracture, mesh_dir / "curvedFrac_fracture_smoke.vtu")


if __name__ == "__main__":
    make_hybrid_property_mesh()
    make_fault_mesh()
    make_curved_fracture_mesh()
