from collections import defaultdict
from dataclasses import dataclass
from glob import glob
import logging
import sys
from typing import (
    Collection,
    Dict,
    Optional,
    FrozenSet,
    Mapping
)
from copy import deepcopy

from vtkmodules.vtkCommonDataModel import (
    vtkUnstructuredGrid
)
from vtkmodules.vtkIOXML import (
    vtkXMLUnstructuredGridReader,
)
from vtkmodules.util.numpy_support import (
    vtk_to_numpy,
    numpy_to_vtk,
)
from vtk import (vtkCell,
                 vtkBoundingBox)

logger = logging.getLogger("ghosting")


def __read_vtu(vtk_input_file: str) -> Optional[vtkUnstructuredGrid]:
    reader = vtkXMLUnstructuredGridReader()
    logger.info(f"Testing file format \"{vtk_input_file}\" using XML format reader...")
    if reader.CanReadFile(vtk_input_file):
        reader.SetFileName(vtk_input_file)
        logger.info(f"Reader matches. Reading file \"{vtk_input_file}\" using XML format reader.")
        reader.Update()
        return reader.GetOutput()
    else:
        logger.info("Reader did not match the input file format.")
        return None


@dataclass(frozen=True)
class Node:
    local: int  # vtk index
    global_: int  # vtk global index


@dataclass(frozen=True)
class Edge:
    nodes: tuple[int, int]  # Local indices in the vtk mesh. Sorted.


@dataclass(frozen=True)
class MeshGraph:
    nodes: Collection[Node]
    edges: Collection[Edge]


def build_edges(mesh: vtkUnstructuredGrid, points_gids) -> Collection[Edge]:
    tmp: set[Edge] = set()
    for c in range(mesh.GetNumberOfCells()):
        cell: vtkCell = mesh.GetCell(c)
        for e in range(cell.GetNumberOfEdges()):
            edge: vtkCell = cell.GetEdge(e)
            ls = edge.GetPointId(0), edge.GetPointId(1)
            # g0, g1 = points_gids[l0], points_gids[l1]
            # tmp.add(tuple(sorted((g0, g1))))
            nodes = tuple(sorted(ls))
            tmp.add(Edge(nodes))
    return tuple(tmp)


def compute_graph(mesh: vtkUnstructuredGrid) -> MeshGraph:
    num_points: int = mesh.GetNumberOfPoints()
    points_gids = vtk_to_numpy(mesh.GetPointData().GetGlobalIds())
    nodes = []
    for l, g in enumerate(points_gids):
        nodes.append(Node(l, g))
    edges: Collection[Edge] = build_edges(mesh, points_gids)
    graph: MeshGraph = MeshGraph(tuple(nodes), edges)
    return graph


def find_overlapping_edges(graphs: Collection[MeshGraph], neighbors: Collection[Collection[int]]):
    # [rank] -> [ [edge sorted global nodes] -> [edge index in the MeshGraph] ]
    tmp: dict[int, dict[tuple[int, int], int]] = dict()
    for rank, graph in enumerate(graphs):
        d: dict[tuple[int, int], int] = {}
        for ie, edge in enumerate(graph.edges):
            gn0: int = graph.nodes[edge.nodes[0]].global_
            gn1: int = graph.nodes[edge.nodes[1]].global_
            gns: tuple[int, int] = tuple(sorted((gn0, gn1)))
            d[gns] = ie
        tmp[rank] = d

    # # Build a global to local map for the nodes
    # ng2l: dict[int, tuple[int, int]] = {}  # local is in fact (rank, local)
    # for rank, graph in enumerate(graphs):
    #     for node in graph.nodes:
    #         ng2l[node.global_] = (rank, node.local)

    # raise ValueError("TODO")  # Maybe `intersections` should be done neighborhood by neighborhood, so is `count`, as a tmp.
    # [rank_a, rank_b] -> [Collection[Edge]]
    # intersections: dict[tuple[int, ...], set[tuple[int, int]]] = defaultdict(set)

    rank_to_intersections: list[dict[tuple[int, ...], set[tuple[int, int]]]] = []
    for i in range(len(graphs)):
        rank_to_intersections.append(defaultdict(set))

    for current_rank, intersection in enumerate(rank_to_intersections):
        count: dict[tuple[int, int], set[int]] = defaultdict(set)
        other_ranks = neighbors[current_rank]
        # for current_rank, other_ranks in enumerate(neighbors):
        all_ranks = [current_rank, ] + other_ranks
        for rank in all_ranks:
            for edge in tmp[rank]:
                count[edge].add(rank)
        for edge, ranks in count.items():
            if len(ranks) > 0 and current_rank in ranks:
                intersection[tuple(sorted(ranks))].add(edge)

    # Check if neighborhood is too wide...
    for current_rank, other_ranks in enumerate(neighbors):
        i = rank_to_intersections[current_rank]
        useful_neighbors = set()
        for ranks, edges in i.items():
            if edges:
                useful_neighbors |= set(ranks)
        print(f"Ranks '{set([current_rank, ] + other_ranks) - useful_neighbors}' are not required for rank {current_rank}.")

    # Now let's build the scan (or equivalent)
    scans: list[Mapping[tuple[int, ...], int]] = []
    # for rank in range(len(graphs)):  # PUT BACK
    for rank, intersections in enumerate(rank_to_intersections):
        scan = {} if len(scans) == 0 else deepcopy(scans[-1])
        for ranks, edges in intersections.items():
            # I want the current rank to have its own information for easy access to the offset,
            # so there's no special case between what the local information and the information from the other ranks.
            # Hence the `+ 1` so it's included in the `range`
            if set(range(rank + 1)) & set(ranks):
                scan[ranks] = len(edges)
        scans.append(scan)

    offset_scans: list[Mapping[tuple[int, ...], int]] = []
    for scan in scans:
        offset_scan: dict[tuple[int, ...], int] = dict()
        i: int = 0
        for ranks in sorted(scan.keys()):
            n: int = scan.get(ranks)
            offset_scan[ranks] = i
            i += n
        offset_scans.append(offset_scan)

    # From `intersections`, let's mimic an intersection for the rank only
    real_intersection: list[Mapping[tuple[int, ...], FrozenSet[tuple[int, int]]]] = []
    # for rank in range(len(graphs)):  # PUT BACK
    for rank, intersections in enumerate(rank_to_intersections):
        i: Dict[tuple[int, ...], set[tuple[int, int]]] = dict()
        for ranks, edges in intersections.items():
            if rank in ranks:
                i[ranks] = edges
        real_intersection.append(i)

    return real_intersection, offset_scans


def do_the_numbering(intersection: Mapping[tuple[int, ...], FrozenSet[tuple[int, int]]],
                     offset_scan: Mapping[tuple[int, ...], int]) -> Mapping[tuple[int, int], int]:
    numbered_edges: dict[tuple[int, int], int] = dict()
    for ranks, edges in sorted(intersection.items()):
        off: int = offset_scan[ranks]
        for edge in sorted(edges):
            numbered_edges[edge] = off
            off += 1
    return numbered_edges


def validation(numberings: Mapping[tuple[int, int], int]) -> int:
    all_nodes = set()
    all_edges = dict()
    for numbering in numberings:
        for nodes, gei in numbering.items():
            all_nodes |= set(nodes)
            res = all_edges.get(nodes)
            if res is None:
                all_edges[nodes] = gei
            elif res != nodes:
                assert res == gei
    expected_num_edges: int = 3 * 10 * 11 * 11
    expected_num_nodes: int = 11 * 11 * 11
    assert len(all_edges) == expected_num_edges
    assert set(all_edges.values()) == set(range(expected_num_edges))
    assert all_nodes == set(range(expected_num_nodes))
    return 0


def build_neighborhood(meshes: Collection[vtkUnstructuredGrid]) -> Collection[Collection[int]]:
    bounding_boxes: list[vtkBoundingBox] = []
    for mesh in meshes:
        bb = vtkBoundingBox()
        bb.AddBox(mesh.GetBounds())
        bb.Inflate(0.1)
        bounding_boxes.append(bb)
    neighbors: list[list[int]] = []
    for rank, bb in enumerate(bounding_boxes):
        n = []
        for other_rank, other_bb in enumerate(bounding_boxes):
            if rank == other_rank:
                continue
            if bb.Intersects(other_bb):
                n.append(other_rank)
        neighbors.append(n)
    return neighbors


def main() -> int:
    # For each rank, contains the raw vtk mesh.
    meshes: list[vtkUnstructuredGrid] = []
    pattern: str = "/Users/j0436735/Downloads/meshes-cube/main-*.vtu"
    # pattern: str = "/Users/j0436735/Downloads/meshes-cube-25/main-*.vtu"
    for file_name in glob(pattern):
        m = __read_vtu(file_name)
        if m:
            meshes.append(m)

    # For each rank, contains the neighbor rank (candidates).
    neighbors: Collection[Collection[int]] = build_neighborhood(meshes)

    # For each rank, contains the graph built upon the vtk mesh.
    graphs = []
    for m in meshes:
        graphs.append(compute_graph(m))

    # Perform the core computation of the intersection.
    # `intersections` will contain the intersection for each rank,
    # while `offset_scans` will contains the offset for the global numbering,
    # as it would be done in during the `MPI_scan`.
    intersections: Collection[Mapping[tuple[int, ...], FrozenSet[tuple[int, int]]]]
    offset_scans: Collection[Mapping[tuple[int, ...], int]]
    intersections, offset_scans = find_overlapping_edges(graphs, neighbors)

    # Finish by doing the global numbering for real.
    numberings: list[Mapping[tuple[int, int], int]] = []
    for i, o in zip(intersections, offset_scans):
        numberings.append(do_the_numbering(i, o))

    # Small validation
    return validation(numberings)


if __name__ == '__main__':
    sys.exit(main())
