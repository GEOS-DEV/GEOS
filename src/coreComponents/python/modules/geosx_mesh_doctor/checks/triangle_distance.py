import itertools
from math import sqrt
from typing import Tuple, Union

import numpy
from numpy.linalg import norm


def __div_clamp(num: float, den :float) -> float:
    """
    Computes the division `num / den`. and clamps the result between 0 and 1.
    If `den` is zero, the result of the division is set to 0.
    :param num: The numerator.
    :param den: The denominator.
    :return: The result between 0 and 1.
    """
    if den == 0.:
        return 0.
    tmp: float = num / den
    if tmp < 0:
        return 0.
    elif tmp > 1:
        return 1.
    else:
        return tmp


def distance_between_two_segments(x0: numpy.ndarray, d0: numpy.ndarray,
                                  x1: numpy.ndarray, d1: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """
    Compute the minimum distance between two segments.
    :param x0: First point of segment 0.
    :param d0: Director vector such that x0 + d0 is the second point of segment 0.
    :param x1: First point of segment 1.
    :param d1: Director vector such that x1 + d1 is the second point of segment 1.
    :return: A tuple containing the two points closest point for segments 0 and 1 respectively.
    """
    # The reference paper is:
    # "On fast computation of distance between line segments" by Vladimir J. Lumelsky.
    # Information Processing Letters, Vol. 21, number 2, pages 55-61, 08/16/1985.

    # In the reference, the indices start at 1, while in this implementation, they start at 0.
    tmp: numpy.ndarray = x1 - x0
    D0: float = numpy.dot(d0, d0)  # As such, this is D1 in the reference paper.
    D1: float = numpy.dot(d1, d1)
    R: float = numpy.dot(d0, d1)
    S0: float = numpy.dot(d0, tmp)
    S1: float = numpy.dot(d1, tmp)

    # `t0` parameterizes line 0:
    #   - when t0 = 0 the point is p0.
    #   - when t0 = 1, the point is p0 + u0, the other side of the segment
    # Same for `t1` and line 1.

    # Step 1 of the algorithm considers degenerate cases.
    # They'll be considered along the line using `div_clamp`.

    # Step 2: Computing t0 using eq (11).
    t0: float = __div_clamp(S0 * D1 - S1 * R, D0 * D1 - R * R)

    # Step 3: compute t1 for point on line 1 closest to point at t0.
    t1: float = __div_clamp(t0 * R - S1, D1)  # Eq (10, right)
    sol_1: numpy.ndarray = x1 + t1 * d1            # Eq (3)
    t0: float = __div_clamp(t1 * R + S0, D0)  # Eq (10, left)
    sol_0: numpy.ndarray = x0 + t0 * d0            # Eq (4)

    return sol_0, sol_1


def __compute_nodes_to_triangle_distance(tri_0, edges_0, tri_1) -> Tuple[Union[float, None], Union[numpy.ndarray, None], Union[numpy.ndarray, None], bool]:
    """
    Computes the distance from nodes of `tri_1` points onto `tri_0`.
    :param tri_0: First triangle.
    :param edges_0: The edges of triangle 0. First element being edge [0, 1], etc.
    :param tri_1: Second triangle
    :return: The distance, the closest point on triangle 0, the closest on triangle 1
    and a boolean indicating of the triangles are disjoint. If nothing was found,
    then the first three arguments are None. The boolean being still defined.
    """
    are_disjoint: bool = False
    tri_0_normal: numpy.ndarray = numpy.cross(edges_0[0], edges_0[1])
    tri_0_normal_norm: float = numpy.dot(tri_0_normal, tri_0_normal)

    # Forget about degenerate cases.
    if tri_0_normal_norm > numpy.finfo(float).eps:
        # Build projection lengths of `tri_1` points.
        tri_1_proj = numpy.empty(3, dtype=float)
        for i in range(3):
            tri_1_proj[i] = numpy.dot(tri_0[0] - tri_1[i], tri_0_normal)

        # Considering `tri_0` separates the space in 2,
        # let's check if `tri_1` is on one side only.
        # If so, let's take the closest point.
        point: int = -1
        if numpy.all(tri_1_proj > 0):
            point = numpy.argmin(tri_1_proj)
        elif numpy.all(tri_1_proj < 0):
            point = numpy.argmax(tri_1_proj)

        # So if `tri_1` is actually "on one side",
        # point `tri_1[point]` is candidate to be the closest point.
        if point > -1:
            are_disjoint = True
            # But we must check that its projection is inside `tri_0`.
            if numpy.dot(tri_1[point] - tri_0[0], numpy.cross(tri_0_normal, edges_0[0])) > 0:
                if numpy.dot(tri_1[point] - tri_0[1], numpy.cross(tri_0_normal, edges_0[1])) > 0:
                    if numpy.dot(tri_1[point] - tri_0[2], numpy.cross(tri_0_normal, edges_0[2])) > 0:
                        # It is!
                        sol_0 = tri_1[point]
                        sol_1 = tri_1[point] + (tri_1_proj[point] / tri_0_normal_norm) * tri_0_normal
                        return norm(sol_1 - sol_0), sol_0, sol_1, are_disjoint
    return None, None, None, are_disjoint


def distance_between_two_triangles(tri_0: numpy.ndarray,
                                   tri_1: numpy.ndarray) -> Tuple[float, numpy.ndarray, numpy.ndarray]:
    """
    Returns the minimum distance between two triangles, and the two points where this minimum occurs.
    If the two triangles touch, then distance is exactly 0.
    But the two points are dummy and cannot be used as contact points (they are still though).
    :param tri_0: The first 3x3 triangle points.
    :param tri_1: The second 3x3 triangle points.
    :return: The distance and the two points.
    """
    # Compute vectors along the 6 sides
    edges_0 = numpy.empty((3, 3), dtype=float)
    edges_1 = numpy.empty((3, 3), dtype=float)
    for i in range(3):
        edges_0[i][:] = tri_0[(i + 1) % 3] - tri_0[i]
        edges_1[i][:] = tri_1[(i + 1) % 3] - tri_1[i]

    min_sol_0 = numpy.empty(3, dtype=float)
    min_sol_1 = numpy.empty(3, dtype=float)
    are_disjoint: bool = False

    min_dist = numpy.inf

    # Looping over all the pair of edges.
    for i, j in itertools.product(range(3), repeat=2):
        # Find the closest points on edges i and j.
        sol_0, sol_1 = distance_between_two_segments(tri_0[i], edges_0[i], tri_1[j], edges_1[j])
        # Computing the distance between the two solutions.
        delta_sol = sol_1 - sol_0
        dist: float = numpy.dot(delta_sol, delta_sol)
        # Update minimum if relevant and check if it's the closest pair of points.
        if dist <= min_dist:
            min_sol_0[:] = sol_0
            min_sol_1[:] = sol_1
            min_dist = dist

            # `tri_0[(i + 2) % 3]` is the points opposite to edges_0[i] where the closest point sol_0 lies.
            # Computing those scalar products and checking the signs somehow let us determine
            # if the triangles are getting closer to each other when approaching the sol_(0|1) nodes.
            # If so, we have a minimum.
            a: float = numpy.dot(tri_0[(i + 2) % 3] - sol_0, delta_sol)
            b: float = numpy.dot(tri_1[(j + 2) % 3] - sol_1, delta_sol)
            if a <= 0 <= b:
                return sqrt(dist), sol_0, sol_1

            if a < 0:
                a = 0
            if b > 0:
                b = 0
            # `dist - a + b` expands to `numpy.dot(tri_1[(j + 2) % 3] - tri_0[(i + 2) % 3], sol_1 - sol_0)`.
            # If the "directions" of the (sol_1 - sol_0) vector and the vector joining the extra points of the triangles
            # (i.e. not involved in the current edge check) re the "same", then the triangles do not intersect.
            if dist - a + b > 0:
                are_disjoint = True
    # No edge pair contained the closest points.
    # Checking the node/face situation.
    distance, sol_0, sol_1, are_disjoint_tmp = __compute_nodes_to_triangle_distance(tri_0, edges_0, tri_1)
    if distance:
        return distance, sol_0, sol_1
    are_disjoint = are_disjoint or are_disjoint_tmp

    distance, sol_0, sol_1, are_disjoint_tmp = __compute_nodes_to_triangle_distance(tri_1, edges_1, tri_0)
    if distance:
        return distance, sol_0, sol_1
    are_disjoint = are_disjoint or are_disjoint_tmp
    # It's not a node/face situation.
    # If the triangles do not overlap, let's return the minimum found during the edges loop.
    # (maybe an edge was parallel to the other face, and we could not decide for a unique closest point).
    if are_disjoint:
        return sqrt(min_dist), min_sol_0, min_sol_1
    else:  # Surely overlapping or degenerate triangles.
        return 0., numpy.zeros(3, dtype=float), numpy.zeros(3, dtype=float)
