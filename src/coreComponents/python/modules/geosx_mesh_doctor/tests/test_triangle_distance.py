from dataclasses import dataclass

import numpy
from numpy.linalg import norm
import pytest

from checks.triangle_distance import distance_between_two_segments, distance_between_two_triangles


@dataclass(frozen=True)
class ExpectedSeg:
    p0: numpy.array
    u0: numpy.array
    p1: numpy.array
    u1: numpy.array
    x: numpy.array
    y: numpy.array

    @classmethod
    def from_tuples(cls, p0, u0, p1, u1, x, y):
        return cls(
            numpy.array(p0),
            numpy.array(u0),
            numpy.array(p1),
            numpy.array(u1),
            numpy.array(x),
            numpy.array(y)
        )


def __get_segments_references():
    # Node to node configuration.
    yield ExpectedSeg.from_tuples(
        p0=(0., 0., 0.),
        u0=(1., 0., 0.),
        p1=(2., 0., 0.),
        u1=(1., 0., 0.),
        x=(1., 0., 0.),
        y=(2., 0., 0.),
    )
    # Node to edge configuration.
    yield ExpectedSeg.from_tuples(
        p0=(0., 0., 0.),
        u0=(1., 0., 0.),
        p1=(2., -1., -1.),
        u1=(0., 1., 1.),
        x=(1., 0., 0.),
        y=(2., 0., 0.),
    )
    # Edge to edge configuration.
    yield ExpectedSeg.from_tuples(
        p0=(0., 0., -1.),
        u0=(0., 0., 2.),
        p1=(1., -1., -1.),
        u1=(0., 2., 2.),
        x=(0., 0., 0.),
        y=(1., 0., 0.),
    )
    # Example from "On fast computation of distance between line segments" by Vladimir J. Lumelsky.
    # Information Processing Letters, Vol. 21, number 2, pages 55-61, 08/16/1985.
    # It's a node to edge configuration.
    yield ExpectedSeg.from_tuples(
        p0=(0., 0., 0.),
        u0=(1., 2., 1.),
        p1=(1., 0., 0.),
        u1=(1., 1., 0.),
        x=(1./6., 2./6., 1./6.),
        y=(1., 0., 0.),
    )
    # Overlapping edges.
    yield ExpectedSeg.from_tuples(
        p0=(0., 0., 0.),
        u0=(2., 0., 0.),
        p1=(1., 0., 0.),
        u1=(2., 0., 0.),
        x=(0., 0., 0.),
        y=(0., 0., 0.),
    )
    # Crossing edges.
    yield ExpectedSeg.from_tuples(
        p0=(0., 0., 0.),
        u0=(2., 0., 0.),
        p1=(1., -1., 0.),
        u1=(0., 2., 0.),
        x=(0., 0., 0.),
        y=(0., 0., 0.),
    )


@pytest.mark.parametrize("expected", __get_segments_references())
def test_segments(expected: ExpectedSeg):
    eps = numpy.finfo(float).eps
    x, y = distance_between_two_segments(expected.p0, expected.u0, expected.p1, expected.u1)
    if norm(expected.x - expected.y) == 0:
        assert norm(x - y) == 0.
    else:
        assert norm(expected.x - x) < eps
        assert norm(expected.y - y) < eps


@dataclass(frozen=True)
class ExpectedTri:
    t0: numpy.array
    t1: numpy.array
    d: float
    p0: numpy.array
    p1: numpy.array

    @classmethod
    def from_tuples(cls, t0, t1, d, p0, p1):
        return cls(
            numpy.array(t0),
            numpy.array(t1),
            float(d),
            numpy.array(p0),
            numpy.array(p1)
        )


def __get_triangles_references():
    # Node to node configuration.
    yield ExpectedTri.from_tuples(
        t0=((0., 0., 0.), (1., 0., 0.), (0., 1., 1.)),
        t1=((2., 0., 0.), (3., 0., 0.), (2., 1., 1.)),
        d=1.,
        p0=(1., 0., 0.),
        p1=(2., 0., 0.)
    )
    # Node to edge configuration.
    yield ExpectedTri.from_tuples(
        t0=((0., 0., 0.), (1., 0., 0.), (0., 1., 1.)),
        t1=((2., -1., 0.), (3., 0., 0.), (2., 1., 0.)),
        d=1.,
        p0=(1., 0., 0.),
        p1=(2., 0., 0.)
    )
    # Edge to edge configuration.
    yield ExpectedTri.from_tuples(
        t0=((0., 0., 0.), (1., 1., 1.), (1., -1., -1.)),
        t1=((2., -1., 0.), (2., 1., 0.), (3., 0., 0.)),
        d=1.,
        p0=(1., 0., 0.),
        p1=(2., 0., 0.)
    )
    # Point to face configuration.
    yield ExpectedTri.from_tuples(
        t0=((0., 0., 0.), (1., 0., 0.), (0., 1., 1.)),
        t1=((2., -1., 0.), (2., 1., -1.), (2, 1., 1.)),
        d=1.,
        p0=(1., 0., 0.),
        p1=(2., 0., 0.)
    )
    # Same triangles configuration.
    yield ExpectedTri.from_tuples(
        t0=((0., 0., 0.), (1., 0., 0.), (0., 1., 1.)),
        t1=((0., 0., 0.), (1., 0., 0.), (0., 1., 1.)),
        d=0.,
        p0=(0., 0., 0.),
        p1=(0., 0., 0.)
    )
    # Crossing triangles configuration.
    yield ExpectedTri.from_tuples(
        t0=((0., 0., 0.), (2., 0., 0.), (2., 0., 1.)),
        t1=((1., -1., 0.), (1., 1., 0.), (1., 1., 1.)),
        d=0.,
        p0=(0., 0., 0.),
        p1=(0., 0., 0.)
    )


@pytest.mark.parametrize("expected", __get_triangles_references())
def test_triangles(expected: ExpectedTri):
    eps = numpy.finfo(float).eps
    d, p0, p1 = distance_between_two_triangles(expected.t0, expected.t1)
    assert abs(d - expected.d) < eps
    if d != 0:
        assert norm(p0 - expected.p0) < eps
        assert norm(p1 - expected.p1) < eps
