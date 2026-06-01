# -*- coding: utf-8 -*-
"""Convenience helpers for GEOS-MPM Lagrangian tracer-point input.

The helpers in this file generate the SolidMechanics_MPM XML attributes used by
lagrangian particle tracers.  They are intentionally small and dependency-free
so that they can be imported directly from a PFW input file.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import math
from collections.abc import Iterable


_AXIS_INDEX = {"x": 0, "y": 1, "z": 2, 0: 0, 1: 1, 2: 2}


def _axis_index(axis):
  try:
    return _AXIS_INDEX[axis]
  except KeyError:
    raise ValueError("axis must be 'x', 'y', 'z', 0, 1, or 2")


def _as_values(value):
  """Return value as a list, treating strings as scalar values."""
  if isinstance(value, str):
    return [value]
  if isinstance(value, Iterable):
    return list(value)
  return [value]


def _format_float(value):
  """Format a numeric value for GEOS XML input."""
  return "{:.17g}".format(float(value))


def linspace(start, stop, count):
  """Return count evenly spaced points from start to stop, inclusive."""
  if count <= 0:
    raise ValueError("count must be positive")
  if count == 1:
    return [float(start)]
  start = float(start)
  stop = float(stop)
  step = (stop - start) / float(count - 1)
  return [start + i * step for i in range(count)]


def tensor_grid(x, y, z):
  """Create a tensor-product list of (x, y, z) tracer coordinates.

  Each argument can be either a scalar or a sequence.  For example,
  tensor_grid(x=x_face, y=linspace(-r, r, 10), z=linspace(-r, r, 10))
  creates a 10x10 plane of tracer points normal to x.
  """
  points = []
  for xi in _as_values(x):
    for yi in _as_values(y):
      for zi in _as_values(z):
        points.append((float(xi), float(yi), float(zi)))
  return points


def plane_grid(center, axis1, axis2, span1, span2, count1, count2):
  """Create a rectangular plane of tracer points around center.

  axis1 and axis2 select the two in-plane coordinate axes.  The remaining axis
  is held fixed at center.  span1/span2 are the total widths along the selected
  axes.
  """
  c = [float(v) for v in center]
  if len(c) != 3:
    raise ValueError("center must contain exactly three coordinates")

  i1 = _axis_index(axis1)
  i2 = _axis_index(axis2)
  if i1 == i2:
    raise ValueError("axis1 and axis2 must be different axes")

  values1 = linspace(c[i1] - 0.5 * float(span1), c[i1] + 0.5 * float(span1), count1)
  values2 = linspace(c[i2] - 0.5 * float(span2), c[i2] + 0.5 * float(span2), count2)

  points = []
  for v1 in values1:
    for v2 in values2:
      p = list(c)
      p[i1] = v1
      p[i2] = v2
      points.append(tuple(p))
  return points


def disk(center, radius, normal_axis="x", radial_count=5, angular_count=32, include_center=True):
  """Create a disk of tracer points in the plane normal to normal_axis.

  radial_count is the number of nonzero radial rings.  angular_count points are
  placed on each ring.  The returned point count is radial_count*angular_count,
  plus one if include_center is True.
  """
  c = [float(v) for v in center]
  if len(c) != 3:
    raise ValueError("center must contain exactly three coordinates")
  if radius < 0.0:
    raise ValueError("radius must be non-negative")
  if radial_count < 0:
    raise ValueError("radial_count must be non-negative")
  if angular_count <= 0:
    raise ValueError("angular_count must be positive")

  normal = _axis_index(normal_axis)
  in_plane_axes = [axis for axis in range(3) if axis != normal]

  points = []
  if include_center:
    points.append(tuple(c))

  if radial_count == 0 or radius == 0.0:
    return points

  for radial_index in range(1, radial_count + 1):
    r = float(radius) * radial_index / float(radial_count)
    for angular_index in range(angular_count):
      theta = 2.0 * math.pi * angular_index / float(angular_count)
      p = list(c)
      p[in_plane_axes[0]] += r * math.cos(theta)
      p[in_plane_axes[1]] += r * math.sin(theta)
      points.append(tuple(p))

  return points


def format_coordinates(points):
  """Format tracer coordinates as a GEOS array2d XML attribute string."""
  rows = []
  for point in points:
    p = list(point)
    if len(p) != 3:
      raise ValueError("each tracer point must contain exactly three coordinates")
    rows.append("{ " + ", ".join(_format_float(value) for value in p) + " }")

  if not rows:
    raise ValueError("at least one tracer point is required")

  return "{ " + ", ".join(rows) + " }"


def format_variables(variables):
  """Format tracer variable names as a GEOS string_array XML attribute."""
  values = [str(variable) for variable in _as_values(variables)]
  if not values:
    raise ValueError("at least one tracer variable is required")
  return "{ " + ", ".join(values) + " }"


def set_tracers(pfw, points, variables, write_interval=None, cycle_interval=None, output_prefix="tracerPoint"):
  """Populate a PFW dictionary with lagrangian tracer-point settings.

  Use either write_interval for time-based output or cycle_interval for
  cycle-based output, not both.  The input points can be generated with
  tensor_grid(), plane_grid(), disk(), or supplied as an iterable of triples.
  """
  if write_interval is not None and cycle_interval is not None:
    raise ValueError("use either write_interval or cycle_interval, not both")

  pfw["tracerHistory"] = 1
  pfw["tracerCoordinates"] = format_coordinates(points)
  pfw["tracerVariables"] = format_variables(variables)
  pfw["tracerOutputPrefix"] = str(output_prefix)

  if write_interval is not None:
    pfw["tracerWriteInterval"] = write_interval
  if cycle_interval is not None:
    pfw["tracerCycleInterval"] = cycle_interval

  return pfw
