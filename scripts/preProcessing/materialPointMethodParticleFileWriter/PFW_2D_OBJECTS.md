# PFW objects for native plane strain

Native XY geometry should not require an artificial extrusion to keep end caps
away from particles. The writer now prepares supported legacy extrusions before
querying occupancy, flags, normals, or surface offsets. Physical thickness is the
single domain scalar `planeStrainThickness`, independent of object geometry.

## Audit of the supplied plane-strain examples

| Examples | Previous object use | Native treatment |
|---|---|---|
| Brazilian disk: B-spline, FLIP, FMPM, PIC, XPIC | z-axis cylinder inside a 3D Voronoi-Weibull wrapper | disk plus dim=2 Voronoi-Weibull; sample volume=cell area*t |
| Borehole collapse | z-axis annular cylinder and 3D Voronoi-Weibull wrapper | disk with ri for the hole, dim=2 Weibull |
| Colliding disks; elastic disk | z-axis cylinders | disk objects with XY centers and no end caps |
| Copper foam compression | foam(dim=2) and polygon platens | Already query in XY; retained and tested |
| Mixed powder compression | packedSphericalBed(dim=2), 2D box/polygon fixtures | Already query in XY; retained and tested |
| Rigid-body jamming | polygons and cylinders with custom disk-normal overrides | Native disks and existing polygons; obsolete normal override removed |

All 11 inputs select PlaneStrain and specify one physical thickness. None needs
zmin/zmax/nK for mesh construction. Existing default thickness values are retained
for comparison, with comments explaining that physical t should be chosen
independently of grid resolution. The new 2D Weibull realizations differ from the
old 3D random fields even with the same seed.

## Native constructors

```python
solid = geom.disk('disk', center=[0.0, 0.0], r=1.0,
                  particleType=2, mat=0, group=0)
ring = geom.disk('ring', center=[0.0, 0.0], r=1.0, ri=0.4,
                 flaggedSurfaces=[True, True])  # outer, inner
block = geom.rectangle('block', x0=[-1.0, -0.5], x1=[1.0, 0.5],
                        flaggedSurfaces=[True, True, True, True])
# rectangle is a convenience constructor for box(..., dim=2).
```

Disk radius must be positive and 0 <= ri < r. The disk center belongs to a solid
disk. Annular normals point outward on the outer circle and into the hole on the
inner circle. Flags select the eligible boundaries for both normals and offsets;
no enabled boundary returns zero surface data. Rectangle surface order follows
box(dim=2): x-min, y-min, x-max, y-max. These are object flags, not solver BC order.
Returned positions/velocities/material directions retain the established
three-component interfaces; geometry queries ignore z. The disk material frame
matches the positive-z legacy cylinder frame, so migrating geometry does not
silently rotate material directions.

## Old PFW inputs

For PlaneStrain only, the preparation traversal adapts:

- Axis-aligned box(dim=3) to its XY rectangle, dropping z caps and preserving
  the x/y flags (old indices 0,1,3,4).
- z-axis cylinders with transverse cap planes to disk/annulus sections. Cylinder
  length, cap flags, and query z no longer affect the section's surface data.
- Owned subobjects of supported wrappers, composites, and transformations that
  preserve the XY plane. Outer material/group/velocity/particle-type attributes
  and property methods remain on the original wrapper.

Consequently, an old extruded disk input works even if its cylinder is very thin;
there is no need to make its depth larger than the radial surface distance.
The previous cylinder normal/offset calculation could select a nearby end cap
even when the cap was not the desired flagged boundary. The native section uses
only its in-plane eligible boundaries.

This adaptation assumes these simple shapes represent intended plane-strain
extrusions. A disjoint box/cylinder positioned away from the old query z plane
is also treated as its XY extrusion. To model an intentional finite 3D slice,
use an explicit cross-section geometry/custom object instead of relying on caps.
3D runs retain their existing object behavior.

Tilted cylinders, oblique end caps, spheres, STL/general 3D objects, and transforms
mixing XY with z retain existing 3D cross-section sampling. They are not
silently converted to disks or polygons. Their surface methods still need to
describe the intended XY section; simply projecting a 3D normal does not make
arbitrary 3D cutting geometry exact. Prefer native polygons/disks or an explicit
section implementation for these cases. Existing 3D stochastic wrappers are
not silently converted to 2D statistical models.

## Thickness and Weibull fields

Native dim=2 versions of voronoiWeibullBoxWrapper,
layeredStrengthScaleBoxWrapper, voronoiWrapperWithLayeredStrengthScale, and
layeredVoronoiWeibullWrapper use physical volume max(A*t,vMin) in the strength
law. The writer binds t after object creation, retaining random draws and
material frames. Rebinding is idempotent. An explicit legacy layered-wrapper
thickness must agree with the domain scalar. sizeEffectVoronoiWeibullBoxWrapper
and custom 3D flaw-file models are not converted to native 2D by this patch.

Built-in scalar/per-particle fields outside these wrappers retain their stated
units. A custom wrapper can implement setPlaneStrainThickness(t); arbitrary
user-provided StrengthScale values are left as supplied.

## Verification limits

Object tests compare membership, flags, normals and offsets with query z ranging
from -1e6 to +1e6; extrusion-adapter tests vary legacy depth from 1e-8 to 100.
Native Weibull tests exercise all four constructors, fixed random draws,
physical minimum-volume clamping and repeat binding. Full-writer checks verify
native quad/triangle records and zero-z output. Actual mechanical and visual
integration must be checked with the target GEOS/VisIt build.
