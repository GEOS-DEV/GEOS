"""PFW layout selection and compact diagnostics; no new GEOS particle state.

The hexTet preset selects existing native simplex cutting and interpolation.
Particle geometry/type statistics count only final emitted records, never
discarded coarse candidates or temporary clipping simplices.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class SurfaceLayout:
    method: str
    geometry: str
    simplex_particle_type: object
    hex_tet: bool


def resolve_surface_layout(method, geometry, reference_geometry, simplex_type,
                           plane_strain, reference_geometry_explicit=False):
    """Resolve hexTet without altering defaults for existing cartesian/simplex inputs.

    An explicit preferred geometry overrides the old reference-geometry alias,
    as before. A contradictory box request is rejected in the new preset.
    Point interpolation on cut tets remains available for comparisons; omitting
    the cut type selects CPTI. Object bulk types are checked separately after
    material/type metadata is collected from every writer rank.
    """
    method = str(method).lower().replace('-conforming', '')
    hex_tet = method == 'hextet'
    if hex_tet:
        if plane_strain:
            raise ValueError("surfaceRefinementMethod='hexTet' requires planeStrain=0; true 2D is separate")
        if geometry is not None and geometry != 'fromRefinement':
            raise ValueError("hexTet requires particleGeometry='fromRefinement'; remove the explicit box setting")
        if geometry is None and reference_geometry_explicit and reference_geometry != 'fromRefinement':
            raise ValueError("hexTet conflicts with particleReferenceGeometry='equivalentBox'; use fromRefinement")
        geometry = 'fromRefinement'
        simplex_type = 3 if simplex_type is None else simplex_type
        method = 'simplex'
    if geometry is None:
        geometry = 'fromRefinement' if reference_geometry == 'fromRefinement' else 'box'
    return SurfaceLayout(method, geometry, simplex_type, hex_tet)


def validate_hex_tet_bulk_type(particle_type):
    """Hex bulk requires an interpolation type that accepts box geometry."""
    if particle_type not in (0, 1, 2):
        raise ValueError("hexTet bulk objects must use particleType 0, 1, or 2 (CPDI); "
                         "particleType=3 would tetrahedralize the bulk")


class ParticleLayoutStatistics:
    """Counts and compensated physical-volume sums for emitted shape/type pairs.

    One small record per pair is gathered from each rank. This does not gather
    particle arrays, retain per-particle metadata, or parse the completed file.
    Aggregate volumes may differ by summation roundoff across rank counts;
    the generated particle records and integer counts remain deterministic.
    """
    def __init__(self):
        self.pairs = {}

    def record(self, particle_type, shape, volume):
        if not math.isfinite(volume) or volume <= 0:
            raise ValueError('Layout statistics require finite positive particle volume')
        key = int(particle_type), int(shape)
        # Count, compensated sum, compensation, smallest/largest particle.
        data = self.pairs.setdefault(key, [0, 0., 0., volume, volume])
        data[0] += 1
        value = volume - data[2]
        total = data[1] + value
        data[2] = (total - data[1]) - value
        data[1] = total
        data[3] = min(data[3], volume)
        data[4] = max(data[4], volume)

    def report(self):
        return [dict(particle_type=key[0], geometry_shape=key[1], count=data[0],
                     volume=data[1], min_particle_volume=data[3], max_particle_volume=data[4])
                for key, data in sorted(self.pairs.items())]


def combine_layout_statistics(reports):
    """Merge rank summaries, preserving empty-rank participation."""
    pairs = {}
    for report in reports:
        for row in report:
            key = row['particle_type'], row['geometry_shape']
            pairs.setdefault(key, []).append(row)
    rows = []
    shape_names = {0: 'hex', 1: 'wedge', 2: 'tet', 3: 'quad', 4: 'tri'}
    for key, values in sorted(pairs.items()):
        rows.append(dict(particle_type=key[0], geometry_shape=key[1],
                         shape_name=shape_names[key[1]], count=sum(v['count'] for v in values),
                         volume=math.fsum(v['volume'] for v in values),
                         min_particle_volume=min(v['min_particle_volume'] for v in values),
                         max_particle_volume=max(v['max_particle_volume'] for v in values)))
    return dict(total_count=sum(row['count'] for row in rows),
                total_volume=math.fsum(row['volume'] for row in rows),
                by_type_and_shape=rows)
