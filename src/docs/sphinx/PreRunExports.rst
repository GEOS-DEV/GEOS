.. _PreRunExports:

Machine-readable pre-run exports
================================

These opt-in commands do not change ordinary simulation, restart, schema, or
``--validate-input`` behavior. They use the selected executable's registered
objects and native mesh writer; they do not run the event/time loop.

Capability discovery
--------------------

.. code-block:: sh

   geosx --capabilities --format=json

The standalone query runs before GEOS runtime/MPI initialization. Its JSON
object has ``formatVersion: 1``, ``version``, ``schemaVersion: 1``, boolean
``input-catalog`` and ``export-mesh`` keys, ``inputCatalogFormatVersion: 1``,
``inputCatalogScope: "global"``, ``meshExportFormatVersion: 1`` and
``meshExportFormat: "vtm"``. Builds without VTK advertise ``export-mesh: false``
and reject that command. Unknown members may be ignored. Consumers must
reject unsupported format versions instead of assuming compatibility.

``schemaVersion`` identifies this schema/catalog contract revision, not a
previously existing GEOS release-wide schema version. ``version`` is the
normal GEOS version/build string. The catalog adds a fingerprint of its exact
expanded XSD. Consumers should also identify the actual executable, since a
local patch may not change the source commit printed by ``getVersion()``.

The query explicitly advertises arbitrary-plane phase initialization,
generated-set phase initialization, and well-neighborhood refinement as false.
These physics/mesh extensions are not implemented by the export commands.

An invalid capability query exits 2 and writes a JSON ``error`` object to
stdout. Some vendor libraries can log *before main()*; route those diagnostics
to stderr in the launch environment. In particular, UCX supports
``UCX_LOG_FILE=/dev/stderr``. This preserves its diagnostics without mixing
them into the JSON response. The capability query is a serial query; do not
launch one copy per MPI rank.

Global input catalog
--------------------

.. code-block:: sh

   geosx --input-catalog new-catalog.json
   geosx -i input.xml --input-catalog new-catalog.json

Without ``-i``, the command expands all registered catalogs using the same
skeleton and schema-deviation code as ``--schema``. With ``-i``, it first
performs normal input setup and initial-condition validation, discards that
instance, and then builds the **same global catalog**. Thus a deck containing
one fluid model does not hide the other models available in that executable.
Both serial and MPI invocations are supported. No events are executed.
``mpiSize`` records the communicator size, since a few authoritative defaults
(such as VTK target processes) depend on it. Reproducibility comparisons must
use the same rank count; the schema fingerprint changes when those defaults do.

The top-level object contains:

* ``formatVersion`` and ``schemaVersion`` (both 1), ``geosVersion`` and
  ``scope: "global"``;
* ``schemaIdentity`` with ``algorithm: "fnv1a64"`` and a hexadecimal ``digest``
  of the generated XSD bytes. This deterministic identity is not a
  cryptographic authenticity claim;
* ``unitMetadata: "unavailable-in-wrapper-registry"`` for this GEOS pin;
* ``elements``, a deterministic flat array of the input schema's paths.

Each element has ``type`` (XML element name), ``schemaType``, a stable ``path``
such as ``/Problem/Solvers/CompositionalMultiphaseFVM``, ``group`` (parent path),
``description``, ``properties``, and ``children``. Paths are schema/type paths,
not names of live XML instances. A recursive schema edge is represented once
with ``recursive: true`` and its children, rather than expanded infinitely.
Child records contain ``name``, ``type``, ``minOccurs`` and ``maxOccurs``,
with XML Schema defaults (1) when an attribute is omitted. ``childChoice``
records the enclosing choice group's occurrence limits; its optional/repeated
semantics apply in addition to the child limits.

Each property includes:

* ``name`` and stable ``path`` ending in ``/@attributeName``;
* authoritative raw runtime ``type`` and XML ``schemaType``;
* ``description``, ``required``, ``inputFlag``, and optional string ``default``;
* ``pattern`` and ``choices``. Only an entire literal alternation is promoted
  to choices. Other validation regexes are preserved without interpretation;
* verbatim ``limits``, ``limitsMode`` (0 indicative, 1 warning, 2 error), and
  finite numeric ``minimum``/``maximum`` with boolean
  ``exclusiveMinimum``/``exclusiveMaximum`` where declared. For arrays these
  bounds apply to each value, matching the registered wrapper. Indicative or
  warning limits must not be presented as hard errors;
* ``units``, ``unitsStatus``, and ``mutability``.

This pin has no structured unit or mutability declaration on input wrappers.
Accordingly it emits ``units: null``, ``unitsStatus: "unknown"`` and
``mutability: "unknown"``. Unknown is **not** dimensionless. A future declared
dimensionless quantity uses ``units: "1"`` and ``unitsStatus: "declared"``.
No units or numeric constraints are guessed from property names, prose,
defaults, or a loaded simulation. Empty descriptions remain empty.

Authoritative pre-run mesh
--------------------------

.. code-block:: sh

   geosx -i input.xml --export-mesh new-preview.vtm
   mpirun -np 2 geosx -i input.xml -x 2 --export-mesh new-preview.vtm

GEOS parses and constructs the normal problem, including its mesh levels,
numerical methods and imported fields. It stops **before applying initial
conditions or running events**. Thus it validates mesh construction, not all
initial-condition physics. This is a geometry export, not a simulation result.

The native GEOS VTK writer supplies the geometry, VTK connectivity convention,
point ordering, region/subregion cell ordering, named body/level/region blocks
and per-rank datasets. Shallow-copy levels are omitted as in normal VTK output.
No rank aggregation, point merging, ghost filtering or decimation is requested.
Ghost cells are retained, with ``ghostRank`` and ``localToGlobalMap`` identity
arrays. Normal output using ``writeGhostCells="1"`` is the appropriate ordering
reference for the same partitioning. Particles/wells/surfaces retain their
normal native VTK representations. The export preserves real64 coordinates
(``coordinatePrecision: "float64"``); normal VTK output retains its previous
precision. This prevents small cells at large coordinates from collapsing
through a float32 conversion.

The result is ``new-preview.vtm`` and ``new-preview.vtm.data/``. The latter
contains VTU sidecars and ``metadata.json`` with ``formatVersion: 1``,
``kind: "geos-pre-run-mesh"``, ``geosVersion``,
``ordering: "native-vtk-writer"``, ``ghostCells: true``,
``timeLoopEntered: false`` and ``initialConditionsApplied: false``.
Coordinate units are explicitly ``null``/``unknown`` because this pin does not
carry a structured coordinate-unit declaration. Consumers must not silently
apply a metres conversion. No time metadata, PVD file, result fields or
simulation timestep is exported. Internal directory names inherited from the
VTK writer do not represent executed timesteps.

Publication and errors
----------------------

The destination must be a **new file in an existing directory**. A mesh must
have the ``.vtm`` suffix and its ``.data`` sibling must not already exist.
Existing outputs are never replaced. Export modes cannot be combined with
each other, schema generation, restart, or ``--validate-input``.

Work is confined to a hidden staging directory in the destination filesystem.
The catalog is published with an exclusive atomic link after completion. For
a mesh, all complete sidecars are moved into place first; an exclusive atomic
link publishes the root VTM last as the **commit marker**. Until that marker
exists, the bundle is not authoritative. Interrupted/failed processes can leave
hidden staging directories or unreferenced sidecars; consumers must never treat
those as completed exports. Successful export removes its staging directory.
Concurrent replacement/overwrite of export destinations is unsupported.

Errors return nonzero and do not publish a success marker. Caught export errors
also emit a one-line JSON object with ``error.code: "pre-run-export-failed"``
and ``error.message`` on stdout alongside normal GEOS diagnostics. Fatal
runtime/OS termination can preclude that diagnostic, so exit status and the
absence of the output remain authoritative. Capability errors use
``invalid-capabilities-arguments``. Never accept partial JSON or a VTM with
missing sidecars as a valid export.

Regression checks
-----------------

.. code-block:: sh

   python3 scripts/testPreRunExports.py --geos /path/to/geosx --mpiexec mpirun

The test requires Python VTK. It constructs a small real mesh and compares
geometry, block identities, point/cell ordering and identity arrays against
normal same-build solver output in serial and with two MPI ranks. It also
checks catalog scope, repeatability with/without a deck, unknown-unit handling,
choices/constraints, rejected modes, invalid input, preservation of existing
outputs, complete sidecars and the absence of simulation metadata. Omitting
``--mpiexec`` explicitly skips the MPI test rather than claiming it passed.
