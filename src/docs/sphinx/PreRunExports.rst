.. _PreRunExports:

Machine-readable capabilities and global input metadata
=======================================================

``geosx --capabilities --format=json`` returns JSON before runtime/MPI setup.
The format/schema revision is 1; ``input-catalog`` is true and ``export-mesh``
is false in this focused PR. Consumers must check the feature booleans and
contract revision. The input catalog is global, never just the loaded instance.
Route vendor startup diagnostics to stderr (UCX_LOG_FILE=/dev/stderr when
applicable) so standalone capability stdout remains JSON.

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

Publication and regression
--------------------------

The destination must be new and its parent directory must exist. Work happens
in a hidden staging directory; an exclusive atomic link publishes the completed
catalog. Existing outputs are never replaced. Schema, restart, and validate-only
modes cannot be combined with catalog export. Failed exports return nonzero and
may leave unpublished staging files; only a complete final file is authoritative.
Caught errors include a JSON error with code ``pre-run-export-failed`` among
diagnostics; invalid standalone capability syntax uses
``invalid-capabilities-arguments`` and exits 2.

Run ``python3 scripts/testInputCatalog.py --geos /path/to/geosx``. This independent
regression checks global XSD parity, defaults/enum/bounds/unknown units,
repeatability with and without an input deck, errors and no-overwrite behavior.
