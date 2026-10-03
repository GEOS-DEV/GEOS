.. _uniform-refinement-scaling:

Uniform refinement scaling
==========================

``benchmarkUniformRefinement`` measures the current volume refinement components:
templates and geometry checks, local incidence, distributed IDs/sharing, and typed
point/cell field transfer. It uses the same implementation as the component tests.
The coupled ``AllMeshes`` controller and positive XML import have separate serial
and MPI correctness tests. This harness does not time that controller, fracture
associations, final GEOS connectivity, ghosts, DOFs, or restart. These measurements
cannot establish the scalability of those parts of the feature.

Build and run
------------

Use a VTK/MPI-enabled Release or RelWithDebInfo build, with the build directory
outside the source tree. Configure with ``ENABLE_BENCHMARKS=ON`` to register the
mesh benchmark targets::

  cmake --build /tmp/geos-uniform-refinement \
    --target benchmarkUniformRefinement --parallel 2
  mpiexec -n 4 /tmp/geos-uniform-refinement/bin/benchmarkUniformRefinement \
    --cells 16 16 16 --levels 2 --point-components 8 > hex.csv

``--cells`` specifies the global coarse Cartesian grid. ``--weak`` makes it the
grid per rank instead, multiplying its dimensions by the process grid. The
default process grid comes from ``MPI_Dims_create``; ``--process-grid PX PY PZ``
can compare compact and slab partitions. Children remain with their parent rank.
``--kind hex|tet|pyramid`` generates conforming grids; each coarse cube becomes
one hex, six tetrahedra, or six pyramids. Thus equal cube counts do not imply
equal refinement work across cell types.

A VTU can exercise a supplied mesh, including supported mixed cells/prisms::

  mpiexec -n 4 /tmp/geos-uniform-refinement/bin/benchmarkUniformRefinement \
    --input "$HOME/projects/geos-dataset/ALM_FIM_MGR/slantedFault/meshes/fault_contact_hex_mesh_ids.vtu" \
    --levels 1 > dataset.csv

Input reading on rank zero and Cartesian coarse scatter are reported separately.
Existing integral IDs remain exact; absent active global-ID arrays are seeded on
the coarse input. The harness excludes VTK extraction provenance arrays
``vtkOriginalPointIds``/``vtkOriginalCellIds`` and does not rebuild lineage. Those
are identifiers requiring production special handlers, not numerical fields.
This harness does not rebuild relational fracture data or apply
production node-set policy configuration. An unsupported field/geometry is an
error. ``--point-components`` and ``--weak`` apply to generated meshes only.

The default ``--discovery boundary --sharing interfaces`` extracts locally
exposed volume faces and their edges/vertices, plus separate ID-only records for
duplicate volume ownership checks. The latter still communicate O(coarse cells)
IDs and are included in ``coarseDiscoveryBoundaryAndVolumeIds``. Fine levels retain
only shared vertices/edges/faces, subdividing their parent traces without
rebuilding all fine volume incidence. Interior point participants are implicit.

``--discovery all --sharing all`` retains the diagnostic reference path that
discovers all coarse entities and rebuilds all local current-level incidence.
Compare both modes using matching partitions and input. Tests check exposed-face
candidate coverage against all-entity discovery on connected partition blocks,
including edge/vertex neighbors. Fine interface metadata matches the full
reference through the mixed triangle/quad and polygon-cap test matrix. The
production controller additionally runs a full coarse-face validation pass and,
when attachments exist, coarse actual-side association discovery. Normal coarse
vertex discovery uses boundary candidates. If any rank retains unused original
points, the compatibility path includes all original main-point IDs so an unused
copy can match an interior vertex on another rank. That exceptional path adds
O(coarse points) metadata; fine levels retain shared traces only. These production communication costs are
absent from this component harness; its boundary
candidate coverage assumes conforming manifold volume topology.

If any coarse triangle exists, the full-face pass also routes four triangle-subset
probes per unique coarse quad. Pure quadrilateral meshes skip these probes.
A triangle with three of that quad's global corners is diagnosed as a
nonmatching volume interface. These probes use the same coarse directory pass
and add no fine-level directory traffic. They are a topology check rather than
a spatial intersection search; conforming input remains a prerequisite.

Runs on Dane or another Slurm system
------------------------------------

Use an existing allocation and the site's GEOS host-config/compiler/TPL stack.
The runner launches job steps; it does not request an allocation or choose an
account/partition. Set ranks per node to a layout appropriate for that allocation.
For example, the following uses a chosen layout of 16 ranks per node::

  export GEOS_DANE_BUILD=/path/to/shared/gcc13-build
  export GEOS_DANE_RESULTS=/path/to/shared/refinement-results
  cmake --build "$GEOS_DANE_BUILD" \
    --target benchmarkUniformRefinement testMeshAdjacency testVTKRefinementTemplates \
             testVTKRefinementFields testVTKRefinementCommunication \
             testVTKUniformRefinement testVTKImport testVTKImport_mpi \
             testVTKRefinedFractureImport testVTKRefinedFractureImport_mpi \
             testVTKRefinedNodalInitialization testVTKRefinedNodalInitialization_mpi \
             testVTKRefinedFieldImport testVTKRefinedFieldImport_mpi \
             testVTKMeshScattering_mpi --parallel 2
  ctest --test-dir "$GEOS_DANE_BUILD" \
    -R '^(testMeshAdjacency|testVTKRefinement(Templates|Fields|Communication(_[248]ranks)?)|testVTKUniformRefinement(_[248]ranks)?|testVTKRefined(FractureImport|NodalInitialization|FieldImport)(_mpi|_[48]ranks)?|testVTKImport(_mpi|_pvtu_4ranks)?|testVTKMeshScattering_mpi)$' \
    --output-on-failure
  srun --nodes=1 --ntasks=4 --cpu-bind=cores \
    "$GEOS_DANE_BUILD/tests/testVTKImport_mpi" -x 4 \
    --gtest_filter=VTKImport.parallelFileWithFewerPiecesThanRanks

  export OMP_NUM_THREADS=1
  export VTK_SMP_MAX_THREADS=1
  python3 benchmarks/runUniformRefinementScaling.py \
    "$GEOS_DANE_BUILD/bin/benchmarkUniformRefinement" \
    "$GEOS_DANE_RESULTS/components" \
    --launcher "srun --cpu-bind=cores --cpus-per-task=1" --ranks-per-node 16 --nodes 1 2 4 \
    --strong-cells 32 32 32 --weak-cells 4 4 4 \
    --kinds hex tet pyramid --levels 1 2 --repeats 3 --point-components 8

The checked-in ``host-configs/LLNL/dane-toss_4_x86_64_ib-gcc@13.3.1.cmake``
selects GCC 13 and VTK 9.7. Local full builds and runtime checks pass with
Clang 23/VTK 9.4, AMD clang 22/ROCm 7.2/VTK 9.7 and a GCC 13.3.0/VTK 9.4
candidate containing separate serial-build prerequisites. Exact recipes and
source manifests accompany the local results. Run the tests with the actual
Dane stack before interpreting its performance results; local checks do not
establish multinode scalability.

The coupled tests include three levels of both supported-cell encodings,
source-separated pyramid descendants, marker and field transfer, remote fracture
sides, two named auxiliary namespaces, replicas, empty local blocks, and collective
failure. Import tests exercise ``uniformRefinement`` through the normal XML and
CellBlock path, including final fracture relations. The expanded fracture suite
checks triangular and all seven polygon-cap traces, separated contact sides,
two named junction blocks, and a connected fracture across rank boundaries.
It checks full region/ghost initialization, face/cell/node-side ordering, normals,
areas and final ID ranges. Isolated ranks exercise graph partition setup and
coloring. These cases pass locally through level two on 1/2/4/8/16 ranks.
The nodal suite calls the real FEM system setup and verifies unique global DOFs,
face/edge numbering, maps/maxima and normal ghost synchronization for edge-only
and vertex-only neighbors. An unrelated coincident-mesh control retains separate
DOFs. The field-import suite checks regular scalar/vector/tensor values and real
material stress/modulus values on connected mixed parents with separate original
block selectors but the same attribute. Both Float32/Float64 sources and
native/polyhedral encodings pass through level two on 1/2/4/8/16 local ranks.
The nodal fixtures also check overlapping named boundary node sets, isolated
corners and ghost memberships. Unselected varying integer/string point arrays
are included to verify that irrelevant metadata does not abort refinement.
The application harness adds
FVM/FEM, ghost/DOF, restart and sparse-ID well
checks. Run these checks on the actual compiler/TPL stack as well.

The persistent application correctness harness requires ``h5py`` and ``numpy``
only for checkpoint inspection::

  python3 benchmarks/verifyUniformRefinementInitialization.py \
    "$GEOS_DANE_BUILD/bin/geosx" "$GEOS_DANE_RESULTS/correctness" \
    --launcher "srun --cpu-bind=cores --cpus-per-task=1" --ranks 1 4 --wells

It runs FVM and FEM at levels zero through two with sparse IDs above ``2^53``.
Checks include ghost ownership, geometric coarse-parent recovery, partition-independent
solutions, an analytic uniaxial elastic patch with a ramped load, matching restart continuation
versus two uninterrupted steps, changed-level rejection, and a zero-level checkpoint with
the new optional wrapper removed to represent an older checkpoint.
The checkpoint for continuation precedes the second solve; the harness requires
solver iterations in the resumed log as well as agreement of final state.
``--wells`` adds two coupled wells at all three levels, including a fractured
reservoir with sparse surface IDs. At every level, well elements start above
the maximum existing element ID, including surface elements.
Checks cover exact element/node ranges and behavior near the limits of the
global-ID type. That optional well
regression uses the normal FGMRES/MGR reservoir/well solver configuration. It
also checks widening of the legacy global-next-element checkpoint array and
compares its values with the original 64-bit array.

The connected mixed-mesh correctness harness uses the committed hex/pentagonal-
prism/hexagonal-prism fixture. It verifies production pressure/porosity import
against coarse parents located from geometry, parent-volume recovery, constant-state
flow, closed-boundary mass conservation with nonzero internal fluxes, and an affine elasticity patch::

  python3 benchmarks/verifyUniformRefinementMixedSolvers.py \
    "$GEOS_DANE_BUILD/bin/geosx" "$GEOS_DANE_RESULTS/mixed-correctness" \
    --launcher "srun --cpu-bind=cores --cpus-per-task=1" --ranks 1 4

It checks levels zero through two, owned and ghost solutions and newly created
boundary nodes. The local 1/2/4/8/16-rank matrix passes all 45 runs. Add
``--verify-only`` to recheck the completed checkpoints without launching GEOS.
This verifies final nodal ownership and solutions; global DOF numbers are not
checkpointed. The separate live nodal test above verifies those numbers.

Linux shared builds provide an optional runtime call audit::

  cmake --build "$GEOS_DANE_BUILD" --target auditUniformRefinementCalls \
    testVTKUniformRefinement testVTKRefinedFractureImport_mpi --parallel 2
  ctest --test-dir "$GEOS_DANE_BUILD" -R '^testUniformRefinementAudit' \
    --output-on-failure

The preload library observes the positive refinement interval and forwards calls
unchanged. It counts root MPI Gather/Gatherv, GEOS scatter/redistribution,
VTK append/clean/merge and kd-tree filters, and enabled GEOS graph-partition
entry points. Tests reject nonzero counts and missing observations. An
injected-call self-test proves detection; named-junction import is audited at
four ranks. Use the audit for correctness, separately from timing runs.
It covers those dynamically linked entry points, not static/non-Linux builds
or arbitrary fine-data communication. Production statistics and source review
remain necessary to evaluate the communication model.

Start with small weak-scaling grids and increase after checking peak RSS. Pyramid
descendant counts grow differently from hex/tet counts. Run each invocation with
an exclusive allocation and fixed binding/thread settings. At least three repeats
help distinguish a regression from timing noise. Add ``--dry-run`` to inspect the
exact commands without launching MPI. Result directories must be new.

For input-mesh strong scaling, use ``--input FILE.vtu --modes strong``. For local
smoke runs use ``--launcher mpiexec --ranks-per-node 1 --nodes 1 2 4 8``; with
``mpiexec`` the runner sets process counts, while host placement is controlled by
the launcher options/environment. These local runs do not represent multiple nodes.

The output directory contains each run's CSV, stderr, command/status JSON, source
and executable hashes, git revision/status, build-cache settings, selected thread
and Slurm environment, and ``summary.csv``. Keep this directory with the matching
host-config, module list, and build log when comparing results. MPI implementation
version and actual node placement should also accompany shared results.

Full application measurements
-----------------------------

``runUniformRefinementApplicationScaling.py`` exercises the coupled controller,
normal GEOS connectivity, ghost construction, DOFs and one FVM solve. It retains
Caliper's timing report in ``caliper-report.txt`` and the existing GEOS
initialization/load-balance report. An explicit Caliper file avoids losing
the report to GEOS' redirected diagnostic stream during shutdown.
Every rank records its hostname, CPU affinity, elapsed process time and GNU
``time`` peak RSS. Checkpoint I/O is excluded. Its summary reports complete
application time, which includes coarse input/scatter and the solve; use the
Caliper scopes and the parsed hierarchy in ``phases.csv`` to assess refinement
separately. Inclusive scope times contain their descendants, so do not sum nested
rows. Level zero supplies the baseline for the existing initialization path.
``summary.csv`` separates refinement
efficiency from application efficiency and includes connectivity, ghost, DOF
and combined system setup timing, plus total linear iterations for the time step.
System setup includes DOF numbering, sparsity and vector allocation during the
first solve. Separate DOF/sparsity columns are populated only when those optional
scopes are present in the supplied build; otherwise they are empty.

The runner sets the mesh ``logLevel`` to two and saves growth forecasts and
per-level production statistics in ``runs.json``. It checks forecast cell counts
against the actual counts and requires zero fine-level directory exchanges.
``summary.csv`` includes final cell imbalance, unique/shared point counts,
maximum-rank payload messages/bytes, and resource estimates. The
``uniformRefinement/statistics`` scope measures reporting overhead within the
refinement timer. Level zero produces no refinement reports or allocations.

Byte estimates use conservative point-copy and retained-field payload bounds.
The refiner peak model includes the parent, child, registries, plans/descriptors,
retained coarse input and worst-case neighbor buffers. GEOS owned connectivity
is modeled separately; the ghost model assumes all other ranks' datasets could
be ghosted. Allocator/runtime overhead and physics/material fields are excluded.
These models are not a physical-memory budget. Compare them with measured
``peak_rank_rss_kib``, which includes process-lifetime allocator/runtime memory.
Models use bytes, while RSS uses KiB. The largest resource estimates and message
counts describe rank maxima; point/cell totals describe unique owners or copies
as indicated by each column.

For the four-node GCC 13 Dane allocation, start with this example layout::

  cmake --build "$GEOS_DANE_BUILD" --target geosx --parallel 8
  export OMP_NUM_THREADS=1
  export VTK_SMP_MAX_THREADS=1
  python3 benchmarks/runUniformRefinementApplicationScaling.py \
    "$GEOS_DANE_BUILD/bin/geosx" "$GEOS_DANE_RESULTS/application" \
    --launcher "srun --cpu-bind=cores --cpus-per-task=1" \
    --nodes 1 2 4 --ranks-per-node 16 \
    --strong-cells 32 32 32 --weak-cells 4 4 4 \
    --levels 0 1 2 --repeats 3

The generated mesh is a connected hex grid with a nonzero pressure-gradient
solve. Weak scaling grows the global coarse grid and physical domain in proportion to
the rank count, keeping coarse cell size and aspect ratio fixed.
An existing VTU, including supported mixed cells or prisms, can be measured with
``--mesh FILE.vtu --modes strong``; that input uses a uniform volume source and
closed boundaries. Named fracture blocks in a VTM can be measured with
``--main-block`` and ``--face-blocks``::

  python3 benchmarks/runUniformRefinementApplicationScaling.py \
    "$GEOS_DANE_BUILD/bin/geosx" "$GEOS_DANE_RESULTS/fracture-hex" \
    --mesh "$HOME/projects/geos-dataset/ALM_FIM_MGR/slantedFault/meshes/fault_contact_hex.vtm" \
    --main-block main --face-blocks fracture --modes strong \
    --launcher "srun --cpu-bind=cores --cpus-per-task=1" \
    --nodes 1 2 4 --ranks-per-node 16 --levels 0 1 --repeats 3

Use shared paths for the dataset as well as the binary/results. Repeat with
``fault_contact_tet.vtm`` in a new result directory. Start at levels zero/one;
the tet fixture has substantially more coarse cells. These runs retain exact
input global IDs for collocation references, initialize named fracture regions
and normal ghosts/DOFs, and solve matrix/fracture single-phase flow with fixed
aperture and parallel-plate permeability. They do not run contact mechanics.
``summary.csv`` separates fracture-association time, and ``phases.csv`` includes
``geos/MeshLevel/expandCollocatedGhostNodes`` within ghost preparation.
Provenance hashes both the VTM and all referenced mesh pieces.

Collocated ghost expansion visits fracture buckets once per neighbor/subregion
and adjacency depth, preserving the existing one-hop and subregion ordering.
It no longer rebuilds and scans the entire bucket set for every boundary node.
``testMeshAdjacency`` checks depth limits, overlapping/duplicate buckets, junctions,
missing remote nodes and named-subregion expansion, without VTK.

Once the initial weak run fits, repeat with ``--modes weak
--weak-cells 8 8 8`` and a new result directory to increase work per rank.
``--partition-refinement 1`` can exercise coarse graph partitioning;
``--partition-method ptscotch`` selects the ordinary PT-Scotch path.

Set both example paths to shared storage visible on every compute node;
node-local ``/tmp`` is unsuitable for these Slurm launches. The binary and its
TPL runtime libraries must also be accessible on every node. Keep the allocation
exclusive. The runner launches sequential steps inside the existing allocation;
it does not request nodes or change compiler modules. Its Python interpreter,
source files and ``/usr/bin/time`` must be available on those nodes. Add
``--dry-run`` to review generated inputs and commands. On a workstation,
``--launcher mpiexec --ranks-per-node 1 --nodes 1 2 4`` measures local process
counts rather than multiple physical nodes.

Reading the measurements
------------------------

CSV rows report minimum/mean/maximum rank time, cell-count imbalance, local point
replicas, unique owned points, shared point replicas, incidence counts, neighbor
degree, maximum rank RSS/high-water RSS, transmitted payload bytes/chunks, and
directory/neighbor exchange counts. Owned/shared counts are available after
incidence has been installed. Point allocation authority is not final GEOS DOF
ownership. RSS includes MPI/VTK/runtime memory and allocator retention; a peak is
a process-lifetime high-water measurement, not a fresh peak for each phase.

Each phase starts after a barrier. Report generation and those barriers are
excluded from its time. ``levelTotal`` sums the measured local phase durations,
then reports their rank maximum. It includes releasing the parent and temporary
state, and excludes the coarse read/scatter/discovery. Later levels must show zero
directory exchanges. Payload statistics omit self-routed records and count-control
bytes are reported separately; scalar reductions and MPI-internal traffic are not
included in the byte counters. Caliper scopes are available through the existing
GEOS profiling macros.

``summary.csv`` uses the median sum of ``levelTotal.seconds_max`` over repeats.
Strong efficiency is ``Tbaseline * Pbaseline / (T * P)``; weak efficiency is
``Tbaseline / T``. Baselines use the smallest process count in the same family.
Interpret throughput together with imbalance, shared replicas, RSS, and phase
times. Initial all-to-all count metadata has O(P) size per rank; later point
payloads use cached neighbors. ID-range count validation uses constant-size
reductions with checked limb sums before Exscan.

Acceptance remains local work/memory proportional to local fine entities and
communication proportional to actual shared entities. The final production gate
also requires integrated boundary-candidate/attachment validation, connected
mixed/prism and fracture measurements and production memory estimates.
Production estimates are available with mesh ``logLevel="2"``; actual multinode
memory and communication measurements remain necessary for acceptance.
Fine shared records use the layouts checked on coarse shared vertices/surface
replicas; coarse records still carry their full schema. The application report
now separates GEOS connectivity/ghost/DOF construction from refinement. Record
measured limitations; do not infer multinode scalability from small MPI tests.
