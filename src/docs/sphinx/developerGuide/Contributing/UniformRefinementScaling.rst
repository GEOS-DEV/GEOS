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
outside the source tree::

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
when attachments exist, coarse actual-side association discovery. Their separate
communication costs are absent from this component harness; its boundary
candidate coverage assumes conforming manifold volume topology.

Runs on Dane or another Slurm system
------------------------------------

Use an existing allocation and the site's GEOS host-config/compiler/TPL stack.
The runner launches job steps; it does not request an allocation or choose an
account/partition. Set ranks per node to a layout appropriate for that allocation.
For example, the following uses a chosen layout of 16 ranks per node::

  cmake --build /tmp/geos-uniform-refinement \
    --target benchmarkUniformRefinement testVTKRefinementTemplates \
             testVTKRefinementFields testVTKRefinementCommunication \
             testVTKUniformRefinement testVTKImport testVTKImport_mpi --parallel 2
  ctest --test-dir /tmp/geos-uniform-refinement \
    -R '^(testVTKRefinement(Templates|Fields|Communication(_[248]ranks)?)|testVTKUniformRefinement(_[248]ranks)?|testVTKImport(_mpi)?)$' \
    --output-on-failure

  export OMP_NUM_THREADS=1
  export VTK_SMP_MAX_THREADS=1
  python3 benchmarks/runUniformRefinementScaling.py \
    /tmp/geos-uniform-refinement/bin/benchmarkUniformRefinement \
    /tmp/uniform-refinement-scaling \
    --launcher "srun --cpu-bind=cores --cpus-per-task=1" --ranks-per-node 16 --nodes 1 2 4 \
    --strong-cells 32 32 32 --weak-cells 4 4 4 \
    --kinds hex tet pyramid --levels 1 2 --repeats 3 --point-components 8

The checked-in ``host-configs/LLNL/dane-toss_4_x86_64_ib-gcc@13.3.1.cmake``
selects GCC 13 and VTK 9.7. The local development build used Clang 23 and VTK 9.4;
run the component tests with the actual Dane stack before interpreting its
performance results. Local GCC 13 syntax checks do not replace that build/run.

The coupled tests include three levels of both supported-cell encodings,
source-separated pyramid descendants, marker and field transfer, remote fracture
sides, two named auxiliary namespaces, replicas, empty local blocks, and collective
failure. Import tests exercise ``uniformRefinement`` through the normal XML and
CellBlock path, including final fracture relations. Full FVM/FEM, ghost/DOF,
restart and sparse-ID well acceptance is still required.

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
mixed/prism and fracture measurements, compact field schemas, and
separate timings for subsequent GEOS connectivity/ghost/DOF construction. Record
measured limitations; do not infer multinode scalability from small MPI tests.
