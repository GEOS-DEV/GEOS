.. _RunningOnSherlock:

[Unsupported] Running GEOS on Sherlock Stanford Cluster
==========================================================

`Sherlock <https://www.sherlock.stanford.edu/docs/>`_ Stanford cluster is a new resource that offer mixed CPU/GPU architecture.
We will then differentiate two profile in the following documentation on how to interact with it

  - *Users* who are only interested in running the head version of the code to benefit from
    its features. They will focus effort on building cases and preprocessed inputs and/or run
    multiple realisations with parameters sensitivities,

  - *Developers* who are interested in testing new feature that they developed in a cluster context
    and/or do limited debugging and tests of their code versions.

Though *Users* can follow developers path on sync to `origin/develop` and get to the same results, there is
quicker way for them to have such a version on Sherlock.

All data and resources on Sherlock will be centrally located at: `/oak/stanford/schools/ees/COLLABORATIONS/geosx`. Note
this is a per-request access but is agnostic of groups.

Prior to code location discussion, we stress the fact that beyond explicit dependencies that are either compiled in the
third party libraries (TPL), e.g. HYPRE, or as submodules, e.g. LvArray, some dependencies are assumed to be system, such
as the compiler choices (gcc, clang, intel), the cuda and mpi implementations and others as git or cmake.

The appropriate pathing for those is set through the module manager thanks to successive module loads as described in
`sherlock-modules.sh <https://github.com/GEOS-DEV/GEOS/blob/develop/host-configs/Stanford/sherlock-modules.sh>`_.
An equivalent result can be obtain in restoring a module state previously saved (cf. `module save <https://lmod.readthedocs.io/en/latest/010_user.html#user-collections>`_).
All these states are login-wise saved at `$HOME/.lmod.d/` if one exist. If not, it is created at the first saved list.
To this purpose, several module list for GEOS are located in `/oak/stanford/schools/ees/COLLABORATIONS/geosx/lua`.

Using auto-deployed version
----------------------------

Under `/oak/stanford/schools/ees/COLLABORATIONS/geosx/CPU`, both *Users* and *Developers* will find binaries sorted
respectively of GEOS under `GEOS-<current_develop_commit_hash>/` and third party libs (TPL) under `GEOSX_TPL-<TPL_tag>-<current_commit_hase>/`.
These are regularly fetched and deployed at this location.

.. note::

   There exists a `/oak/stanford/schools/ees/COLLABORATIONS/geosx/GPU` that will serve the same purpose but for CUDA aware binaries.

These binaries are generated and uploaded into a Google Could bucket as artifacts of successful CI/CD stages.
Precompiled binaries for third party libraries can be useful for both *Users* and *Developers* as they offers
easy to reached certified versions of the external dependencies. Hence the *Developer* can safely link against,
provided that their feature branch is in sync with the target version of the TPL.

The *Users* can directly use the binary at `/oak/stanford/schools/ees/COLLABORATIONS/geosx/CPU/GEOS-<current_develop_commit_hash>/bin/geos`
which is a deployed version precompiled for Sherlock.

Using Singularity containers
-----------------------------

Another option for *Users* that want to solely use GEOS is to use containerized environment. The one available on Sherlock
is *Singularity/Apptainer*. We refer interested users to :ref:`UsingSingularity`. The generated images for CPU-only and GPU
environement are located at `/oak/stanford/schools/ees/COLLABORATIONS/geosx/sif`. Examples of batch run files for those are
stored in `/oak/stanford/schools/ees/COLLABORATIONS/geosx/batch-examples`.

Using self-compiled version
----------------------------

Whether they chose to use the precompiled TPL or compile them themselves, *Developers* will have to go through the cloning
steps explained at :ref:`QuickStart`. The appropriate module to load can be imported inspired from the shell script at
`sherlock-modules.sh <https://github.com/GEOS-DEV/GEOS/blob/develop/host-configs/Stanford/sherlock-modules.sh>`_ or the module
savelist lua file at `/oak/stanford/schools/ees/COLLABORATIONS/geosx/lua`.

