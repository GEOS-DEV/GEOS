.. _UsingSingularity:

[Unsupported] Using Singularity containers to deploy GEOS to users
====================================================================


An option to use containerized version of GEOS on cluster is to use `Singularity/Apptainer <https://apptainer.org/>`_ containers.
An interface allow user to pull and convert from docker repository.


CPU only container
-------------------

.. code-block::

    srun -c 4 singularity pull docker://jafranc/geosx_u22-g11-omp41

The docker image located at `docker://jafranc/geosx_u22-g11-omp41 <https://hub.docker.com/repository/docker/jafranc/geosx_u22-g11-omp41>`_
has been generated from `geosx/ubuntu22.04-gcc11:231-985 <https://hub.docker.com/r/geosx/ubuntu22.04-gcc11>`_
CI container with an extra layer cloning and compiling the develop version of GEOS. An example of a _Dockerfile_ producing such
a container can be found in sources at *src/docs/sphinx/developerGuide/Contributing/Docker_to_singularity.Dockerfile*. It is organized as:
Once this pull/convert operation done, a new `geosx_u22-g11-omp41_latest.sif` should have been created.


.. note::

   User can choose once generated, to set a central location for images in order to reuse.



.. code-block::

    /opt
        /GEOSX
            /GEOSX-version
            /GEOSX_TPL-231-985-9498edc/
            /inputFiles

respectively for GEOS installation (GEOSX-version), GEOSX TPL as generated in the root image  *geosx/ubuntu22.04-gcc11:231-985*
(/GEOSX_TPL-231-985-9498edc) and a copy of input files that contains all xml test files (inputFiles).

Then launching a shell inside the singularity image,

.. code-block::

    singularity shell geosx_u22-g11-omp41_latest.sif

and launching a 2-cores example.

.. code-block::

    Apptainer> mpirun -n 2 --oversubscribe /opt/GEOSX/GEOSX-version/bin/geosx -i /opt/GEOSX/inputFiles/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml -x 2

then, as you local directory is automatically mounted, if the results are produced locally, you should be able to inspect them once the virtualization ended.

.. note::
    It might be required to `export OMPI_MCA_opal_cuda_support=false` in order to disable cuda support loading of the ompi version as the base TPL image
might have CUDA compiled in it and might perturb the CPU only test at stake here.

.. note::

    During pull/convert operation Singularity will refuse to overwrite an existing image.
   Moreover, Singularity images are read-only container by default.

Now let us use that to launch `sbatch` computation.

.. code-block::

    #!/bin/sh
    #SBATCH --job-name=singularity_test-cpu
    #SBATCH --time=00:10:00
    #SBATCH -n 2

    module restore <your-module-list>

    export OMPI_MCA_opal_cuda_support=false
    export INPUT=/opt/GEOSX/inputFiles/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml

    #parallel executable
    singularity exec geosx_u22-g11-omp41_latest.sif /opt/GEOSX/GEOSX-version/bin/geosx -i ${INPUT} -x 2 -y 1 -z 1

GPU enable container
---------------------

Now if the target is a GPU-run then it will requires to use `jafranc/geos-gpu-test <https://hub.docker.com/repository/docker/jafranc/geos-gpu-test/general>`_
which is full 3-stage docker image based on `nvidia/cuda:11.5.2-devel-centos <https://hub.docker.com/r/nvidia/cuda>`_, building TPL then GEOS with proper settings.
The pull/convert command then is formulated as,

.. code-block::

    srun -c 12 singularity pull docker://jafranc/geos-gpu-test:11.5.2-sm80-devel

We will use either `11.5.2-sm80-devel` tag for targeting *sm_80* tag cards (e.g. A100-SXM4).

A batch file example to lauch them is then,

.. code-block::

    #!/bin/sh
    #BATCH --job-name=singularity_test-gpu
    #SBATCH --time=00:10:00
    #SBATCH -n 2
    #SBATCH -G 1
    #SBATCH -C GPU_GEN:AMP

    export IMAGE=geos-gpu-test_11.5.2-sm80-devel.sif

    module restore <your-module-list>

    ## either use in-image copied around xml input
    #export INPUT=/alt/geos/src/inputFiles/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
    ## or locally mounted
    export INPUT=compressible_1d.xml

    singularity exec --nv ${IMAGE} /usr/bin/nvidia-smi
    singularity exec --nv ${IMAGE} /alt/geos/src/build/bin/geosx -i ${INPUT}

.. note::

   The `--nv` option is required as it allows gpu-enable run in singularity