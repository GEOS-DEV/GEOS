.. _UsingSingularity:

[Unsupported] Using Singularity containers to deploy GEOSX to users
====================================================================


An option to use containerized version of GEOSX on cluster is to use *singularity* containers.
An interface allow user to pull and convert from docker repository.

.. code-block::
   :caption: An example of singularity

    srun -c 4 --pty bash
    singularity pull docker://jafranc/geosx_u22-g11-omp41

The docker image located at *docker://jafranc/geosx_u22-g11-omp41* has been generated from *geosx/ubuntu22.04-gcc11:213-913*
with an extra layer cloning and compiling the develop version of GEOSX. It is organized as:

.. code-block::
    :caption: Organization of the docker/singularity image

    /opt
        /GEOSX
            /GEOSX-version
            /GEOSX_TPL-213-913-4dd5a33
            /inputFiles

respectively for GEOSX installation (GEOSX-version), GEOSX TPL as generated in the root image  *geosx/ubuntu22.04-gcc11:213-913*
(/GEOSX_TPL-213-913-4dd5a33) and a copy of input files that contains all xml test files (inputFiles).

Then launching a shell inside the singularity image,

.. code-block::
    :caption:

    singularity shell geosx_u22-g11-omp41_latest.sif

and launching a 2-cores example.

.. code-block::
    :caption:

    Apptainer> mpirun -n 2 --oversubscribe /opt/GEOSX/GEOSX-version/bin/geosx -i /opt/GEOSX/inputFiles/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml -x 2

then secure copying results back locally in order to process them.

.. note::
    It might be required to `export OMPI_MCA_opal_cuda_support=false` in order to disable cuda support loading of the ompi version.