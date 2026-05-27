.. _ExternalMaterialModels:

External Material Models
========================

GEOS can optionally compile material models that live outside the public source tree.
This is intended for site-local or proprietary models that should not be committed to the repository.

Configure GEOS with the CMake cache variable ``GEOS_EXTERNAL_CONSTITUTIVE_MODELS_DIR`` set to the directory that contains the external model manifest.
The ``setupMPM`` helper also accepts ``--external-material-models DIR`` and forwards the directory to CMake for either Dane or Tuolumne builds.

The external directory should contain either ``GEOSExternalConstitutiveModels.cmake`` or ``CMakeLists.txt``.
That file can register model source files, include directories, dependencies, and MPM dispatch types with ``geos_register_external_constitutive_models``.
For example:

.. code-block:: cmake

  geos_register_external_constitutive_models(
    SOURCES
      MyInternalModel.cpp
    HEADERS
      MyInternalModel.hpp
    MPM_INCLUDES
      MyInternalModel.hpp
    MPM_MODELS
      MyInternalModel )

For explicit MPM builds, every external continuum model used by the solver must be listed in ``MPM_MODELS`` and every corresponding header must be listed in ``MPM_INCLUDES``.
External MPM model types are placed before GEOS built-in types in the pass-through dispatch list so derived site-local models are checked before their public base classes.

External cohesive-zone models can be registered the same way with ``MPM_COHESIVE_ZONE_MODELS`` and ``MPM_COHESIVE_ZONE_INCLUDES``.

The model implementation still needs to register itself with the normal GEOS constitutive catalog, for example by using ``REGISTER_CATALOG_ENTRY( ConstitutiveBase, MyInternalModel, std::string const &, Group * const )`` in the model source file.
