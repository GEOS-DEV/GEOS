# Specify solver headers
set( physicsPackages_headers
     ${physicsPackages_headers}
     surfaceGeneration/EmbeddedSurfaceGenerator.hpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.hpp
     surfaceGeneration/ParallelTopologyChange.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     surfaceGeneration/SurfaceGeneratorFields.hpp
     surfaceGeneration/kernels/surfaceGenerationKernels.hpp
     surfaceGeneration/kernels/surfaceGenerationKernelsHelpers.hpp
     PARENT_SCOPE )

# Specify solver sources
set( physicsPackages_sources
     ${physicsPackages_sources}
     surfaceGeneration/EmbeddedSurfaceGenerator.cpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.cpp
     surfaceGeneration/ParallelTopologyChange.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     PARENT_SCOPE )