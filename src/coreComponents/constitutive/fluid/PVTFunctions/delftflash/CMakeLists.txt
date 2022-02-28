cmake_minimum_required(VERSION 3.10)

#project(Flash)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++17")

#find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set(SOURCES
    main.cpp
    flash/flash.cpp flash/rr.cpp flash/ssi.cpp
    phase-correlations/eos_aq1.cpp phase-correlations/eos_aq2.cpp phase-correlations/eos_hydrate.cpp phase-correlations/eos_pr.cpp phase-correlations/eos_srk.cpp
    phase-correlations/k_initial.cpp phase-correlations/phase.cpp
    global/global.cpp global/misc.cpp
    )


target_include_directories(PVTPackage PRIVATE
    ${PROJECT_SOURCE_DIR}/flash
    ${PROJECT_SOURCE_DIR}/phase-correlations
    ${PROJECT_SOURCE_DIR}/global
"/data/gpfs/ixdev/rustem/GEOSX_MASTER/GEOSX/src/coreComponents/constitutive/PVTPackage/PVTPackage/eigen-3.3.9"
    )

