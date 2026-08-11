# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 TotalEnergies
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved.
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

set( coverage_smoke_inputs
     materialPointMethod/singleParticle/mpm_singleParticle.xml
     materialPointMethod/particleRepartition/mpm_particleRepartition.xml
     materialPointMethod/dfgMovingGrid/mpm_dfgMovingGrid.xml
     lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_smoke.xml
     lagrangianContactMechanics/ContactMechanics_SimpleCubes_smoke.xml
     hydraulicFracturing/kgdViscosityDominated_smoke.xml
     hydraulicFracturing/kgdViscosityDominated_poroelastic_smoke.xml
     proppant/FlowProppantTransport2d_smoke.xml
     proppant/FlowProppantBedTransport2d_smoke.xml
     singlePhaseReactiveTransport/mixedReactionCarbonateSystemPorosityUpdateThermal_1D_smoke.xml
     phaseField/PhaseFieldFracture_Nucleation_smoke.xml
     phaseField/PhaseFieldPoromechanics_Nucleation_smoke.xml )

foreach( coverage_smoke_input IN LISTS coverage_smoke_inputs )
  set( input_path "${CMAKE_SOURCE_DIR}/../inputFiles/${coverage_smoke_input}" )
  if( NOT EXISTS "${input_path}" )
    message( FATAL_ERROR "Coverage smoke input does not exist: ${input_path}" )
  endif()

  get_filename_component( test_name "${coverage_smoke_input}" NAME_WE )
  set( output_directory "${CMAKE_CURRENT_BINARY_DIR}/${test_name}" )
  file( MAKE_DIRECTORY "${output_directory}" )

  add_test( NAME coverageSmoke_${test_name}
            COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/geosx"
                    -i "${input_path}"
                    -o "${output_directory}" )
  set_tests_properties( coverageSmoke_${test_name}
                        PROPERTIES
                        ENVIRONMENT "OMP_NUM_THREADS=1"
                        LABELS coverage_smoke
                        TIMEOUT 180
                        WORKING_DIRECTORY "${output_directory}" )
endforeach()
