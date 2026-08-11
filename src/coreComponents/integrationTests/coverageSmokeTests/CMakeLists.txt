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
     phaseField/PhaseFieldPoromechanics_Nucleation_smoke.xml
     inducedSeismicity/SeismicityRate_analytical_verification_smoke.xml
     inducedSeismicity/SpringSliderExplicit_A_smoke.xml
     inducedSeismicity/SpringSliderImplicit_A_smoke.xml
     efemFractureMechanics/EmbFrac_Compression_CoulombFriction_smoke.xml
     poromechanics/PoroDelftEggWellbore_smoke.xml
     poromechanics/PoroViscoDruckerPrager_smoke.xml
     poromechanicsFractures/Contact/LagrangianMultipliers/PoroElastic_conformingFracture_2d_faultSlip_smoke.xml
     poromechanicsFractures/EDFM/PoroElastic_efem-edfm_verticalFrac_smoke.xml
     thermoPoromechanics/ThermoPoroDruckerPrager_consolidation_smoke.xml
     thermoPoromechanicsFractures/ThermoPoroElastic_conforming_smoke.xml
     simplePDE/10x10x10Hex_LaplaceVEM_smoke.xml
     solidMechanics/elasticHollowCylinder_orthotropic_smoke.xml
     wellbore/CasedThermoElasticWellbore_ImperfectInterfaces_smoke.xml
     compositionalMultiphaseFlow/2ph_cap_1d_ihu.xml
     compositionalMultiphaseFlow/soreideWhitson/lockExchange/lockExchange_Z.xml
     poromechanics/PoroElastic_staircase_co2_3d_fim_stab.xml
     compositionalMultiphaseWell/staircase_co2_wells_3d.xml
     constitutiveDriver/testPVT.xml
     constitutiveDriver/testRelperm_Table_2_phase.xml
     constitutiveDriver/testTriaxial_druckerPragerExtended.xml
     thermalMultiphaseFlow/co2_thermal_2d.xml
     compositionalMultiphaseFlow/4comp_2ph_1d.xml
     wavePropagation/acous3D_pml_smoke.xml
     multipleMeshBodies/testMultipleBodies.xml
     singlePhaseWell/compressible_single_phase_wells_1d.xml
     simplePDE/10x10x10Hex_LaplaceFEM_smoke.xml
     poromechanics/PoroElastic_staircase_singlephase_3d_fim_stab.xml
     phaseField/PhaseFieldFracture_CohesiveModel.xml
     singlePhaseFlow/incompressible_pebi3d.xml
     meshGeneration/multiBodyMeshGen_smoke.xml
     phaseField/PhaseFieldFracture_DamageOnly.xml )

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
