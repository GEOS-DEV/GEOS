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
     phaseField/PhaseFieldFracture_DamageOnly.xml
     compositionalMultiphaseFlow/deadoil_3ph_staircase_obl_3d.xml
     constitutiveDriver/testTriaxial_delftEggCase1.xml
     constitutiveDriver/testRelperm_VanGenuchten.xml
     singlePhaseReactiveTransport/mixedReactionCarbonateSystem_1D_smoke.xml
     wavePropagation/acouselas3D_Q2_abc_smoke.xml
     compositionalMultiphaseWell/resvol_constraint.xml
     compositionalMultiphaseWell/isothm_mass_inj_table.xml
     singlePhaseFlow/PressurePerm_1d_smoke.xml
     solidMechanics/elasticHollowCylinder_transverseIsotropic_smoke.xml
     constitutiveDriver/testRelperm_TableHyst2ph.xml
     thermalMultiphaseFlow/co2_thermal_obl_3d.xml
     compositionalMultiphaseFlow/c1-ppu/grav_seg_c1ppu_hyst.xml
     singlePhaseWell/perf_status_test.xml
     solidMechanics/sedov_finiteStrain_smoke.xml
     compositionalMultiphaseWell/benchmarks/Class09Pb3/class09_pb3_smoke_3d.xml
     compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_smoke_3d.xml
     wavePropagation/acous3D_vti_smoke.xml
     compositionalMultiphaseFlow/co2_flux_dirichlet.xml
     compositionalMultiphaseFlow/dbc/grav_seg_1d/grav_seg_1d_DBC.xml
     multiphaseFlowFractures/deadOil_fractureMatrixFlow_edfm_horizontalFrac_smoke.xml
     compositionalMultiphaseFlow/co2_flux_3d.xml
     solidMechanics/sedov_ssle_smoke.xml
     singlePhaseFlow/FieldCaseTutorial3_Thermal_smoke.xml
     constitutiveDriver/testTriaxial_elasticIsotropicPressureDependent.xml
     singlePhaseWell/staircase_single_phase_wells_3d.xml
     compositionalMultiphaseWell/staircase_co2_wells_hybrid_3d.xml
     thermoPoromechanics/ThermoPoroElastic_staircase_co2_smoke.xml
     thermoPoromechanicsFractures/ThermoPoroElastic_efem-edfm_verticalFrac_smoke.xml
     poromechanicsFractures/Contact/AugmentedLagrangianMultipliers/ALM_inclinedFault_ISG_smoke.xml
     efemFractureMechanics/Sneddon_embeddedFrac_staticCondensation_smoke.xml
     thermalSinglePhaseFlowFractures/fractureMatrixThermalFlow_edfm_smoke.xml
     phaseField/PhaseFieldFracture_VolDevSplit.xml
     inducedSeismicity/SeismicityRate_poromechanics_smoke.xml
     constitutiveDriver/testPVT_CO2BrineTables.xml
     singlePhaseFlow/sourceFlux_1d.xml
     singlePhaseFlowFractures/fractureFlow_conforming_2d.xml
     simplePDE/hybridHexPrism_LaplaceVEM_smoke.xml
     singlePhaseFlow/vtk/3D_10x10x10_compressible_hex_gravity_smoke.xml
     singlePhaseFlow/vtk/3D_10x10x10_compressible_tetra_gravity_smoke.xml
     singlePhaseFlowFractures/fractureMatrixFlow_edfm_horizontalFrac_smoke.xml
     solidMechanics/gravity.xml
     lagrangianContactMechanics/ContactMechanics_PEBICrack_smoke.xml
     singlePhaseFlow/FieldCaseTutorial3_Isothermal_smoke.xml
     lagrangianContactMechanics/ContactMechanics_UnstructuredCrack_smoke.xml
     lagrangianContactMechanics/ContactMechanics_Sneddon_smoke.xml
     lagrangianContactMechanics/ContactMechanics_TFrac_smoke.xml
     lagrangianContactMechanics/ContactMechanics_SingleFracCompression_smoke.xml
     lagrangianContactMechanics/ContactMechanics_PassingCrack_smoke.xml
     lagrangianContactMechanics/ALM_UnstructuredCrack_smoke.xml
     lagrangianContactMechanics/ALM_Sneddon_smoke.xml
     lagrangianContactMechanics/ALM_TFrac_smoke.xml
     lagrangianContactMechanics/ALM_SingleFracCompression_smoke.xml
     lagrangianContactMechanics/ALM_PassingCrack_smoke.xml )

function( add_coverage_smoke_test )
  set( options )
  set( one_value_args NAME MPI_RANKS SCHEMA TIMEOUT )
  set( multi_value_args INPUTS GEOS_ARGS )
  cmake_parse_arguments( arg
                         "${options}"
                         "${one_value_args}"
                         "${multi_value_args}"
                         ${ARGN} )

  if( arg_UNPARSED_ARGUMENTS )
    message( FATAL_ERROR "Unknown arguments for coverage smoke test: ${arg_UNPARSED_ARGUMENTS}" )
  endif()
  if( NOT arg_NAME )
    message( FATAL_ERROR "Coverage smoke test requires NAME" )
  endif()
  if( arg_INPUTS AND arg_SCHEMA )
    message( FATAL_ERROR "Coverage smoke test ${arg_NAME} cannot use both INPUTS and SCHEMA" )
  endif()
  if( NOT arg_INPUTS AND NOT arg_SCHEMA )
    message( FATAL_ERROR "Coverage smoke test ${arg_NAME} requires INPUTS or SCHEMA" )
  endif()

  set( ctest_name "coverageSmoke_${arg_NAME}" )
  if( TEST "${ctest_name}" )
    message( FATAL_ERROR "Duplicate coverage smoke test name: ${ctest_name}" )
  endif()

  set( input_arguments )
  foreach( coverage_smoke_input IN LISTS arg_INPUTS )
    set( input_path "${PROJECT_SOURCE_DIR}/../inputFiles/${coverage_smoke_input}" )
    if( NOT EXISTS "${input_path}" )
      message( FATAL_ERROR "Coverage smoke input does not exist: ${input_path}" )
    endif()
    list( APPEND input_arguments -i "${input_path}" )
  endforeach()

  set( output_directory "${CMAKE_CURRENT_BINARY_DIR}/${arg_NAME}" )
  file( MAKE_DIRECTORY "${output_directory}" )

  set( test_command "$<TARGET_FILE:geosx>" )
  if( arg_MPI_RANKS )
    if( NOT MPIEXEC_EXECUTABLE )
      message( FATAL_ERROR "Coverage smoke test ${arg_NAME} requires MPI" )
    endif()
    set( test_command
         "${MPIEXEC_EXECUTABLE}"
         "${MPIEXEC_NUMPROC_FLAG}"
         "${arg_MPI_RANKS}"
         ${MPIEXEC_PREFLAGS}
         "$<TARGET_FILE:geosx>"
         ${MPIEXEC_POSTFLAGS} )
  endif()

  if( arg_SCHEMA )
    set( run_arguments -s "${output_directory}/${arg_SCHEMA}" ${arg_GEOS_ARGS} )
  else()
    set( run_arguments ${input_arguments} ${arg_GEOS_ARGS} -o "${output_directory}" )
  endif()

  add_test( NAME "${ctest_name}"
            COMMAND ${test_command} ${run_arguments} )

  if( NOT arg_TIMEOUT )
    set( arg_TIMEOUT 180 )
  endif()
  set_tests_properties( "${ctest_name}"
                        PROPERTIES
                        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1"
                        LABELS coverage_smoke
                        TIMEOUT "${arg_TIMEOUT}"
                        WORKING_DIRECTORY "${output_directory}" )
  if( arg_MPI_RANKS )
    set_tests_properties( "${ctest_name}"
                          PROPERTIES
                          PROCESSORS "${arg_MPI_RANKS}" )
  elseif( GEOS_ENABLE_LLVM_SOURCE_COVERAGE )
    # Each instrumented process flushes every component DSO profile at exit
    # (about 70 MiB per serial smoke). Account for that I/O weight so a
    # --parallel 12 run does not issue twelve simultaneous profile flushes.
    set_tests_properties( "${ctest_name}"
                          PROPERTIES
                          PROCESSORS 4 )
  endif()
endfunction()

foreach( coverage_smoke_input IN LISTS coverage_smoke_inputs )
  get_filename_component( test_name "${coverage_smoke_input}" NAME_WE )
  set( coverage_smoke_overlays )

  if( test_name STREQUAL "ALM_inclinedFault_ISG_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/ALM_inclinedFault_small_mesh.xml )
  elseif( test_name STREQUAL "staircase_co2_wells_hybrid_3d" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/staircase_co2_wells_hybrid_small_mesh.xml )
  elseif( test_name STREQUAL "acouselas3D_Q2_abc_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/acouselas3D_Q2_small_mesh.xml )
  elseif( test_name STREQUAL "ALM_TFrac_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/TFrac_small_mesh.xml )
  elseif( test_name STREQUAL "ThermoPoroElastic_staircase_co2_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/ThermoPoroElastic_staircase_small_mesh.xml )
  elseif( test_name STREQUAL "PoroElastic_staircase_co2_3d_fim_stab" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/PoroElastic_staircase_fim_small_mesh.xml )
  elseif( test_name STREQUAL "ThermoPoroDruckerPrager_consolidation_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/ThermoPoroDruckerPrager_fast_events.xml )
  elseif( test_name STREQUAL "LagrangeContactBubbleStab_FixedSlip_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/LagrangeContactBubbleStab_small_mesh.xml )
  elseif( test_name STREQUAL "PhaseFieldFracture_Nucleation_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/PhaseFieldFracture_Nucleation_fast_events.xml )
  elseif( test_name STREQUAL "FieldCaseTutorial3_Thermal_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/FieldCaseTutorial3_Thermal_small_mesh.xml )
  elseif( test_name STREQUAL "FieldCaseTutorial3_Isothermal_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/FieldCaseTutorial3_Isothermal_small_mesh.xml )
  elseif( test_name STREQUAL "SeismicityRate_poromechanics_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/SeismicityRate_small_mesh.xml )
  elseif( test_name STREQUAL "co2_thermal_obl_3d" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/co2_thermal_obl_small_mesh.xml )
  elseif( test_name STREQUAL "deadOil_fractureMatrixFlow_edfm_horizontalFrac_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/deadOil_edfm_small_mesh.xml )
  elseif( test_name STREQUAL "SpringSliderImplicit_A_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/SpringSliderImplicit_fast_events.xml )
  elseif( test_name STREQUAL "staircase_co2_wells_3d" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/staircase_co2_wells_small_mesh.xml )
  elseif( test_name STREQUAL "deadOilEgg_smoke_3d" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/deadOilEgg_small_mesh.xml )
  elseif( test_name STREQUAL "ContactMechanics_SimpleCubes_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/ContactMechanics_SimpleCubes_small_mesh.xml )
  elseif( test_name STREQUAL "grav_seg_1d_DBC" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/grav_seg_1d_DBC_small_mesh.xml )
  elseif( test_name STREQUAL "co2_flux_dirichlet" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/co2_flux_dirichlet_small_mesh.xml )
  elseif( test_name STREQUAL "co2_flux_3d" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/co2_flux_3d_small_mesh.xml )
  elseif( test_name STREQUAL "kgdViscosityDominated_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/kgdViscosityDominated_small_mesh.xml )
  elseif( test_name STREQUAL "FlowProppantBedTransport2d_smoke" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/FlowProppantBedTransport2d_fast_outputs.xml )
  elseif( test_name STREQUAL "isothm_mass_inj_table" )
    list( APPEND coverage_smoke_overlays
          coverageSmokes/isothm_mass_inj_table_fast_outputs.xml )
  endif()

  add_coverage_smoke_test( NAME "${test_name}"
                           INPUTS "${coverage_smoke_input}" ${coverage_smoke_overlays} )

  if( GEOS_ENABLE_LLVM_SOURCE_COVERAGE AND
      ( test_name STREQUAL "ContactMechanics_SimpleCubes_smoke" OR
        test_name STREQUAL "LagrangeContactBubbleStab_FixedSlip_smoke" ) )
    # These nonlinear contact matrices remain sensitive to concurrent
    # instrumented solves and profile flushes.
    set_tests_properties( "coverageSmoke_${test_name}"
                          PROPERTIES
                          PROCESSORS 12 )
  endif()
endforeach()

# Higher-order DG tetrahedra have very high template multiplicity. Initialize order 3 separately
# so the order-2 timestep and order-3 setup form a bounded function-coverage milestone.
add_coverage_smoke_test(
  NAME acous3D_dg_order2
  INPUTS wavePropagation/acous3D_dg_order2_coverage.xml )

add_coverage_smoke_test(
  NAME acous3D_dg_order3_setup
  INPUTS wavePropagation/acous3D_dg_order2_coverage.xml
         wavePropagation/acous3D_dg_order3_overlay.xml
         wavePropagation/acous3D_dg_setup_only_overlay.xml )

# One order-3 timestep atomically exercises the generated mass, stiffness, and surface kernels.
add_coverage_smoke_test(
  NAME acous3D_dg_order3
  INPUTS wavePropagation/acous3D_dg_order2_coverage.xml
         wavePropagation/acous3D_dg_order3_overlay.xml )

# Schema generation expands every registered object catalog into a complete data-repository tree.
add_coverage_smoke_test(
  NAME schema_catalog
  SCHEMA catalog_schema.xsd )

# The second input is an ordered overlay that reduces each cube problem to a one-cycle smoke test.
foreach( multiscale_solver IN ITEMS singlePhaseFlow solidMechanics singlePhasePoromechanics )
  add_coverage_smoke_test(
    NAME "multiscale_${multiscale_solver}_msrsb"
    INPUTS "multiscalePreconditioner/${multiscale_solver}/cube_msrsb.xml"
           "multiscalePreconditioner/cube_coverage_smoke.xml" )
endforeach()

# These decompositions cover topology migration and partition-only kernels. PROCESSORS prevents
# CTest --parallel from oversubscribing the twelve CPU slots used for the coverage build.
add_coverage_smoke_test(
  NAME SurfaceGenerator_mpi_2x2x2
  MPI_RANKS 8
  INPUTS surfaceGeneration/SurfaceGenerator.xml
  GEOS_ARGS -x 2 -y 2 -z 2 )

add_coverage_smoke_test(
  NAME mpm_particleRepartition_mpi_2x2x1
  MPI_RANKS 4
  INPUTS materialPointMethod/particleRepartition/mpm_particleRepartition.xml
  GEOS_ARGS -x 2 -y 2 -z 1 )

add_coverage_smoke_test(
  NAME mpm_dfgMovingGrid_mpi_2x2x1
  MPI_RANKS 4
  INPUTS materialPointMethod/dfgMovingGrid/mpm_dfgMovingGrid.xml
  GEOS_ARGS -x 2 -y 2 -z 1 )
