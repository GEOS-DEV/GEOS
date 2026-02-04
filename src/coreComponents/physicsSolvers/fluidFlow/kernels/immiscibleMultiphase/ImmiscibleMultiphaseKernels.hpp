/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ImmiscibleMultiphaseKernels.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_MULTIPHASEKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_MULTIPHASEKERNELS_HPP


 #include "codingUtilities/Utilities.hpp"
 #include "common/DataLayouts.hpp"
 #include "common/DataTypes.hpp"
 #include "common/GEOS_RAJA_Interface.hpp"
 #include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluid.hpp"
 #include "constitutive/solid/CoupledSolidBase.hpp"
 #include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluidFields.hpp"
 #include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
 #include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
 #include "constitutive/capillaryPressure/JFunctionCapillaryPressure.hpp"
 #include "constitutive/capillaryPressure/TableCapillaryPressure.hpp"
 #include "constitutive/permeability/PermeabilityBase.hpp"
 #include "constitutive/permeability/PermeabilityFields.hpp"
 #include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
 #include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
 #include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"

#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluid.hpp"


 #include "fieldSpecification/AquiferBoundaryCondition.hpp"
 #include "finiteVolume/BoundaryStencil.hpp"
 #include "finiteVolume/CellElementStencilTPFA.hpp"
 #include "finiteVolume/FluxApproximationBase.hpp"
 #include "linearAlgebra/interfaces/InterfaceTypes.hpp"
 #include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
 #include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
 #include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlowFields.hpp"
 #include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
 #include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
 #include "physicsSolvers/fluidFlow/kernels/compositional/RelativePermeabilityUpdateKernel.hpp"
 #include "physicsSolvers/fluidFlow/kernels/compositional/CapillaryPressureUpdateKernel.hpp"
 #include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
 #include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/KernelLaunchSelectors.hpp"
 #include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/KernelLaunchSelectors.hpp"


namespace geos
{
namespace immiscibleMultiphaseKernels
{

using namespace constitutive;

GEOS_HOST_DEVICE
inline
static void local_solver( real64 uT, stdVector< real64 > saturations, stdVector< real64 > pressures, stdVector< real64 > JFMultipliers, stdVector< real64 > trappedSats1,
                          stdVector< real64 > trappedSats2, stdVector< real64 > transHat, stdVector< real64 > dTransHat_dP, stdVector< real64 > gravCoefHat, stdVector< real64 > gravCoef,
                          stdVector< real64 >  cellCenterDuT, stdVector< real64 > cellCenterDens, stdVector< real64 > cellCenterDens_dP,
                          std::vector< RelativePermeabilityBase * > relPerms, std::vector< CapillaryPressureBase * > capPressures,
                          std::vector< TwoPhaseImmiscibleFluid * > fluids, std::vector< real64 > &phi, std::vector< real64 > &grad_phi_P, std::vector< real64 > &grad_phi_S, bool &converged )
{

  // getting wrappers:

  constitutive::constitutiveUpdatePassThru( *capPressures[0], [&] ( auto & castedCapPres1 )
  {
    auto capPresWrapper1 = castedCapPres1.createKernelWrapper();


    constitutive::constitutiveUpdatePassThru( *capPressures[1], [&] ( auto & castedCapPres2 )
    {
      auto capPresWrapper2 = castedCapPres2.createKernelWrapper();


      constitutive::constitutiveUpdatePassThru( *relPerms[0], [&] ( auto & castedRelPerm1 )
      {
        auto relPermWrapper1 = castedRelPerm1.createKernelWrapper();


        constitutive::constitutiveUpdatePassThru( *relPerms[1], [&] ( auto & castedRelPerm2 )
        {
          auto relPermWrapper2 = castedRelPerm2.createKernelWrapper();

          auto fluidWrapper1 = fluids[0]->createKernelWrapper();
          auto fluidWrapper2 = fluids[1]->createKernelWrapper();

          // Create an output file stream object (ofstream) for analyzing the local solver's performance
          std::ofstream outFile( "iterations2.csv" );


          // Write data to the file
          outFile << "Jacobian";
          outFile << ",";
          outFile << "residual";
          outFile << ",";
          outFile << "Fw_alpha";
          outFile << ",";
          outFile << "Fw_beta";
          outFile << ",";
          outFile << "Pc_int";
          outFile << ",";
          outFile << "Pc_int1";
          outFile << ",";
          outFile << "Pc_int2";
          outFile << ",";
          outFile << "Fn_alpha";
          outFile << ",";
          outFile << "Fn_beta";
          outFile << ",";
          outFile << "Vw_alpha";
          outFile << ",";
          outFile << "Vn_alpha";
          outFile << ",";
          outFile << "Vw_beta";
          outFile << ",";
          outFile << "Vn_beta";
          outFile << ",";
          outFile << "Gw_alpha";
          outFile << ",";
          outFile << "Gn_alpha";
          outFile << ",";
          outFile << "Gw_beta";
          outFile << ",";
          outFile << "Gn_beta";
          outFile << ",";
          outFile << "Cw_alpha";
          outFile << ",";
          outFile << "Cn_alpha";
          outFile << ",";
          outFile << "Cw_beta";
          outFile << ",";
          outFile << "Cn_beta";
          outFile << ",";
          outFile << "Pc1_ip0";
          outFile << ",";
          outFile << "Pc1_ip1";
          outFile << ",";
          outFile << "S_alpha";
          outFile << ",";
          outFile << "S_beta";
          outFile << std::endl;

          // nonlinear solver's parameters
          real64 tol = 1.0e-9;
          int max_iter = 50;
          converged = 0;
          bool damping = true;
          bool bisection = false;
          // bool lineSearch = false;
          bool newton_path = false;

          // Local newton loop:

          // Use of the capillary pressure kernel wrapper

          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseVolFrac1( 1, 2 );
          StackArray< real64, 3, 2, constitutive::cappres::LAYOUT_CAPPRES > capPres1( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dPhaseVolFrac_dCapPres1( 1, 1, 2, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dCapPres1_dPhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > facePhaseVolFrac1( 1, 2 );
          StackArray< real64, 3, 2, constitutive::cappres::LAYOUT_CAPPRES > faceCapPres1( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dfacePhaseVolFrac_dCapPres1( 1, 1, 2, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dCapPres1_dfacePhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseVolFrac2( 1, 2 );
          StackArray< real64, 3, 2, constitutive::cappres::LAYOUT_CAPPRES > capPres2( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dPhaseVolFrac_dCapPres2( 1, 1, 2, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dCapPres2_dPhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > facePhaseVolFrac2( 1, 2 );
          StackArray< real64, 3, 2, constitutive::cappres::LAYOUT_CAPPRES > faceCapPres2( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dfacePhaseVolFrac_dCapPres2( 1, 1, 2, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dCapPres2_dfacePhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 1, 2 > JFunc1( 2 );
          StackArray< real64, 1, 2 > JFunc2( 2 );

          phaseVolFrac1[0][0] = saturations[0];
          phaseVolFrac1[0][1] = 1.0 - saturations[0];

          phaseVolFrac2[0][0] = saturations[1];
          phaseVolFrac2[0][1] = 1.0 - saturations[1];

          real64 Pc1_min = 0.0;
          real64 Pc2_min = 0.0;
          real64 Pc1_max = 0.0;
          real64 Pc2_max = 0.0;

          real64 density2[2]{};
          real64 dDens_dP2[2][2]{};

          density2[0] =  cellCenterDens[0];  
          density2[1] =  cellCenterDens[1];         

          dDens_dP2[0][0] = cellCenterDens_dP[0];
          dDens_dP2[0][1] = cellCenterDens_dP[1];
          dDens_dP2[1][0] = cellCenterDens_dP[2];
          dDens_dP2[1][1] = cellCenterDens_dP[3];

          JFunc1[0] = JFMultipliers[0];
          JFunc2[0] = JFMultipliers[1];

          using T1 = std::decay_t< decltype(castedCapPres1) >;
          if constexpr (std::is_same_v< T1, JFunctionCapillaryPressure >) {
            capPresWrapper1.compute( phaseVolFrac1[0],
                                     JFunc1.toSliceConst(),
                                     capPres1[0][0],
                                     dCapPres1_dPhaseVolFrac[0][0] );

            facePhaseVolFrac1[0][1] = 0.0;
            facePhaseVolFrac1[0][0] = 1.0;
            capPresWrapper1.compute( facePhaseVolFrac1[0],
                                     JFunc1.toSliceConst(),
                                     faceCapPres1[0][0],
                                     dCapPres1_dfacePhaseVolFrac[0][0] );
            Pc1_min = faceCapPres1[0][0][0];
            facePhaseVolFrac1[0][1] = 1.0;
            facePhaseVolFrac1[0][0] = 0.0;
            capPresWrapper1.compute( facePhaseVolFrac1[0],
                                     JFunc1.toSliceConst(),
                                     faceCapPres1[0][0],
                                     dCapPres1_dfacePhaseVolFrac[0][0] );
            Pc1_max = faceCapPres1[0][0][0];

          }
          else
          {
            capPresWrapper1.compute( phaseVolFrac1[0],
                                     capPres1[0][0],
                                     dCapPres1_dPhaseVolFrac[0][0] );

            facePhaseVolFrac1[0][1] = 0.0;
            facePhaseVolFrac1[0][0] = 1.0;
            capPresWrapper1.compute( facePhaseVolFrac1[0],
                                     faceCapPres1[0][0],
                                     dCapPres1_dfacePhaseVolFrac[0][0] );
            Pc1_min = faceCapPres1[0][0][0];
            facePhaseVolFrac1[0][1] = 1.0;
            facePhaseVolFrac1[0][0] = 0.0;
            capPresWrapper1.compute( facePhaseVolFrac1[0],
                                     faceCapPres1[0][0],
                                     dCapPres1_dfacePhaseVolFrac[0][0] );
            Pc1_max = faceCapPres1[0][0][0];
          }

          using T2 = std::decay_t< decltype(castedCapPres2) >;
          if constexpr (std::is_same_v< T2, JFunctionCapillaryPressure >) {
            // evaluating cell-center Pc:

            capPresWrapper2.compute( phaseVolFrac2[0],
                                     JFunc2.toSliceConst(),
                                     capPres2[0][0],
                                     dCapPres2_dPhaseVolFrac[0][0] );

// finding endpoints:

            facePhaseVolFrac2[0][1] = 0.0;
            facePhaseVolFrac2[0][0] = 1.0;
            capPresWrapper2.compute( facePhaseVolFrac2[0],
                                     JFunc2.toSliceConst(),
                                     faceCapPres2[0][0],
                                     dCapPres2_dfacePhaseVolFrac[0][0] );
            Pc2_min = faceCapPres2[0][0][0];
            facePhaseVolFrac2[0][1] = 1.0;
            facePhaseVolFrac2[0][0] = 0.0;
            capPresWrapper2.compute( facePhaseVolFrac2[0],
                                     JFunc2.toSliceConst(),
                                     faceCapPres2[0][0],
                                     dCapPres2_dfacePhaseVolFrac[0][0] );
            Pc2_max = faceCapPres2[0][0][0];

          }
          else
          {
            // evaluating cell-center Pc:

            capPresWrapper2.compute( phaseVolFrac2[0],
                                     capPres2[0][0],
                                     dCapPres2_dPhaseVolFrac[0][0] );

// finding endpoints:

            facePhaseVolFrac2[0][1] = 0.0;
            facePhaseVolFrac2[0][0] = 1.0;
            capPresWrapper2.compute( facePhaseVolFrac2[0],
                                     faceCapPres2[0][0],
                                     dCapPres2_dfacePhaseVolFrac[0][0] );
            Pc2_min = faceCapPres2[0][0][0];
            facePhaseVolFrac2[0][1] = 1.0;
            facePhaseVolFrac2[0][0] = 0.0;
            capPresWrapper2.compute( facePhaseVolFrac2[0],
                                     faceCapPres2[0][0],
                                     dCapPres2_dfacePhaseVolFrac[0][0] );
            Pc2_max = faceCapPres2[0][0][0];
          }

          // Use of the relative permeability kernel wrapper

          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > faceTrappedVolFrac1( 1, 1, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > faceRelPerm1( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::relperm::LAYOUT_RELPERM_DS > dfacePhaseRelPerm1_dPhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > trappedVolFrac1( 1, 1, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > relPerm1( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::relperm::LAYOUT_RELPERM_DS > dPhaseRelPerm1_dPhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > faceTrappedVolFrac2( 1, 1, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > faceRelPerm2( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::relperm::LAYOUT_RELPERM_DS > dfacePhaseRelPerm2_dPhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > trappedVolFrac2( 1, 1, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > relPerm2( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::relperm::LAYOUT_RELPERM_DS > dPhaseRelPerm2_dPhaseVolFrac( 1, 1, 2, 2 );

          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseMaxHistoricalVolFraction1( 1, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseMinHistoricalVolFraction1( 1, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseMaxHistoricalVolFraction2( 1, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseMinHistoricalVolFraction2( 1, 2 );

          // compute relative permeability for both cell centers:

          trappedVolFrac1[0][0][0] = trappedSats1[0];
          trappedVolFrac1[0][0][1] = trappedSats1[1];

          trappedVolFrac2[0][0][0] = trappedSats2[0];
          trappedVolFrac2[0][0][1] = trappedSats2[1];

          faceTrappedVolFrac1[0][0][0] = trappedSats1[0];
          faceTrappedVolFrac1[0][0][1] = trappedSats1[1];

          faceTrappedVolFrac2[0][0][0] = trappedSats2[0];
          faceTrappedVolFrac2[0][0][1] = trappedSats2[1];

          using T5 = std::decay_t< decltype(castedRelPerm1) >;
          if constexpr (std::is_same_v< T5, constitutive::TableRelativePermeabilityHysteresis >) {
            relPermWrapper1.compute( phaseVolFrac1[0],
                                     phaseMaxHistoricalVolFraction1[0],
                                     phaseMinHistoricalVolFraction1[0],
                                     trappedVolFrac1[0][0],
                                     relPerm1[0][0],
                                     dPhaseRelPerm1_dPhaseVolFrac[0][0] );

          }
          else
          {

            relPermWrapper1.compute( phaseVolFrac1[0],
                                     trappedVolFrac1[0][0],
                                     relPerm1[0][0],
                                     dPhaseRelPerm1_dPhaseVolFrac[0][0] );
          }

          using T6 = std::decay_t< decltype(castedRelPerm2) >;
          if constexpr (std::is_same_v< T6, constitutive::TableRelativePermeabilityHysteresis >) {
            relPermWrapper2.compute( phaseVolFrac2[0],
                                     phaseMaxHistoricalVolFraction2[0],
                                     phaseMinHistoricalVolFraction2[0],
                                     trappedVolFrac2[0][0],
                                     relPerm2[0][0],
                                     dPhaseRelPerm2_dPhaseVolFrac[0][0] );

          }
          else
          {

            relPermWrapper2.compute( phaseVolFrac2[0],
                                     trappedVolFrac2[0][0],
                                     relPerm2[0][0],
                                     dPhaseRelPerm2_dPhaseVolFrac[0][0] );
          }


          // Use of the fluid model kernel wrapper
          StackArray< real64, 3, 2, constitutive::multifluid::LAYOUT_PHASE > phaseDensity1( 1, 1, 2 );
          StackArray< real64, 3, 2, constitutive::multifluid::LAYOUT_PHASE > phaseViscosity1( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::multifluid::LAYOUT_PHASE_DC > dPhaseDens1_dP( 1, 1, 2, 2 );
          StackArray< real64, 4, 4, constitutive::multifluid::LAYOUT_PHASE_DC > dPhaseVisc1_dP( 1, 1, 2, 2 );
          StackArray< real64, 3, 2, constitutive::multifluid::LAYOUT_PHASE > phaseDensity2( 1, 1, 2 );
          StackArray< real64, 3, 2, constitutive::multifluid::LAYOUT_PHASE > phaseViscosity2( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::multifluid::LAYOUT_PHASE_DC > dPhaseDens2_dP( 1, 1, 2, 2 );
          StackArray< real64, 4, 4, constitutive::multifluid::LAYOUT_PHASE_DC > dPhaseVisc2_dP( 1, 1, 2, 2 );

          // Declare the MultiFluidVar (PhaseProp) type
          TwoPhaseImmiscibleFluid::PhaseProp phaseDensity_temp1;
          TwoPhaseImmiscibleFluid::PhaseProp phaseViscosity_temp1;
          phaseDensity_temp1.value.resize( 1, 1, 2 ); // or whatever sizes you need
          phaseDensity_temp1.derivs.resize( 1, 1, 2, 2 ); // make sure all dims > 0
          phaseViscosity_temp1.value.resize( 1, 1, 2 ); // or whatever sizes you need
          phaseViscosity_temp1.derivs.resize( 1, 1, 2, 2 ); // make sure all dims > 0

          auto phaseDensitySlice0 = TwoPhaseImmiscibleFluid::PhaseProp::SliceType(
            phaseDensity_temp1.value[0][0], phaseDensity_temp1.derivs[0][0] );
          auto phaseViscositySlice0 = TwoPhaseImmiscibleFluid::PhaseProp::SliceType(
            phaseViscosity_temp1.value[0][0], phaseViscosity_temp1.derivs[0][0] );
          fluidWrapper1.compute( pressures[0], phaseDensitySlice0, phaseViscositySlice0 );

          // Temporary:
          phaseDensity1[0][0][0] = phaseDensitySlice0.value[0];
          phaseDensity1[0][0][1] = phaseDensitySlice0.value[1];
          dPhaseDens1_dP[0][0][0][0] = phaseDensitySlice0.derivs[0][0];
          dPhaseDens1_dP[0][0][0][1] = 0;
          dPhaseDens1_dP[0][0][1][0] = phaseDensitySlice0.derivs[1][0];
          dPhaseDens1_dP[0][0][1][1] = 0;

          phaseViscosity1[0][0][0] = phaseViscositySlice0.value[0];
          phaseViscosity1[0][0][1] = phaseViscositySlice0.value[1];
          dPhaseVisc1_dP[0][0][0][0] = phaseViscositySlice0.derivs[0][0];
          dPhaseVisc1_dP[0][0][0][1] = 0;
          dPhaseVisc1_dP[0][0][1][0] = phaseViscositySlice0.derivs[1][0];
          dPhaseVisc1_dP[0][0][1][1] = 0;

          auto phaseDensitySlice1 = TwoPhaseImmiscibleFluid::PhaseProp::SliceType(
            phaseDensity_temp1.value[0][0], phaseDensity_temp1.derivs[0][0] );
          auto phaseViscositySlice1 = TwoPhaseImmiscibleFluid::PhaseProp::SliceType(
            phaseViscosity_temp1.value[0][0], phaseViscosity_temp1.derivs[0][0] );
          fluidWrapper2.compute( pressures[1], phaseDensitySlice1, phaseViscositySlice1 );


          phaseDensity2[0][0][0] = phaseDensitySlice1.value[0];
          phaseDensity2[0][0][1] = phaseDensitySlice1.value[1];
          dPhaseDens2_dP[0][0][0][0] = phaseDensitySlice1.derivs[0][0];
          dPhaseDens2_dP[0][0][0][1] = 0;
          dPhaseDens2_dP[0][0][1][0] = phaseDensitySlice1.derivs[1][0];
          dPhaseDens2_dP[0][0][1][1] = 0;

          phaseViscosity2[0][0][0] = phaseViscositySlice1.value[0];
          phaseViscosity2[0][0][1] = phaseViscositySlice1.value[1];
          dPhaseVisc2_dP[0][0][0][0] = phaseViscositySlice1.derivs[0][0];
          dPhaseVisc2_dP[0][0][0][1] = 0;
          dPhaseVisc2_dP[0][0][1][0] = phaseViscositySlice1.derivs[1][0];
          dPhaseVisc2_dP[0][0][1][1] = 0;

          // clear working arrays
          real64 halfFluxVal[2][2]{};
          real64 dhalfFlux1_dP[2][2]{};
          real64 dhalfFlux1_dS[2][2]{};
          real64 dhalfFlux2_dP[2][2]{};
          real64 dhalfFlux2_dS[2][2]{};
          real64 dhalfFlux_duT[2][2]{};
          real64 dhalfFlux_dpc[2][2]{};

          //new
        real64 fluxVal[2]{};
        real64 dFlux_dP[2][2]{};
        real64 dFlux_dS[2][2]{};

        real64 duT_dP[2]{};
        real64 duT_dS[2]{};

        duT_dP[0] = cellCenterDuT[0];
        duT_dP[1] = cellCenterDuT[1];

        duT_dS[0] = cellCenterDuT[2];
        duT_dS[1] = cellCenterDuT[3];

          // initial guess:

          real64 const Pc1 = capPres1[0][0][0];
          real64 const Pc2 = capPres2[0][0][0];

          real64 Pc_int = ( Pc1 + Pc2 ) / 2.0;

          real64 Pc_min_all = fmax( Pc1_min, Pc2_min );
          real64 Pc_max_all = fmin( Pc1_max, Pc2_max );

          if( Pc_int < Pc_min_all || Pc_int > Pc_max_all )
          {
            Pc_int = ( Pc_min_all + Pc_max_all ) / 2.0;
          }

          // While loop (newton loop)
          int iter = 0;
          int div = 0;
          // int ext_iter0 = 0;
          // int ext_iter1 = 0;
          real64 next_Pc_int = 0.0;
          real64 old_Pc_int = 0.0;
          real64 old_residual = 0.0;
          
          if (bisection) { 
            Pc_int = fmax( Pc1_max, Pc2_max );
            next_Pc_int = fmin( Pc1_min, Pc2_min );
          }

          if (newton_path){ 
            Pc_int = fmax( Pc1_max, Pc2_max );
             Pc_int = 2.0e5;
             next_Pc_int = (fmax( Pc1_max, Pc2_max ) - fmin( Pc1_min, Pc2_min )) / (max_iter - 1);
             next_Pc_int = (2.0e5 - 5.0e4) / (max_iter - 1);
          }
  
          real64 Pc_int_iterate = Pc_int;

          while( iter < max_iter )
          {

            Pc_int_iterate = Pc_int;

            // clear working arrays
            real64 density[2]{};
            real64 dDens_dP[2][2]{};
            real64 gravityCof[2]{};
            real64 viscosity[2]{};
            real64 dVisc_dP[2][2]{};

            real64 viscous[2][2]{};
            real64 bouyancy[2][2]{};
            real64 capillarity[2][2]{};

            real64 dV1_dS[2][2]{};
            real64 dG1_dS[2][2]{};
            real64 dC1_dS[2][2]{};
            real64 dV2_dS[2][2]{};
            real64 dG2_dS[2][2]{};
            real64 dC2_dS[2][2]{};
            real64 dV1_dpc[2][2]{};
            real64 dG1_dpc[2][2]{};
            real64 dC1_dpc[2][2]{};
            real64 dV2_dpc[2][2]{};
            real64 dG2_dpc[2][2]{};
            real64 dC2_dpc[2][2]{};

            real64 local_residual = 0;
            real64 local_jacobian = 0;

            // truncate the capillary pressure iterate (ensures the inverse will compute a saturation bounded between 0 and 1):

            Pc_int = fmin( fmax( Pc1_max, Pc2_max ), fmax( Pc_int, fmin( Pc1_min, Pc2_min ) ));
            
            faceCapPres1[0][0][0] = fmin( Pc1_max, fmax( Pc_int, Pc1_min));
            faceCapPres2[0][0][0] = fmin( Pc2_max, fmax( Pc_int, Pc2_min));

            // Compute the inverse:

            JFunc1[0] = JFMultipliers[0];
            JFunc2[0] = JFMultipliers[1];

            using T3 = std::decay_t< decltype(castedCapPres1) >;
            if constexpr (std::is_same_v< T3, JFunctionCapillaryPressure >) {
              capPresWrapper1.computeInv( facePhaseVolFrac1[0],
                                          JFunc1.toSliceConst(),
                                          faceCapPres1[0][0],
                                          dfacePhaseVolFrac_dCapPres1[0][0] );
              facePhaseVolFrac1[0][0] = fmin( 1.0, fmax( facePhaseVolFrac1[0][0], 0.0 ));
              facePhaseVolFrac1[0][1] = fmin( 1.0, fmax( facePhaseVolFrac1[0][1], 0.0 ));
              //get derivatives:
              capPresWrapper1.compute( facePhaseVolFrac1[0],
                                       JFunc1.toSliceConst(),
                                       faceCapPres1[0][0],
                                       dCapPres1_dfacePhaseVolFrac[0][0] );

            }
            else
            {
              capPresWrapper1.computeInv( facePhaseVolFrac1[0],
                                       faceCapPres1[0][0],
                                       dfacePhaseVolFrac_dCapPres1[0][0] );
              facePhaseVolFrac1[0][0] = fmin( 1.0, fmax( facePhaseVolFrac1[0][0], 0.0 ));
              facePhaseVolFrac1[0][1] = fmin( 1.0, fmax( facePhaseVolFrac1[0][1], 0.0 ));
              //get derivatives:
              capPresWrapper1.compute( facePhaseVolFrac1[0],
                                       faceCapPres1[0][0],
                                       dCapPres1_dfacePhaseVolFrac[0][0] );
            }

            using T4 = std::decay_t< decltype(castedCapPres2) >;
            if constexpr (std::is_same_v< T4, JFunctionCapillaryPressure >) {
              // evaluating cell-center Pc:
              capPresWrapper2.computeInv( facePhaseVolFrac2[0],
                                          JFunc2.toSliceConst(),
                                          faceCapPres2[0][0],
                                          dfacePhaseVolFrac_dCapPres2[0][0] );
              facePhaseVolFrac2[0][0] = fmin( 1.0, fmax( facePhaseVolFrac2[0][0], 0.0 ));
              facePhaseVolFrac2[0][1] = fmin( 1.0, fmax( facePhaseVolFrac2[0][1], 0.0 ));

//get derivatives:

              capPresWrapper2.compute( facePhaseVolFrac2[0],
                                       JFunc2.toSliceConst(),
                                       faceCapPres2[0][0],
                                       dCapPres2_dfacePhaseVolFrac[0][0] );

            }
            else
            {
              // evaluating cell-center Pc:
              capPresWrapper2.computeInv( facePhaseVolFrac2[0],
                                       faceCapPres2[0][0],
                                       dfacePhaseVolFrac_dCapPres2[0][0] );
              facePhaseVolFrac2[0][0] = fmin( 1.0, fmax( facePhaseVolFrac2[0][0], 0.0 ));
              facePhaseVolFrac2[0][1] = fmin( 1.0, fmax( facePhaseVolFrac2[0][1], 0.0 ));

//get derivatives:

              capPresWrapper2.compute( facePhaseVolFrac2[0],
                                       faceCapPres2[0][0],
                                       dCapPres2_dfacePhaseVolFrac[0][0] );
            }
            
            // compute relative permeability for both faces:

            using T7 = std::decay_t< decltype(castedRelPerm1) >;
            if constexpr (std::is_same_v< T7, constitutive::TableRelativePermeabilityHysteresis >) {

              relPermWrapper1.compute( facePhaseVolFrac1[0],
                                       phaseMaxHistoricalVolFraction1[0],
                                       phaseMinHistoricalVolFraction1[0],
                                       faceTrappedVolFrac1[0][0],
                                       faceRelPerm1[0][0],
                                       dfacePhaseRelPerm1_dPhaseVolFrac[0][0] );

            }
            else
            {

              relPermWrapper1.compute( facePhaseVolFrac1[0],
                                       faceTrappedVolFrac1[0][0],
                                       faceRelPerm1[0][0],
                                       dfacePhaseRelPerm1_dPhaseVolFrac[0][0] );
            }

            using T8 = std::decay_t< decltype(castedRelPerm2) >;
            if constexpr (std::is_same_v< T8, constitutive::TableRelativePermeabilityHysteresis >) {

              relPermWrapper2.compute( facePhaseVolFrac2[0],
                                       phaseMaxHistoricalVolFraction2[0],
                                       phaseMinHistoricalVolFraction2[0],
                                       faceTrappedVolFrac2[0][0],
                                       faceRelPerm2[0][0],
                                       dfacePhaseRelPerm2_dPhaseVolFrac[0][0] );

            }
            else
            {

              relPermWrapper2.compute( facePhaseVolFrac2[0],
                                       faceTrappedVolFrac2[0][0],
                                       faceRelPerm2[0][0],
                                       dfacePhaseRelPerm2_dPhaseVolFrac[0][0] );
            }


            // Brenier and Jaffre's PPU upwinding:
            bool check = false;
            bool k_up_0_w = 1;
            bool k_up_1_w = 1;
            bool k_up_0_n = 1;
            bool k_up_1_n = 1;

            if( uT < 0.0 )
            {
              k_up_0_w = 0;
              k_up_1_w = 0;
              k_up_0_n = 0;
              k_up_1_n = 0;
            }

            if( std::fabs( uT ) < 1e-20 )
            {
              k_up_0_w = 1;
              k_up_1_w = 1;
              k_up_0_n = 0;
              k_up_1_n = 0;
            }


            localIndex k_up_0[2] = {static_cast< localIndex >(!k_up_0_w), static_cast< localIndex >(!k_up_0_n)};
            localIndex k_up_1[2] = {static_cast< localIndex >(k_up_1_w), static_cast< localIndex >(k_up_1_n)};
            bool k_up_0_check[2] = {false, false};
            bool k_up_1_check[2] = {false, false};

            while( !check )
            {
              k_up_0_check[0] = false;
              k_up_0_check[1] = false;
              k_up_1_check[0] = false;
              k_up_1_check[1] = false;
              for( integer ix = 0; ix < 2; ++ix ) // for loop over each half flux
              {

                // clear working arrays for each half flux:
                real64 densMean[2]{};
                real64 dDensMean_dP[2][2]{};

                real64 presGrad[2]{};
                real64 dPresGrad_dP[2][2]{};

                real64 gravHead[2]{};
                real64 dGravHead_dP[2][2]{};

                real64 capGrad[2]{};
                // real64 capPresIC[2][2]{};
                // real64 jFMultiplier[2][2]{};
                real64 dCapGrad_dP[2][2]{};
                real64 dCapGrad_dS[2][2]{};

                real64 mobility[2]{};
                real64 dMob_dP[2][2]{};
                real64 dMob_dS[2][2]{};

                real64 total_mobility = 0;
                gravityCof[0] = 0;
                gravityCof[1] = 0;

                for( integer ip = 0; ip < 2; ++ip ) // loop over phases
                {
                  // calculate quantities on primary connected cells
                  if( ix == 0 )
                  {
                    density[ip]  = phaseDensity1[0][0][ip];
                    dDens_dP[ip][ix] = dPhaseDens1_dP[0][0][ip][ip];

                    viscosity[ip]  = phaseViscosity1[0][0][ip];
                    dVisc_dP[ip][ix] = dPhaseVisc1_dP[0][0][ip][ip];
                  }
                  else
                  {
                    density[ip]  = phaseDensity2[0][0][ip];
                    dDens_dP[ip][ix] = dPhaseDens2_dP[0][0][ip][ip];

                    viscosity[ip]  = phaseViscosity2[0][0][ip];
                    dVisc_dP[ip][ix] = dPhaseVisc2_dP[0][0][ip][ip];
                  }

                  densMean[ip] = density[ip];
                  dDensMean_dP[ip][0] = dDens_dP[ip][ix];
                  dDensMean_dP[ip][1] = dDens_dP[ip][ix];

                  //***** calculation of flux *****

                  // compute potential difference
                  real64 potScale = 0.0;
                  real64 dPresGrad_dTrans = 0.0;
                  real64 dGravHead_dTrans = 0.0;
                  real64 dCapGrad_dTrans = 0.0;
                  constexpr int signPotDiff[2] = {1, -1};
                  constexpr int signTix[2] = {1, -1};

                  for( integer ke = 0; ke < 2; ++ke )
                  {

                    real64 const pressure = pressures[ix];
                    presGrad[ip] += signTix[ke] * transHat[ix] * pressure;
                    dPresGrad_dTrans += signPotDiff[ke] * pressure;
                    dPresGrad_dP[ip][ke] = signTix[ke] * transHat[ix];

                    real64 gravD = 0.0;

                    if( ke == 0 )
                    {
                      gravD += signTix[ke] * transHat[ix] * gravCoef[ix];
                    }
                    else
                    {
                      gravD += signTix[ke] * transHat[ix] * gravCoefHat[ix];
                    }

                    real64 pot = signTix[ke] * transHat[ix] * pressure - densMean[ip] * gravD;

                    gravHead[ip] += densMean[ip] * gravD;
                    gravityCof[ip] += gravD;

                    dGravHead_dTrans += signPotDiff[ke] * densMean[ip] * gravCoefHat[ix];

                    for( integer i = 0; i < 2; ++i )
                    {
                      dGravHead_dP[ip][i] += dDensMean_dP[ip][i] * gravD;
                    }

                    real64 capPres = capPres1[0][0][ip];

                    if( ke == 1 && ix == 0 )
                    {
                      capPres = faceCapPres1[0][0][ip];
                    }
                    else if( ke == 1 && ix == 1 )
                    {
                      capPres = faceCapPres2[0][0][ip];
                    }
                    else if( ke == 0 && ix == 1 )
                    {

                      capPres = capPres2[0][0][ip];
                    }

                    dCapGrad_dTrans -= signPotDiff[ke] * capPres;
                    pot -= signTix[ke] * transHat[ix] * capPres;

                    capGrad[ip] -= signTix[ke] * transHat[ix] * capPres;

                    potScale = fmax( potScale, fabs( pot ) );
                  }

                  for( integer ke = 0; ke < 2; ++ke )
                  {
                    dPresGrad_dP[ip][ke] += signTix[ke] * dTransHat_dP[ix] * dPresGrad_dTrans;

                    dGravHead_dP[ip][ke] += signTix[ke] * dTransHat_dP[ix] * dGravHead_dTrans;

                    // real64 constexpr eps = 1e-18;
                    real64 dCapPres_dS = dCapPres1_dPhaseVolFrac[0][0][ip][ip];

                    if( ke == 1 && ix == 0 )
                    {
                      dCapPres_dS = dCapPres1_dfacePhaseVolFrac[0][0][ip][ip];
                    }
                    else if( ke == 1 && ix == 1 )
                    {
                      dCapPres_dS = dCapPres2_dfacePhaseVolFrac[0][0][ip][ip];
                    }
                    else if( ke == 0 && ix == 1 )
                    {
                      dCapPres_dS = dCapPres2_dPhaseVolFrac[0][0][ip][ip];
                    }

                    dCapGrad_dP[ip][ke] += signTix[ke] * dTransHat_dP[ix] * dCapGrad_dTrans;
                    dCapGrad_dS[ip][ke] -= signTix[ke] * transHat[ix] * dCapPres_dS;
                  }

                  // *** upwinding ***
                  // compute potential gradient
                  real64 potGrad = presGrad[ip] - gravHead[ip];

                  potGrad += capGrad[ip];

                  // choose upstream cell
                  constexpr int sign[2] = {1, -1};

                  if( k_up_0[ip] == 1 && ix == 0 )
                  {
                    mobility[ip] = faceRelPerm1[0][0][ip] / viscosity[ip];
                    dMob_dP[ip][k_up_0[ip]] = mobility[ip] * (-dVisc_dP[ip][ix] / viscosity[ip]);

                    dMob_dS[ip][k_up_0[ip]] = sign[ip] * dfacePhaseRelPerm1_dPhaseVolFrac[0][0][ip][ip] / viscosity[ip];
                  }
                  else if( k_up_1[ip] == 0 && ix == 1 )
                  {
                    mobility[ip] = relPerm2[0][0][ip] / viscosity[ip];
                    dMob_dP[ip][0] = mobility[ip] * (-dVisc_dP[ip][ix] / viscosity[ip]);
                    dMob_dS[ip][0] = sign[ip] * dPhaseRelPerm2_dPhaseVolFrac[0][0][ip][ip] / viscosity[ip];
                  }
                  else if( k_up_1[ip] == 1 && ix == 1 )
                  {
                    mobility[ip] = faceRelPerm2[0][0][ip] / viscosity[ip];
                    dMob_dP[ip][1] = mobility[ip] * (-dVisc_dP[ip][ix] / viscosity[ip]);
                    dMob_dS[ip][1] = sign[ip] * dfacePhaseRelPerm2_dPhaseVolFrac[0][0][ip][ip] / viscosity[ip];
                  }
                  else
                  {
                    mobility[ip] = relPerm1[0][0][ip] / viscosity[ip];
                    dMob_dP[ip][ix] = mobility[ip] * (-dVisc_dP[ip][ix] / viscosity[ip]);
                    dMob_dS[ip][ix] = sign[ip] * dPhaseRelPerm1_dPhaseVolFrac[0][0][ip][ip] / viscosity[ip];
                  }
                  real64 constexpr eps = 0.0;
                  total_mobility += mobility[ip] + eps;
                } // loop over phases

                /// Three Forces Flux Contribution: 1- Viscous 2- Gravitational 3- Capillary
                constexpr int sign[2] = {1, -1};
                // real64 constexpr eps = 0.0;

                // loop over phases
                for( integer ip = 0; ip < 2; ++ip )
                {
                  // 1- Viscous: pressure gradient depends on all points in the stencil
                  viscous[ip][ix] = mobility[ip] / total_mobility * uT;
                  halfFluxVal[ip][ix] = viscous[ip][ix];

                  for( integer ke = 0; ke < 2; ++ke )
                  {
                    real64 dV_dP = sign[ip] * (dMob_dP[0][ke] * mobility[1] - dMob_dP[1][ke] * mobility[0]) / (total_mobility * total_mobility) * uT;
                    real64 dV_dS = sign[ip] * (dMob_dS[0][ke] * mobility[1] - dMob_dS[1][ke] * mobility[0]) / (total_mobility * total_mobility) * uT;
                    real64 dV_du = mobility[ip] / total_mobility;

                    if( ix == 0 )
                    {
                      dV1_dS[ip][ke] = dV_dS;
                      dhalfFlux1_dP[ip][ke] = dV_dP;
                      dhalfFlux1_dS[ip][ke] = dV1_dS[ip][ke];
                      dhalfFlux_duT[ip][ix] = dV_du;
                      dV1_dpc[ip][ke] = dV1_dS[ip][ke] * dfacePhaseVolFrac_dCapPres1[0][0][0][0];
                      // GEOS_UNUSED_VAR( dV1_dpc[ip][ke] );
                      GEOS_UNUSED_VAR( dhalfFlux_duT[ip][ix] );
                    }
                    else
                    {
                      dV2_dS[ip][ke] = dV_dS;
                      dhalfFlux2_dP[ip][ke] =  dV_dP;
                      dhalfFlux2_dS[ip][ke] =  dV2_dS[ip][ke];
                      dhalfFlux_duT[ip][ix]= dV_du;
                      dV2_dpc[ip][ke] = dV2_dS[ip][ke] * dfacePhaseVolFrac_dCapPres2[0][0][0][0];
                      // GEOS_UNUSED_VAR( dV2_dpc[ip][ke] );
                      GEOS_UNUSED_VAR( dhalfFlux_duT[ip][ix] );
                    }
                  }

                  // 2- Gravitational: gravitational head depends only on the two cells connected (same as mean density)
                  bouyancy[ip][ix] = -1.0 * sign[ix] * sign[ip] * mobility[0] * mobility[1] / total_mobility * gravityCof[0] * (density[0] - density[1]);
                  halfFluxVal[ip][ix] += bouyancy[ip][ix];

                  for( integer ke = 0; ke < 2; ++ke )
                  {
                    real64 dG_dP = sign[ix] * sign[ip] * (dMob_dP[0][ke] * mobility[1] * mobility[1] + dMob_dP[1][ke] * mobility[0] * mobility[0]) / (total_mobility * total_mobility) * gravityCof[0] *
                                   (density[0] - density[1]) +
                                   sign[ix] * (mobility[0] * mobility[1]) / total_mobility * (dDens_dP[0][ix] - dDens_dP[1][ix]);

                    real64 dG_dS =  sign[ix] * sign[ip] * (dMob_dS[0][ke] * mobility[1] * mobility[1] + dMob_dS[1][ke] * mobility[0] * mobility[0]) / (total_mobility * total_mobility) *
                                   gravityCof[0] * (density[0] - density[1]);

                    if( ix == 0 )
                    {
                      dG1_dS[ip][ke] = dG_dS;
                      dhalfFlux1_dP[ip][ke] -= dG_dP;
                      dhalfFlux1_dS[ip][ke] -= dG1_dS[ip][ke];
                      dG1_dpc[ip][ke] = dG1_dS[ip][ke] * dfacePhaseVolFrac_dCapPres1[0][0][0][0];
                      // GEOS_UNUSED_VAR( dG1_dpc[ip][ke] );
                    }
                    else
                    {
                      dG2_dS[ip][ke] = dG_dS;
                      dhalfFlux2_dP[ip][ke] -= dG_dP;
                      dhalfFlux2_dS[ip][ke] -= dG2_dS[ip][ke];
                      dG2_dpc[ip][ke] = dG2_dS[ip][ke] * dfacePhaseVolFrac_dCapPres2[0][0][0][0];
                      // GEOS_UNUSED_VAR( dG2_dpc[ip][ke] );
                    }
                  }

                  // 3- Capillary: capillary pressure contribution
                  capillarity[ip][ix] = -1.0 * sign[ix] * sign[ip] * mobility[0] * mobility[1] / total_mobility * (capGrad[1] -capGrad[0]);
                  halfFluxVal[ip][ix] += capillarity[ip][ix];

                  for( integer ke = 0; ke < 2; ++ke )
                  {
                    real64 dC_dP = sign[ix] * sign[ip] * (dMob_dP[0][ke] * mobility[1] * mobility[1] + dMob_dP[1][ke] * mobility[0] * mobility[0]) / (total_mobility * total_mobility) *
                                   (capGrad[1] -capGrad[0]) +
                                   sign[ix] * sign[ip] * mobility[0] * mobility[1] / total_mobility * (dCapGrad_dP[1][ke] - dCapGrad_dP[0][ke]);

                    real64 dC_dS_term1 = sign[ix] * sign[ip] * (dMob_dS[0][ke] * mobility[1] * mobility[1] + dMob_dS[1][ke] * mobility[0] * mobility[0]) / (total_mobility * total_mobility) *
                                   (capGrad[1] -capGrad[0]);

                    real64 dC_dS_term2 = sign[ix] * sign[ip] * (mobility[0] * mobility[1]) / total_mobility;

                    real64 dC_dS_term3 = (dCapGrad_dS[1][ke] - dCapGrad_dS[0][ke]);

                    if( ix == 0 )
                    {
                      dC1_dS[ip][ke] = dC_dS_term1 + dC_dS_term2 * dC_dS_term3;
                      dhalfFlux1_dP[ip][ke] -= dC_dP;
                      dhalfFlux1_dS[ip][ke] -= dC1_dS[ip][ke];
                      if (std::fabs(facePhaseVolFrac1[0][0] - 1.0) > 1e-8){
                        dC1_dpc[ip][ke] =  dC1_dS[ip][ke] * dfacePhaseVolFrac_dCapPres1[0][0][0][0];
                      } else {
                        dC1_dpc[ip][ke] = dC_dS_term1 * dfacePhaseVolFrac_dCapPres1[0][0][0][0];
                      }
                      
                      // GEOS_UNUSED_VAR( dC1_dpc[ip][ke] );
                    }
                    else
                    {
                      dC2_dS[ip][ke] = dC_dS_term1 + dC_dS_term2 * dC_dS_term3;
                      dhalfFlux2_dP[ip][ke] -= dC_dP;
                      dhalfFlux2_dS[ip][ke] -= dC2_dS[ip][ke];
                      if (std::fabs(facePhaseVolFrac2[0][0] - 1.0) > 1e-8){
                        dC2_dpc[ip][ke] =  dC2_dS[ip][ke] * dfacePhaseVolFrac_dCapPres2[0][0][0][0];
                      } else {
                        dC2_dpc[ip][ke] = dC_dS_term1 * dfacePhaseVolFrac_dCapPres2[0][0][0][0];
                      }
                      // GEOS_UNUSED_VAR( dC2_dpc[ip][ke] );
                    }

                  }


                  if( halfFluxVal[ip][0] > 0.0 )
                  {
                    k_up_0_check[ip] = true;
                  }
                  if( halfFluxVal[ip][1] > 0.0 )
                  {
                    k_up_1_check[ip] = true;
                  }

                } // loop over phases

              } // loop over half fluxes

              check = true;

              // Brenier and Jaffre's PPU check:
              bool flip0[2] = {0, 0};
              bool flip0_k_up[2] = {0, 0};
              bool flip1[2] = {0, 0};
              bool flip1_k_up[2] = {0, 0};
              // loop over phases
              for( integer ip = 0; ip < 2; ++ip )
              {
                bool k_up_0_b = !static_cast< bool >(k_up_0[ip]);
                bool k_up_1_b = static_cast< bool >(k_up_1[ip]);
                flip0_k_up[ip] =  k_up_0_b;
                flip1_k_up[ip] =  k_up_1_b;

                if( std::fabs( uT ) < 1e-20 )
                {

                  if((std::fabs( halfFluxVal[ip][0] ) < 1e-20) && (std::fabs( halfFluxVal[ip][1] ) < 1e-20))
                  {
                    k_up_0_check[ip] = !k_up_0_b;
                    k_up_1_check[ip] = !k_up_1_b;
                  }
                  else
                  {
                    k_up_0_check[ip] = k_up_0_b;
                    k_up_1_check[ip] = k_up_1_b;
                  }

                  if( k_up_0_check[ip] != k_up_0_b )
                  {
                    flip0[ip] = 1;
                    check = false;
                  }

                  if( k_up_1_check[ip] != k_up_1_b )
                  {
                    flip1[ip] = 1;
                    check = false;
                  }

                }
                else
                {
                  if( std::fabs( halfFluxVal[ip][0] ) < 1e-20 )                                                 
                  {
                    k_up_0_check[ip] = k_up_0_b;
                  }

                  if( std::fabs( halfFluxVal[ip][1] ) < 1e-20 )
                  {
                    k_up_1_check[ip] = k_up_1_b;
                  }

                  if( k_up_0_check[ip] != k_up_0_b )
                  {
                    k_up_0[ip] = static_cast< localIndex >(k_up_0_b);
                    check = false;
                  }

                  if( k_up_1_check[ip] != k_up_1_b )
                  {
                    k_up_1[ip] = static_cast< localIndex >(!k_up_1_b);
                    check = false;
                  }
                }

              }

              if( flip0[0] || flip0[1] )
              {
                k_up_0[0] = static_cast< localIndex >(flip0_k_up[0]);
                k_up_0[1] = static_cast< localIndex >(flip0_k_up[1]);
              }

              if( flip1[0] || flip1[1] )
              {
                k_up_1[0] = static_cast< localIndex >(!flip1_k_up[0]);
                k_up_1[1] = static_cast< localIndex >(!flip1_k_up[1]);
              }

            } // while check for BJ PPU

            real64 constexpr eps2 = 1.0e-18;
            // newton update
            dhalfFlux_dpc[0][0] = dhalfFlux1_dS[0][1]*dfacePhaseVolFrac_dCapPres1[0][0][0][0];
            real64 dhalfFlux_dpc00 = dV1_dpc[0][1] - dG1_dpc[0][1] - dC1_dpc[0][1];
            dhalfFlux_dpc[0][1] = dhalfFlux2_dS[0][1]*dfacePhaseVolFrac_dCapPres2[0][0][0][0];
            real64 dhalfFlux_dpc01 = dV2_dpc[0][1] - dG2_dpc[0][1] - dC2_dpc[0][1];

            dhalfFlux_dpc[1][0] = dhalfFlux1_dS[1][1]*dfacePhaseVolFrac_dCapPres1[0][0][0][0];
            real64 dhalfFlux_dpc10 = dV1_dpc[1][1] - dG1_dpc[1][1] - dC1_dpc[1][1];
            dhalfFlux_dpc[1][1] = dhalfFlux2_dS[1][1]*dfacePhaseVolFrac_dCapPres2[0][0][0][0];
            real64 dhalfFlux_dpc11 = dV2_dpc[1][1] - dG2_dpc[1][1] - dC2_dpc[1][1];

            dhalfFlux_dpc[0][0] = dhalfFlux_dpc00;
            dhalfFlux_dpc[0][1] = dhalfFlux_dpc01;
            dhalfFlux_dpc[1][0] = dhalfFlux_dpc10;
            dhalfFlux_dpc[1][1] = dhalfFlux_dpc11;

            local_jacobian = dhalfFlux_dpc[0][0] - dhalfFlux_dpc[0][1] + eps2;
            local_residual = halfFluxVal[0][0] - halfFluxVal[0][1];

            // if (iter == 1) {
            //   grad_phi[0] = local_residual;
            // }

            // Check convergence
            if( std::fabs( local_residual ) < tol )
            {
              // if( std::fabs( dhalfFlux_dpc[0][0] - dhalfFlux_dpc[0][1] ) <eps2 ) {
              // std::cout << "**********************ZeroJacobian*******************" << std::endl;                
              // }
              converged = 1;
              break; // Converged
            }
            // converged = 0;
            if( (!converged) && ( iter > (max_iter - 2)) )
            {
              if( div == 0 )
              {
                iter = 0;
                div++;
                Pc_int = Pc_min_all;
              }
              else if( div == 1 )
              {
                iter = 0;
                div++;
                Pc_int = Pc_max_all;
              }
              else if( div > 1 )
              {
              local_jacobian = 0.0;
              std::cout << "**********************Diverged*******************" << std::endl;
              iter = max_iter;
              }
            }
            else
            {

              real64 deltaPc = local_residual/local_jacobian;

              if( std::fabs( local_residual ) < tol )
              {
              // if( std::fabs( dhalfFlux_dpc[0][0] - dhalfFlux_dpc[0][1] ) < eps2 ) {
              // std::cout << "**********************ZeroJacobian*******************" << std::endl;                
              // }
                converged = 1;  
                break; // Converged
              }

              // Damping option:
              if (damping) {
                real64 max_dpc = fmax( fabs(dCapPres1_dfacePhaseVolFrac[0][0][0][0]), fabs(dCapPres2_dfacePhaseVolFrac[0][0][0][0]));
  
                  real64 sign = std::copysign(1.0, deltaPc);
  
                  deltaPc = fmin( fabs(deltaPc), max_dpc * 0.2 );
                  deltaPc *= sign;
  
                }

                if (bisection && iter < 7){
                  if ( iter == 0 ){

                  old_Pc_int = Pc_int;
                  Pc_int = next_Pc_int;
                  old_residual = local_residual;

                  } else if (old_residual * local_residual < 0.0 ) { 
                  
                  Pc_int = (next_Pc_int + old_Pc_int) / 2.0;
                  old_residual = local_residual;
                  old_Pc_int = next_Pc_int;
                  next_Pc_int = Pc_int;

                  } else if (old_residual * local_residual > 0.0 ) {
                  
                  Pc_int = old_Pc_int;
                  old_residual = local_residual;
                  old_Pc_int = next_Pc_int;
                  next_Pc_int = Pc_int;

                  } else {

                  Pc_int = old_Pc_int;
                  old_residual = local_residual;
                  old_Pc_int = next_Pc_int;
                  next_Pc_int = Pc_int;

                }

                } else if (newton_path){
                  Pc_int -= next_Pc_int;

                } else {

                  Pc_int -= deltaPc;

                }
              

              // truncate the updated capillary pressure (extended capillary pressure condition) for reporting/plotting:

              real64 faceCapPres1_plot = fmin( Pc1_max, fmax( Pc_int, Pc1_min ));
              real64 faceCapPres2_plot = fmin( Pc2_max, fmax( Pc_int, Pc2_min ));
              faceCapPres1_plot = fmin( Pc2_max, fmax( faceCapPres1_plot, Pc2_min ));
              faceCapPres2_plot = fmin( Pc1_max, fmax( faceCapPres2_plot, Pc1_min ));
              
              // Write data to the file
              outFile << GEOS_FMT( "{:10.10e}", local_jacobian );
              outFile << GEOS_FMT( ",{:10.10e}", local_residual );
              outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[0][0] );
              outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[0][1] );
              outFile << GEOS_FMT( ",{:10.10e}", Pc_int_iterate );
              outFile << GEOS_FMT( ",{:10.10e}", faceCapPres1[0][0][0] );
              outFile << GEOS_FMT( ",{:10.10e}", faceCapPres2[0][0][0] );
              outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[1][0] );
              outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[1][1] );
              outFile << GEOS_FMT( ",{:10.10e}", viscous[0][0] );
              outFile << GEOS_FMT( ",{:10.10e}", viscous[1][0] );
              outFile << GEOS_FMT( ",{:10.10e}", viscous[0][1] );
              outFile << GEOS_FMT( ",{:10.10e}", viscous[1][1] );
              outFile << GEOS_FMT( ",{:10.10e}", bouyancy[0][0] );
              outFile << GEOS_FMT( ",{:10.10e}", bouyancy[1][0] );
              outFile << GEOS_FMT( ",{:10.10e}", bouyancy[0][1] );
              outFile << GEOS_FMT( ",{:10.10e}", bouyancy[1][1] );
              outFile << GEOS_FMT( ",{:10.10e}", capillarity[0][0] );
              outFile << GEOS_FMT( ",{:10.10e}", capillarity[1][0] );
              outFile << GEOS_FMT( ",{:10.10e}", capillarity[0][1] );
              outFile << GEOS_FMT( ",{:10.10e}", capillarity[1][1] );
              outFile << GEOS_FMT( ",{:10.10e}", capPres1[0][0][0] );
              outFile << GEOS_FMT( ",{:10.10e}", capPres1[0][0][1] );
              outFile << GEOS_FMT( ",{:10.10e}", facePhaseVolFrac1[0][0] );
              outFile << GEOS_FMT( ",{:10.10e}", facePhaseVolFrac2[0][0] );
              outFile << std::endl;

              iter++;

            }

            
          } // while loop

          if( converged )
          {

            // Global derivatives:
            real64 constexpr eps3 = 1.0e-18;

            real64 const dPc_int_dS1 =(-1.0) * (dhalfFlux1_dS[0][0] +  dhalfFlux_duT[0][0] * duT_dS[0] - dhalfFlux_duT[0][1] * duT_dS[0]) / (dhalfFlux_dpc[0][0] - dhalfFlux_dpc[0][1] + eps3);
            real64 const dPc_int_dS2 =(-1.0) * (dhalfFlux_duT[0][0] * duT_dS[1] - dhalfFlux2_dS[0][0] - dhalfFlux_duT[0][1] * duT_dS[1]) / (dhalfFlux_dpc[0][0] - dhalfFlux_dpc[0][1] + eps3);
            real64 const dPc_int_du =(-1.0) * (dhalfFlux_duT[0][0]  - dhalfFlux_duT[0][1]) / (dhalfFlux_dpc[0][0] - dhalfFlux_dpc[0][1] + eps3);

            dFlux_dP[0][0] =  (dhalfFlux_duT[0][0] * duT_dP[0] + dhalfFlux_dpc[0][0] * dPc_int_du * duT_dP[0]) * density2[0] + halfFluxVal[0][0] * dDens_dP2[0][0];
            dFlux_dS[0][0] = (dhalfFlux1_dS[0][0] +  dhalfFlux_duT[0][0] * duT_dS[0] + dhalfFlux_dpc[0][0] * dPc_int_dS1) * density2[0];

            dFlux_dP[0][1] = (dhalfFlux_duT[0][1] * duT_dP[1] + dhalfFlux_dpc[0][1] * dPc_int_du * duT_dP[1]) * density2[0] + halfFluxVal[0][1] * dDens_dP2[0][1];
            dFlux_dS[0][1] = (dhalfFlux2_dS[0][0] +  dhalfFlux_duT[0][1] * duT_dS[1] + dhalfFlux_dpc[0][1] * dPc_int_dS2) * density2[0];

            dFlux_dP[1][0] =  (dhalfFlux_duT[1][0] * duT_dP[0] + dhalfFlux_dpc[1][0] * dPc_int_du * duT_dP[0]) * density2[1] + halfFluxVal[1][0] * dDens_dP2[1][0];
            dFlux_dS[1][0] = (dhalfFlux1_dS[1][0] +  dhalfFlux_duT[1][0] * duT_dS[0] + dhalfFlux_dpc[1][0] * dPc_int_dS1) * density2[1];

            dFlux_dP[1][1] = (dhalfFlux_duT[1][1] * duT_dP[1] + dhalfFlux_dpc[1][1] * dPc_int_du * duT_dP[1]) * density2[1] + halfFluxVal[1][1] * dDens_dP2[1][1];
            dFlux_dS[1][1] = (dhalfFlux2_dS[1][0] +  dhalfFlux_duT[1][1] * duT_dS[1] + dhalfFlux_dpc[1][1] * dPc_int_dS2) * density2[1];

            fluxVal[0] = halfFluxVal[0][0] * density2[0];
            fluxVal[1] = halfFluxVal[1][0] * density2[1];

          } else {
            std::cout << "**********************Diverged*******************" << std::endl;

          }
          
          phi[0] = fluxVal[0];
          phi[1] = fluxVal[1];

          grad_phi_P[0] = dFlux_dP[0][0];
          grad_phi_P[1] = dFlux_dP[0][1];
          grad_phi_P[2] = dFlux_dP[1][0];
          grad_phi_P[3] = dFlux_dP[1][1];

          grad_phi_S[0] = dFlux_dS[0][0];
          grad_phi_S[1] = dFlux_dS[0][1];
          grad_phi_S[2] = dFlux_dS[1][0];
          grad_phi_S[3] = dFlux_dS[1][1];

          converged = 1;
          GEOS_UNUSED_VAR( converged );
// Close the file after writing
          outFile.close();

        } );

      } );

    } );

  } );


}

/******************************** FluxComputeKernelBase ********************************/

/**
 * @brief Base class for FluxComputeKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of dofs).
 */
class FluxComputeKernelBase
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using ImmiscibleMultiphaseFlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::flow::pressure,
                      fields::flow::gravityCoefficient,
                      fields::immiscibleMultiphaseFlow::phaseVolumeFraction,
                      fields::immiscibleMultiphaseFlow::phaseMobility,
                      fields::immiscibleMultiphaseFlow::phaseMass_n,
                      fields::immiscibleMultiphaseFlow::dPhaseMobility >;

  using MultiphaseFluidAccessors =
    StencilMaterialAccessors< constitutive::TwoPhaseImmiscibleFluid,
                              fields::twophaseimmisciblefluid::phaseDensity,
                              fields::twophaseimmisciblefluid::dPhaseDensity,
                              fields::twophaseimmisciblefluid::phaseViscosity,
                              fields::twophaseimmisciblefluid::dPhaseViscosity >;

  using CapPressureAccessors =
    StencilMaterialAccessors< BrooksCoreyCapillaryPressure,
                              fields::cappres::phaseCapPressure,
                              fields::cappres::dPhaseCapPressure_dPhaseVolFraction
                              >;
                              // ,fields::cappres::jFuncMultiplier >;


  using PermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;

  using Deriv = immiscibleFlow::DerivativeOffset;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumberAccessor accessor for the dof numbers
   * @param[in] multiPhaseFlowAccessors accessor for wrappers registered by the solver
   * @param[in] fluidAccessors accessor for wrappers registered by the fluid model
   * @param[in] capPressureAccessors accessor for wrappers registered by the capillary pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[inout] hasCapPressure flag to indicate whether problem includes capillarity
   * @param[inout] useTotalMassEquation flag to indicate whether to use the total mass formulation
   */
  FluxComputeKernelBase( integer const numPhases,
                         globalIndex const rankOffset,
                         DofNumberAccessor const & dofNumberAccessor,
                         ImmiscibleMultiphaseFlowAccessors const & multiPhaseFlowAccessors,
                         MultiphaseFluidAccessors const & fluidAccessors,
                         CapPressureAccessors const & capPressureAccessors,
                         PermeabilityAccessors const & permeabilityAccessors,
                         real64 const & dt,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         integer const hasCapPressure,
                         integer const useTotalMassEquation,
                         integer const checkPhasePresenceInGravity )
    : m_numPhases ( numPhases ),
    m_rankOffset( rankOffset ),
    m_dt( dt ),
    m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_permeability( permeabilityAccessors.get( fields::permeability::permeability {} ) ),
    m_dPerm_dPres( permeabilityAccessors.get( fields::permeability::dPerm_dPressure {} ) ),
    m_ghostRank( multiPhaseFlowAccessors.get( fields::ghostRank {} ) ),
    m_gravCoef( multiPhaseFlowAccessors.get( fields::flow::gravityCoefficient {} ) ),
    m_pres( multiPhaseFlowAccessors.get( fields::flow::pressure {} ) ),
    m_phaseVolFrac( multiPhaseFlowAccessors.get( fields::immiscibleMultiphaseFlow::phaseVolumeFraction {} ) ),
    m_phaseMass_n( multiPhaseFlowAccessors.get( fields::immiscibleMultiphaseFlow::phaseMass_n {} ) ),
    m_mob( multiPhaseFlowAccessors.get( fields::immiscibleMultiphaseFlow::phaseMobility {} ) ),
    m_dMob( multiPhaseFlowAccessors.get( fields::immiscibleMultiphaseFlow::dPhaseMobility {} ) ),
    m_dens( fluidAccessors.get( fields::twophaseimmisciblefluid::phaseDensity {} ) ),
    m_visc( fluidAccessors.get( fields::twophaseimmisciblefluid::phaseViscosity {} ) ),
    m_dDens_dPres( fluidAccessors.get( fields::twophaseimmisciblefluid::dPhaseDensity {} ) ),
    m_dVisc_dPres( fluidAccessors.get( fields::twophaseimmisciblefluid::dPhaseViscosity {} ) ),
    m_phaseCapPressure( capPressureAccessors.get( fields::cappres::phaseCapPressure {} ) ),
    m_dPhaseCapPressure_dPhaseVolFrac( capPressureAccessors.get( fields::cappres::dPhaseCapPressure_dPhaseVolFraction {} ) ),
    // m_jFuncMultiplier( capPressureAccessors.get( fields::cappres::jFuncMultiplier {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs ),
    m_hasCapPressure ( hasCapPressure ),
    m_useTotalMassEquation ( useTotalMassEquation ),
    m_checkPhasePresenceInGravity ( checkPhasePresenceInGravity )
  {}

protected:

  /// Number of fluid phases
  integer const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Time step size
  real64 const m_dt;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > m_permeability;
  ElementViewConst< arrayView3d< real64 const > > m_dPerm_dPres;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;

  // Primary and secondary variables
  /// Views on pressure and phase volume fraction
  ElementViewConst< arrayView1d< real64 const > > const m_pres;
  ElementViewConst< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > > const m_phaseVolFrac;
  ElementViewConst< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > > const m_phaseMass_n;

  /// Views on fluid mobility
  ElementViewConst< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > > const m_mob;
  ElementViewConst< arrayView3d< real64 const, immiscibleFlow::USD_PHASE_DS > > const m_dMob;

  /// Views on fluid density and viscosity
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_dens;
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_visc;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dDens_dPres;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dVisc_dPres;

  /// Views on capillary pressure
  ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const m_phaseCapPressure;
  ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const m_dPhaseCapPressure_dPhaseVolFrac;
  // ElementViewConst< arrayView2d< real64 const > > const m_jFuncMultiplier;

  // Residual and jacobian


  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  // Flags
  integer const m_hasCapPressure;
  integer const m_useTotalMassEquation;
  integer const m_checkPhasePresenceInGravity;
};

/***************************************** */

/**
 * @class FluxComputeKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER >
class FluxComputeKernel : public FluxComputeKernelBase
{
public:

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_EQN;

  /// Maximum number of elements at the face
  static constexpr localIndex maxNumElems = STENCILWRAPPER::maxNumPointsInFlux;

  /// Maximum number of connections at the face
  static constexpr localIndex maxNumConns = STENCILWRAPPER::maxNumConnections;

  /// Maximum number of points in the stencil
  static constexpr localIndex maxStencilSize = STENCILWRAPPER::maxStencilSize;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] multiPhaseFlowAccessors
   * @param[in] fluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] hasCapPressure flags for capillary pressure
   * @param[in] useTotalMassEquation flags for using total velocity formulation
   */
  FluxComputeKernel( integer const numPhases,
                     globalIndex const rankOffset,
                     STENCILWRAPPER const & stencilWrapper,
                     DofNumberAccessor const & dofNumberAccessor,
                     ImmiscibleMultiphaseFlowAccessors const & multiPhaseFlowAccessors,
                     MultiphaseFluidAccessors const & fluidAccessors,
                     CapPressureAccessors const & capPressureAccessors,
                     PermeabilityAccessors const & permeabilityAccessors,
                     real64 const & dt,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs,
                     integer const hasCapPressure,
                     integer const useTotalMassEquation,
                     integer const checkPhasePresenceInGravity )
    : FluxComputeKernelBase( numPhases,
                             rankOffset,
                             dofNumberAccessor,
                             multiPhaseFlowAccessors,
                             fluidAccessors,
                             capPressureAccessors,
                             permeabilityAccessors,
                             dt,
                             localMatrix,
                             localRhs,
                             hasCapPressure,
                             useTotalMassEquation,
                             checkPhasePresenceInGravity ),
    m_stencilWrapper( stencilWrapper ),
    m_seri( stencilWrapper.getElementRegionIndices() ),
    m_sesri( stencilWrapper.getElementSubRegionIndices() ),
    m_sei( stencilWrapper.getElementIndices() )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : stencilSize( size ),
      numFluxElems( numElems ),
      dofColIndices( size * numDof ),
      localFlux( numElems * numEqn ),
      localFluxJacobian( numElems * numEqn, size * numDof )
    {}

    // Stencil information

    /// Stencil size for a given connection
    localIndex const stencilSize;

    /// Number of elements for a given connection
    localIndex const numFluxElems;

    // Transmissibility and derivatives

    /// Transmissibility
    real64 transmissibility[maxNumConns][2]{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTrans_dPres[maxNumConns][2]{};

    /// Transmissibility
    real64 transmissibilityHat[maxNumConns][2]{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTransHat_dPres[maxNumConns][2]{};

    // Local degrees of freedom and local residual/jacobian

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    stackArray1d< globalIndex, maxNumElems * numDof > dofColIndices;

    /// Storage for the face local residual vector (all equations except volume balance)
    stackArray1d< real64, maxNumElems * numEqn > localFlux;
    /// Storage for the face local Jacobian matrix
    stackArray2d< real64, maxNumElems * numEqn * maxStencilSize * numDof > localFluxJacobian;
  };

  /**
   * @brief Getter for the stencil size at this connection
   * @param[in] iconn the connection index
   * @return the size of the stencil at this connection
   */
  GEOS_HOST_DEVICE
  localIndex stencilSize( localIndex const iconn ) const
  { return m_sei[iconn].size(); }

  /**
   * @brief Getter for the number of elements at this connection
   * @param[in] iconn the connection index
   * @return the number of elements at this connection
   */
  GEOS_HOST_DEVICE
  localIndex numPointsInFlux( localIndex const iconn ) const
  { return m_stencilWrapper.numPointsInFlux( iconn ); }

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */

  GEOS_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    // set degrees of freedom indices for this face
    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      globalIndex const offset = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];

      for( integer jdof = 0; jdof < numDof; ++jdof )
      {
        stack.dofColIndices[i * numDof + jdof] = offset + jdof;
      }
    }
  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] NoOpFunc the function used to customize the computation of the flux
   */
  template< typename FUNC = NoOpFunc >    // should change to multiphase
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && kernelOp = NoOpFunc{} ) const
  {
    // first, compute the transmissibilities at this face                                             // get k and dk/dP from global arrays
    // and place in stack
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );

    localIndex k[2];
    localIndex connectionIndex = 0;

    for( k[0] = 0; k[0] < stack.numFluxElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numFluxElems; ++k[1] )
      {
        // clear working arrays
        real64 densMean[numEqn]{};
        real64 dDensMean_dP[numEqn][2]{};

        real64 presGrad[numEqn]{};
        real64 dPresGrad_dP[numEqn][2]{};

        real64 gravHead[numEqn]{};
        real64 dGravHead_dP[numEqn][2]{};

        real64 capGrad[numEqn]{};
        real64 dCapGrad_dP[numEqn][2]{};
        real64 dCapGrad_dS[numEqn][2]{};

        real64 fluxVal[numEqn]{};
        real64 dFlux_dP[numEqn][2]{};
        real64 dFlux_dS[numEqn][2]{};

        real64 mobility[numEqn]{};
        real64 dMob_dP[numEqn][2]{};
        real64 dMob_dS[numEqn][2]{};

        real64 const trans[2] = { stack.transmissibility[connectionIndex][0], stack.transmissibility[connectionIndex][1] };
        real64 const dTrans_dP[2] = { stack.dTrans_dPres[connectionIndex][0], stack.dTrans_dPres[connectionIndex][1] };

        // cell indices
        localIndex const seri[2]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        // loop over phases
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {
          // calculate quantities on primary connected cells
          integer denom = 0;
          for( integer ke = 0; ke < 2; ++ke )
          {
            // density
            bool const phaseExists = (m_phaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][ip] > 0);
            if( m_checkPhasePresenceInGravity && !phaseExists )
            {
              continue;
            }

            real64 const density  = m_dens[seri[ke]][sesri[ke]][sei[ke]][0][ip];                     // r = rho1 || rho2
            real64 const dDens_dP = m_dDens_dPres[seri[ke]][sesri[ke]][sei[ke]][0][ip][Deriv::dP];   // dr/dP = dr1/dP1 || dr2/dP

            // average density and derivatives
            densMean[ip] += density;           // rho = (rho1 + rho2)
            dDensMean_dP[ip][ke] = dDens_dP;   // drho/dP = { (dr1/dP1) , (dr2/dP2) }

            denom++;
          }

          if( denom > 1 )
          {
            densMean[ip] /= denom;  // rho = (rho1 + rho2) / denom
            for( integer ke = 0; ke < 2; ++ke )
            {
              dDensMean_dP[ip][ke] /= denom;  // drho/dP = { (dr1/dP1) / denom , (dr2/dP2) / denom }
            }
          }

          //***** calculation of flux *****

          // compute potential difference
          real64 potScale = 0.0;
          real64 dPresGrad_dTrans = 0.0;
          real64 dGravHead_dTrans = 0.0;
          real64 dCapGrad_dTrans = 0.0;
          constexpr int signPotDiff[2] = {1, -1};

          for( integer ke = 0; ke < 2; ++ke )
          {
            localIndex const er  = seri[ke];
            localIndex const esr = sesri[ke];
            localIndex const ei  = sei[ke];

            real64 const pressure = m_pres[er][esr][ei];       // P = P1 || P2
            presGrad[ip] += trans[ke] * pressure;              // DPv = T (P1 - P2)
            dPresGrad_dTrans += signPotDiff[ke] * pressure;    // dDPv/dT = (P1 - P2)
            dPresGrad_dP[ip][ke] = trans[ke];                  // dDPv/dP = { T , -T }

            real64 const gravD = trans[ke] * m_gravCoef[er][esr][ei];        // D = T g z1 || -T g z2
            real64 pot = trans[ke] * pressure - densMean[ip] * gravD;  // Phi = T P1 - rho T g z1 || -T P2 + rho T g z2

            gravHead[ip] += densMean[ip] * gravD;                                          // DPg = rho (T g z1 - T g z2) = T rho g (z1 -
                                                                                           // z2)
            dGravHead_dTrans += signPotDiff[ke] * densMean[ip] * m_gravCoef[er][esr][ei];  // dDPg/dT = rho g z1 - rho g z2 = rho g (z1 -
                                                                                           // z2)

            for( integer i = 0; i < 2; ++i )
            {
              dGravHead_dP[ip][i] += dDensMean_dP[ip][i] * gravD;  // dDPg/dP = {drho/dP1 * T g (z1 - z2) , drho/dP2 * T g (z1 - z2)}
            }

            if( m_hasCapPressure )   // check sign convention
            {
              real64 const capPres = m_phaseCapPressure[er][esr][ei][0][ip];   // Pc = Pc1 || Pc2
              dCapGrad_dTrans -= signPotDiff[ke] * capPres;                    // dDPc/dT = (-Pc1 + Pc2)
              pot -= trans[ke] * capPres;                                      // Phi = T P1 - rho T g z1 - T Pc1 || -T P2 + rho T g z2 + T
                                                                               // Pc2
              capGrad[ip] -= trans[ke] * capPres;                              // DPc = T (-Pc1 + Pc2)
            }

            potScale = fmax( potScale, fabs( pot ) );  // maxPhi = Phi1 > Phi2 ? Phi1 : Phi2
          }

          for( integer ke = 0; ke < 2; ++ke )
          {
            dPresGrad_dP[ip][ke] += dTrans_dP[ke] * dPresGrad_dTrans;    // dDPv/dP = { T + dT/dP1 * (P1 - P2) , -T + dT/dP2 * (P1 - P2)}
            dGravHead_dP[ip][ke] += dTrans_dP[ke] * dGravHead_dTrans;    // dDPg/dP = { drho/dP1 * T g (z1 - z2) + dT/dP1 * rho g (z1 - z2)
                                                                         // ,
                                                                         //             drho/dP2 * T g (z1 - z2) + dT/dP2 * rho g (z1 - z2)
                                                                         // }
            if( m_hasCapPressure )
            {
              real64 const dCapPres_dS = m_dPhaseCapPressure_dPhaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][0][ip][ip];  // dPc/dS = dPc1/dS1
                                                                                                                      // ||
                                                                                                                      // dPc2/dS2
              dCapGrad_dP[ip][ke] += dTrans_dP[ke] * dCapGrad_dTrans;                                                 // dDPc/dP = { dT/dP1
                                                                                                                      // *
                                                                                                                      // (-Pc1 + Pc2) ,
                                                                                                                      //             dT/dP2
                                                                                                                      // *
                                                                                                                      // (-Pc1 + Pc2) }
              dCapGrad_dS[ip][ke] -= trans[ke] * dCapPres_dS;                                                         // dDPc/dS = { -T *
                                                                                                                      // dPc1/dS1 , T *
                                                                                                                      // dPc2/dS2 }
            }
          }

          // *** upwinding ***

          // compute potential gradient
          real64 potGrad = presGrad[ip] - gravHead[ip];  // DPhi = T (P1 - P2) - T rho g (z1 - z2)
          if( m_hasCapPressure )
          {
            potGrad += capGrad[ip];  // DPhi = T (P1 - P2) - T rho g (z1 - z2) + T (-Pc1 + Pc2)
          }

          // compute upwinding tolerance
          real64 constexpr upwRelTol = 1e-8;
          real64 const upwAbsTol = fmax( potScale * upwRelTol, LvArray::NumericLimits< real64 >::epsilon );  // abstol = maxPhi * tol > eps
                                                                                                             // ?
                                                                                                             // maxPhi * tol : eps

          // decide mobility coefficients - smooth variation in [-upwAbsTol; upwAbsTol]
          real64 const alpha = ( potGrad + upwAbsTol ) / ( 2 * upwAbsTol );    // alpha = (DPhi + abstol) / abstol / 2

          // choose upstream cell
          if( alpha <= 0.0 || alpha >= 1.0 )   // no smoothing needed
          {
            localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );  // 1 upwind -> k_up = 0 || 2 upwind -> k_up = 1

            mobility[ip] = m_mob[seri[k_up]][sesri[k_up]][sei[k_up]][ip];                      // M = Mupstream
            dMob_dP[ip][k_up] = m_dMob[seri[k_up]][sesri[k_up]][sei[k_up]][ip][Deriv::dP];     // dM/dP = {dM/dP1 , 0} OR {0 , dM/dP2}
            dMob_dS[ip][k_up] = m_dMob[seri[k_up]][sesri[k_up]][sei[k_up]][ip][Deriv::dS];     // dM/dS = {dM/dS1 , 0} OR {0 , dM/dS2}
          }
          else   // perform smoothing
          {
            real64 const mobWeights[2] = { alpha, 1.0 - alpha };
            for( integer ke = 0; ke < 2; ++ke )
            {
              mobility[ip] += mobWeights[ke] * m_mob[seri[ke]][sesri[ke]][sei[ke]][ip];                // M = alpha * M1 + (1 - alpha) * M2
              dMob_dP[ip][ke] = mobWeights[ke] * m_dMob[seri[ke]][sesri[ke]][sei[ke]][ip][Deriv::dP];  // dM/dP = {alpha * dM1/dP1 , (1 -
                                                                                                       // alpha) * dM2/dP2}
              dMob_dS[ip][ke] = mobWeights[ke] * m_dMob[seri[ke]][sesri[ke]][sei[ke]][ip][Deriv::dS];  // dM/dP = {alpha * dM1/dS1 , (1 -
                                                                                                       // alpha) * dM2/dS2}
            }
          }

          // pressure gradient depends on all points in the stencil
          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] += dPresGrad_dP[ip][ke];  // dF/dP = { T + dT/dP1 * (P1 - P2) ,
                                                       //          -T + dT/dP2 * (P1 - P2) }
          }

          // gravitational head depends only on the two cells connected (same as mean density)
          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] -= dGravHead_dP[ip][ke];  // dF/dP = { T + dT/dP1 * (P1 - P2) - drho/dP1 * T g (z1 - z2) - dT/dP1 * rho g (z1 -
                                                       // z2) ,
                                                       //          -T + dT/dP2 * (P1 - P2) - drho/dP2 * T g (z1 - z2) - dT/dP2 * rho g (z1 -
                                                       // z2) }
          }

          // capillary pressure contribution
          if( m_hasCapPressure )
          {
            for( integer ke = 0; ke < 2; ++ke )
            {
              dFlux_dP[ip][ke] += dCapGrad_dP[ip][ke];   // dF/dP = { T + dT/dP1 * (P1 - P2) - drho/dP1 * T g (z1 - z2) - dT/dP1 * rho g (z1
                                                         // - z2) + dT/dP1 * (-Pc1 + Pc2) ,
                                                         //          -T + dT/dP2 * (P1 - P2) - drho/dP2 * T g (z1 - z2) - dT/dP2 * rho g (z1
                                                         // - z2) + dT/dP2 * (-Pc1 + Pc2) }

              dFlux_dS[ip][ke] += dCapGrad_dS[ip][ke];   // dF/dS = { T * -dPc/dS1 , T * dPc/dS2 }
            }
          }

          // compute the flux and derivatives using upstream cell mobility
          fluxVal[ip] = mobility[ip] * potGrad;  // F = M * DPhi

          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] *= mobility[ip];    // dF/dP = { M [ T + dT/dP1 * (P1 - P2) - drho/dP1 * T g (z1 - z2) - dT/dP1 * rho g (z1 -
                                                 // z2) + dT/dP1 * (-Pc1 + Pc2)] ,
                                                 //           M [-T + dT/dP2 * (P1 - P2) - drho/dP2 * T g (z1 - z2) - dT/dP2 * rho g (z1 -
                                                 // z2) + dT/dP2 * (-Pc1 + Pc2)] }

            dFlux_dS[ip][ke] *= mobility[ip];    // dF/dS = { M [T * -dPc/dS1] , M [T * dPc/dS2] }
          }

          // add contribution from upstream cell mobility derivatives
          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] += dMob_dP[ip][ke] * potGrad;   // dF/dP = { M [ T + dT/dP1 * (P1 - P2) - drho1/dP * T g (z1 - z2) - dT/dP1 *
                                                             // rho g (z1 - z2) + dT/dP1 * (-Pc1 + Pc2)] + dM/dP1 * DPhi ,
                                                             //           M [-T + dT/dP2 * (P1 - P2) - drho2/dP * T g (z1 - z2) - dT/dP2 *
                                                             // rho g (z1 - z2) + dT/dP2 * (-Pc1 + Pc2)] + dM/dP2 * DPhi }

            dFlux_dS[ip][ke] += dMob_dS[ip][ke] * potGrad;   // dF/dS = { M [T * -dPc/dS1] + dM/dS1 * DPhi , M [T * dPc/dS2] + dM/dS2 * DPhi
                                                             // }
          }

          // populate local flux vector and derivatives
          stack.localFlux[k[0]*numEqn + ip] += m_dt * fluxVal[ip];
          stack.localFlux[k[1]*numEqn + ip] -= m_dt * fluxVal[ip];

          for( integer ke = 0; ke < 2; ++ke )
          {
            // pressure
            localIndex const localDofIndexPres = k[ke] * numDof;
            stack.localFluxJacobian[k[0]*numEqn + ip][localDofIndexPres] += m_dt * dFlux_dP[ip][ke];
            stack.localFluxJacobian[k[1]*numEqn + ip][localDofIndexPres] -= m_dt * dFlux_dP[ip][ke];

            // saturation (hard-coded for 2-phase currently)
            localIndex const localDofIndexSat = k[ke] * numDof + 1;
            stack.localFluxJacobian[k[0]*numEqn + ip][localDofIndexSat] += m_dt * dFlux_dS[ip][ke];
            stack.localFluxJacobian[k[1]*numEqn + ip][localDofIndexSat] -= m_dt * dFlux_dS[ip][ke];
          }

          // Customize the kernel with this lambda
          kernelOp( k, seri, sesri, sei, connectionIndex, alpha, mobility, potGrad, fluxVal, dFlux_dP, dFlux_dS );   // Not sure what this
                                                                                                                     // does

        }  // loop over phases

        connectionIndex++;
      }
    }  // loop over connection elements
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = NoOpFunc >                                   // should change to multiphase
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && kernelOp = NoOpFunc{} ) const
  {
    using namespace compositionalMultiphaseUtilities;

    if( m_useTotalMassEquation )
    {
      // Apply equation/variable change transformation(s)
      stackArray1d< real64, maxStencilSize * numDof > work( stack.stencilSize * numDof );
      shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numPhases, numEqn, numDof * stack.stencilSize, stack.numFluxElems,
                                                               stack.localFluxJacobian, work );
      shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( m_numPhases, numEqn, stack.numFluxElems,
                                                                 stack.localFlux );
    }

    // add contribution to residual and jacobian into:
    // - the mass balance equation
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    for( integer i = 0; i < stack.numFluxElems; ++i )
    {
      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOS_ASSERT_GE( localRow, 0 );

        GEOS_ASSERT_GT( m_localMatrix.numRows(), localRow + numEqn - 1 );

        for( integer ic = 0; ic < numEqn; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow + ic],
                           stack.localFlux[i * numEqn + ic] );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
            ( localRow + ic,
            stack.dofColIndices.data(),
            stack.localFluxJacobian[i * numEqn + ic].dataIfContiguous(),
            stack.stencilSize * numDof );
        }

        // call the lambda to assemble additional terms, such as thermal terms
        kernelOp( i, localRow );
      }
    }
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numConnections the number of connections
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numConnections,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numConnections, [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      typename KERNEL_TYPE::StackVariables stack( kernelComponent.stencilSize( iconn ),
                                                  kernelComponent.numPointsInFlux( iconn ) );

      kernelComponent.setup( iconn, stack );
      kernelComponent.computeFlux( iconn, stack );
      kernelComponent.complete( iconn, stack );
    } );
  }

protected:

  // Stencil information

  /// Reference to the stencil wrapper
  STENCILWRAPPER const m_stencilWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;
};


/**
 * @class FaceBasedAssemblyInterfaceConditionKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @tparam CAPPRESWRAPPER the type of the capillary pressure wrapper
 * @tparam RELPERMWRAPPER the type of the realtive permeability wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
// template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER, typename CAPPRESWRAPPER, typename RELPERMWRAPPER >
template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER >
class FluxComputeInterfaceConditionKernel : public FluxComputeKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >
{
public:

  using AbstractBase = FluxComputeKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using ImmiscibleMultiphaseFlowAccessors = AbstractBase::ImmiscibleMultiphaseFlowAccessors;
  using MultiphaseFluidAccessors = AbstractBase::MultiphaseFluidAccessors;
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using Base = FluxComputeKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
  using Deriv = typename Base::Deriv;
  using StackVariables = typename Base::StackVariables;
  using Base::numEqn;
  using Base::numDof;
  using Base::m_stencilWrapper;
  using Base::m_dofNumber;
  using Base::m_rankOffset;
  using Base::m_localMatrix;
  using Base::m_localRhs;
  using Base::m_numPhases;
  using Base::m_permeability;
  using Base::m_dPerm_dPres;
  using Base::m_phaseVolFrac;
  using Base::m_phaseMass_n;
  using Base::m_phaseCapPressure;
  // using Base::m_jFuncMultiplier;
  using Base::m_dPhaseCapPressure_dPhaseVolFrac;
  using Base::m_dens;
  using Base::m_dDens_dPres;
  using Base::m_visc;
  using Base::m_dVisc_dPres;
  using Base::m_dMob;
  using Base::m_mob;
  using Base::m_gravCoef;
  using Base::m_ghostRank;
  using Base::m_dt;
  using Base::m_hasCapPressure;
  using Base::m_useTotalMassEquation;
  using Base::m_checkPhasePresenceInGravity;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;
  using Base::m_pres;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] capPressureWrapper reference to the capillary pressure wrapper
   * @param[in] dofNumberAccessor
   * @param[in] multiPhaseFlowAccessors
   * @param[in] fluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] hasCapPressure flags for capillary pressure
   * @param[in] useTotalMassEquation flags for using total velocity formulation
   */
  FluxComputeInterfaceConditionKernel( integer const numPhases,
                                       globalIndex const rankOffset,
                                       STENCILWRAPPER const & stencilWrapper,
                                      //  CAPPRESWRAPPER const & capPressureWrapper,
                                      //  RELPERMWRAPPER const & relPermWrapper,
                                       DofNumberAccessor const & dofNumberAccessor,
                                       ImmiscibleMultiphaseFlowAccessors const & multiPhaseFlowAccessors,
                                       MultiphaseFluidAccessors const & fluidAccessors,
                                       CapPressureAccessors const & capPressureAccessors,
                                       PermeabilityAccessors const & permeabilityAccessors,
                                       real64 const & dt,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs,
                                       integer const hasCapPressure,
                                       integer const useTotalMassEquation,
                                       integer const checkPhasePresenceInGravity,
                                       string_array const & interfaceFaceSetNames,
                                       stdVector< std::array< std::tuple< constitutive::RelativePermeabilityBase *,
                                                                          constitutive::CapillaryPressureBase *,
                                                                          constitutive::TwoPhaseImmiscibleFluid * >, 2 > > const & interfaceConstitutivePairs,
                                       unordered_map< localIndex, localIndex > const & interfaceRegionByConnector,
                                      //  std::tuple< constitutive::RelativePermeabilityBase *,
                                      //              constitutive::CapillaryPressureBase *,
                                      //              constitutive::TwoPhaseImmiscibleFluid * > const & interfaceConstitutivePairs_temp,
                                       localIndex const GEOS_UNUSED_PARAM( domainSize ) )
    : Base( numPhases,
            rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            multiPhaseFlowAccessors,
            fluidAccessors,
            capPressureAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs,
            hasCapPressure,
            useTotalMassEquation,
            checkPhasePresenceInGravity ),
    // m_capPressureWrapper( capPressureWrapper ),
    // m_relPermWrapper( relPermWrapper ),
    m_interfaceFaceSetNames( interfaceFaceSetNames ),
    m_interfaceConstitutivePairs( interfaceConstitutivePairs ),
    m_interfaceRegionByConnector( interfaceRegionByConnector )
    // ,
    // m_interfaceConstitutivePairs_temp( interfaceConstitutivePairs_temp )
  {}

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] NoOpFunc the function used to customize the computation of the flux
   */

  template< typename FUNC = NoOpFunc >    // should change to multiphase
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && kernelOp = NoOpFunc{} ) const
  {

     bool connectorHasInterfaceConditionQ = false;
     bool anyInterfaceConditionsQ = not m_interfaceConstitutivePairs.empty();
     if (anyInterfaceConditionsQ) {
         connectorHasInterfaceConditionQ =
             m_interfaceRegionByConnector.find(iconn) != m_interfaceRegionByConnector.end();
     }



    //  if (connectorHasInterfaceConditionQ){
    //    // Improved transmission conditions
    //    int ammar_code = 0;
    //  }else{
    //    // Regular contribution
    //    int standard_code = 0;
    //  }

    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );

    localIndex k[2];
    localIndex connectionIndex = 0;

    // one-sided transmissibility
    m_stencilWrapper.computeHalfWeights( iconn,
                                         m_permeability,
                                         m_dPerm_dPres,
                                         stack.transmissibilityHat,
                                         stack.dTransHat_dPres );



    for( k[0] = 0; k[0] < stack.numFluxElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numFluxElems; ++k[1] )
      {
        // clear working arrays
        real64 densMean[numEqn]{};
        real64 dDensMean_dP[numEqn][2]{};

        real64 presGrad[numEqn]{};
        real64 dPresGrad_dP[numEqn][2]{};

        real64 gravHead[numEqn]{};
        real64 dGravHead_dP[numEqn][2]{};

        real64 capGrad[numEqn]{};
        // real64 capPresIC[numEqn][2]{};
        // real64 jFMultiplier[numEqn][2]{};
        real64 dCapGrad_dP[numEqn][2]{};
        real64 dCapGrad_dS[numEqn][2]{};

        real64 fluxVal[numEqn]{};
        real64 dFlux_dP[numEqn][2]{};
        real64 dFlux_dS[numEqn][2]{};

        real64 mobility[numEqn]{};
        real64 dMob_dP[numEqn][2]{};
        real64 dMob_dS[numEqn][2]{};

        real64 density2[numEqn]{};
        real64 dDens_dP2[numEqn][2]{};
        real64 viscosity[numEqn]{};
        real64 dVisc_dP[numEqn][2]{};
        real64 gravCoef2[numEqn]{};
        real64 gravCoefHat[numEqn]{};

        real64 uT = 0;
        // real64 total_mobility = 0;
        real64 duT_dP[numEqn]{};
        real64 duT_dS[numEqn]{};

        real64 potGrad_ip[numEqn]{};
        real64 alpha_ip[numEqn]{};

        real64 const trans[2] = { stack.transmissibility[connectionIndex][0], stack.transmissibility[connectionIndex][1] };
        real64 const dTrans_dP[2] = { stack.dTrans_dPres[connectionIndex][0], stack.dTrans_dPres[connectionIndex][1] };

        real64 const transHat[2] = { stack.transmissibilityHat[connectionIndex][0], stack.transmissibilityHat[connectionIndex][1] * -1.0};
        real64 const dTransHat_dP[2] = { stack.dTransHat_dPres[connectionIndex][0], stack.dTransHat_dPres[connectionIndex][1] * -1.0};

        // cell indices
        localIndex const seri[2]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        stdVector< real64 > saturations = {m_phaseVolFrac[seri[0]][sesri[0]][sei[0]][0], m_phaseVolFrac[seri[1]][sesri[1]][sei[1]][0] };
        stdVector< real64 > pressures = {m_pres[seri[0]][sesri[0]][sei[0]], m_pres[seri[1]][sesri[1]][sei[1]] };
        // bool isJfunction = 0;

        // loop over phases
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {
          // calculate quantities on primary connected cells
          integer denom = 0;
          for( integer ke = 0; ke < 2; ++ke )
          {
            // density
            bool const phaseExists = (m_phaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][ip] > 0);
            if( m_checkPhasePresenceInGravity && !phaseExists )
            {
              continue;
            }

            real64 const density  = m_dens[seri[ke]][sesri[ke]][sei[ke]][0][ip];         // r = rho1 || rho2
            real64 const dDens_dP = m_dDens_dPres[seri[ke]][sesri[ke]][sei[ke]][0][ip][Deriv::dP]; // dr/dP = dr1/dP1 || dr2/dP
            // average density and derivatives
            densMean[ip] += density; // rho = (rho1 + rho2)
            dDensMean_dP[ip][ke] = dDens_dP; // drho/dP = { (dr1/dP1) , (dr2/dP2) }
            denom++;
          }

          if( denom > 1 )
          {
            densMean[ip] /= denom; // rho = (rho1 + rho2) / denom
            for( integer ke = 0; ke < 2; ++ke )
            {
              dDensMean_dP[ip][ke] /= denom; // drho/dP = { (dr1/dP1) / denom , (dr2/dP2) / denom }
            }
          }

          //***** calculation of flux *****

          // compute potential difference
          real64 potScale = 0.0;
          real64 dPresGrad_dTrans = 0.0;
          real64 dGravHead_dTrans = 0.0;
          real64 dCapGrad_dTrans = 0.0;
          gravCoefHat[0] = 0;
          gravCoefHat[1] = 0;
          constexpr int signPotDiff[2] = {1, -1};

          for( integer ke = 0; ke < 2; ++ke )
          {
            localIndex const er  = seri[ke];
            localIndex const esr = sesri[ke];
            localIndex const ei  = sei[ke];

            real64 const pressure = m_pres[er][esr][ei]; // P = P1 || P2
            presGrad[ip] += trans[ke] * pressure;  // DPv = T (P1 - P2)
            dPresGrad_dTrans += signPotDiff[ke] * pressure; // dDPv/dT = (P1 - P2)
            dPresGrad_dP[ip][ke] = trans[ke];      // dDPv/dP = { T , -T }

            real64 const gravD = trans[ke] * m_gravCoef[er][esr][ei]; // D = T g z1 || -T g z2
            real64 pot = trans[ke] * pressure - densMean[ip] * gravD; // Phi = T P1 - rho T g z1 || -T P2 + rho T g z2
            gravCoefHat[0] += m_gravCoef[er][esr][ei] * 0.5;
            gravCoefHat[1] += m_gravCoef[er][esr][ei] * 0.5;

            gravCoef2[ke] = m_gravCoef[er][esr][ei];

            gravHead[ip] += densMean[ip] * gravD;                              // DPg = rho (T g z1 - T g z2) = T rho g (z1 - z2)
            dGravHead_dTrans += signPotDiff[ke] * densMean[ip] * m_gravCoef[er][esr][ei]; // dDPg/dT = rho g z1 - rho g z2 = rho g (z1 - z2)

            for( integer i = 0; i < 2; ++i )
            {
              dGravHead_dP[ip][i] += dDensMean_dP[ip][i] * gravD; // dDPg/dP = {drho/dP1 * T g (z1 - z2) , drho/dP2 * T g (z1 - z2)}
            }

            if( m_hasCapPressure ) // check sign convention
            {
              real64 const capPres = m_phaseCapPressure[er][esr][ei][0][ip]; // Pc = Pc1 || Pc2
              // jFMultiplier[ip][ke] = m_jFuncMultiplier[er][esr][ei][0];

              // capPresIC[ip][ke] = capPres;
              dCapGrad_dTrans -= signPotDiff[ke] * capPres;      // dDPc/dT = (-Pc1 + Pc2)
              pot -= trans[ke] * capPres;                        // Phi = T P1 - rho T g z1 - T Pc1 || -T P2 + rho T g z2 + T
              // Pc2
              capGrad[ip] -= trans[ke] * capPres;                // DPc = T (-Pc1 + Pc2)
            }

            potScale = fmax( potScale, fabs( pot ) ); // maxPhi = Phi1 > Phi2 ? Phi1 : Phi2
          }


          for( integer ke = 0; ke < 2; ++ke )
          {
            dPresGrad_dP[ip][ke] += dTrans_dP[ke] * dPresGrad_dTrans; // dDPv/dP = { T + dT/dP1 * (P1 - P2) , -T + dT/dP2 * (P1 - P2)}
            dGravHead_dP[ip][ke] += dTrans_dP[ke] * dGravHead_dTrans; // dDPg/dP = { drho/dP1 * T g (z1 - z2) + dT/dP1 * rho g (z1 - z2) ,
            //             drho/dP2 * T g (z1 - z2) + dT/dP2 * rho g (z1 - z2) }
            if( m_hasCapPressure )
            {
              real64 const dCapPres_dS = m_dPhaseCapPressure_dPhaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][0][ip][ip]; // dPc/dS = dPc1/dS1 ||
              // dPc2/dS2
              dCapGrad_dP[ip][ke] += dTrans_dP[ke] * dCapGrad_dTrans;                                   // dDPc/dP = { dT/dP1 *
              // (-Pc1 + Pc2) ,
              //             dT/dP2 *
              // (-Pc1 + Pc2) }
              dCapGrad_dS[ip][ke] -= trans[ke] * dCapPres_dS;                                           // dDPc/dS = { -T *
              // dPc1/dS1 , T *
              // dPc2/dS2 }
            }
          }

          // *** upwinding ***
          // compute potential gradient
          real64 potGrad = presGrad[ip] - gravHead[ip]; // DPhi = T (P1 - P2) - T rho g (z1 - z2)
          if( m_hasCapPressure )
          {
            potGrad += capGrad[ip]; // DPhi = T (P1 - P2) - T rho g (z1 - z2) + T (-Pc1 + Pc2)
          }

          // compute upwinding tolerance
          real64 constexpr upwRelTol = 1e-8;
          real64 const upwAbsTol = fmax( potScale * upwRelTol, LvArray::NumericLimits< real64 >::epsilon ); // abstol = maxPhi * tol > eps ?
          // maxPhi * tol : eps

          // decide mobility coefficients - smooth variation in [-upwAbsTol; upwAbsTol]
          real64 alpha = ( potGrad + upwAbsTol ) / ( 2 * upwAbsTol ); // alpha = (DPhi + abstol) / abstol / 2

          // choose upstream cell
          if( alpha <= 0.0 || alpha >= 1.0 ) // no smoothing needed
          {
            localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) ); // 1 upwind -> k_up = 0 || 2 upwind -> k_up = 1

            mobility[ip] = m_mob[seri[k_up]][sesri[k_up]][sei[k_up]][ip];          // M = Mupstream
            density2[ip] = m_dens[seri[k_up]][sesri[k_up]][sei[k_up]][0][ip];         // r = rho1 || rho2
            dDens_dP2[ip][k_up] = m_dDens_dPres[seri[k_up]][sesri[k_up]][sei[k_up]][0][ip][Deriv::dP]; // dr/dP = dr1/dP1 || dr2/dP
            viscosity[ip]  = m_visc[seri[k_up]][sesri[k_up]][sei[k_up]][0][ip];
            dVisc_dP[ip][k_up] = m_dVisc_dPres[seri[k_up]][sesri[k_up]][sei[k_up]][0][ip][Deriv::dP];
            dMob_dP[ip][k_up] = m_dMob[seri[k_up]][sesri[k_up]][sei[k_up]][ip][Deriv::dP]; // dM/dP = {dM/dP1 , 0} OR {0 , dM/dP2}
            dMob_dS[ip][k_up] = m_dMob[seri[k_up]][sesri[k_up]][sei[k_up]][ip][Deriv::dS]; // dM/dS = {dM/dS1 , 0} OR {0 , dM/dS2}
          }
          else // perform smoothing
          {
            real64 const mobWeights[2] = { alpha, 1.0 - alpha };
            for( integer ke = 0; ke < 2; ++ke )
            {
              mobility[ip] += mobWeights[ke] * m_mob[seri[ke]][sesri[ke]][sei[ke]][ip];  // M = alpha * M1 + (1 - alpha) * M2
              density2[ip] += mobWeights[ke] * m_dens[seri[ke]][sesri[ke]][sei[ke]][0][ip];       // r = rho1 || rho2
              dDens_dP2[ip][ke] = mobWeights[ke] * m_dDens_dPres[seri[ke]][sesri[ke]][sei[ke]][0][ip][Deriv::dP]; // dr/dP = dr1/dP1 ||
              // dr2/dP
              viscosity[ip]  = m_visc[seri[ke]][sesri[ke]][sei[ke]][0][ip];
              dVisc_dP[ip][ke] = m_dVisc_dPres[seri[ke]][sesri[ke]][sei[ke]][0][ip][Deriv::dP];
              dMob_dP[ip][ke] = mobWeights[ke] * m_dMob[seri[ke]][sesri[ke]][sei[ke]][ip][Deriv::dP]; // dM/dP = {alpha * dM1/dP1 , (1 -
              // alpha) * dM2/dP2}
              dMob_dS[ip][ke] = mobWeights[ke] * m_dMob[seri[ke]][sesri[ke]][sei[ke]][ip][Deriv::dS]; // dM/dP = {alpha * dM1/dS1 , (1 -
              // alpha) * dM2/dS2}
            }
          }

          // pressure gradient depends on all points in the stencil
          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] += dPresGrad_dP[ip][ke]; // dF/dP = { T + dT/dP1 * (P1 - P2) ,
            //          -T + dT/dP2 * (P1 - P2) }
          }

          // gravitational head depends only on the two cells connected (same as mean density)
          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] -= dGravHead_dP[ip][ke]; // dF/dP = { T + dT/dP1 * (P1 - P2) - drho/dP1 * T g (z1 - z2) - dT/dP1 * rho g (z1 -
            // z2) ,
            //          -T + dT/dP2 * (P1 - P2) - drho/dP2 * T g (z1 - z2) - dT/dP2 * rho g (z1 -
            // z2) }
          }

          // capillary pressure contribution
          if( m_hasCapPressure )
          {
            for( integer ke = 0; ke < 2; ++ke )
            {
              dFlux_dP[ip][ke] += dCapGrad_dP[ip][ke]; // dF/dP = { T + dT/dP1 * (P1 - P2) - drho/dP1 * T g (z1 - z2) - dT/dP1 * rho g (z1
              // - z2) + dT/dP1 * (-Pc1 + Pc2) ,
              //          -T + dT/dP2 * (P1 - P2) - drho/dP2 * T g (z1 - z2) - dT/dP2 * rho g (z1
              // - z2) + dT/dP2 * (-Pc1 + Pc2) }

              dFlux_dS[ip][ke] += dCapGrad_dS[ip][ke]; // dF/dS = { T * -dPc/dS1 , T * dPc/dS2 }
            }
          }

          // compute the flux and derivatives using upstream cell mobility
          fluxVal[ip] = mobility[ip] * potGrad; // F = M * DPhi

          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] *= mobility[ip]; // dF/dP = { M [ T + dT/dP1 * (P1 - P2) - drho/dP1 * T g (z1 - z2) - dT/dP1 * rho g (z1 -
            // z2) + dT/dP1 * (-Pc1 + Pc2)] ,
            //           M [-T + dT/dP2 * (P1 - P2) - drho/dP2 * T g (z1 - z2) - dT/dP2 * rho g (z1 -
            // z2) + dT/dP2 * (-Pc1 + Pc2)] }

            dFlux_dS[ip][ke] *= mobility[ip]; // dF/dS = { M [T * -dPc/dS1] , M [T * dPc/dS2] }
          }

          // add contribution from upstream cell mobility derivatives
          for( integer ke = 0; ke < 2; ++ke )
          {

            real64 dMob_dP2 = mobility[ip] / density2[ip] * (-dVisc_dP[ip][ke] / viscosity[ip]);

            duT_dP[ke] += dFlux_dP[ip][ke] / density2[ip] + dMob_dP2 * potGrad;

            dFlux_dP[ip][ke] += dMob_dP[ip][ke] * potGrad; // dF/dP = { M [ T + dT/dP1 * (P1 - P2) - drho1/dP * T g (z1 - z2) - dT/dP1 *
            // rho g (z1 - z2) + dT/dP1 * (-Pc1 + Pc2)] + dM/dP1 * DPhi ,
            //           M [-T + dT/dP2 * (P1 - P2) - drho2/dP * T g (z1 - z2) - dT/dP2 *
            // rho g (z1 - z2) + dT/dP2 * (-Pc1 + Pc2)] + dM/dP2 * DPhi }
            dFlux_dS[ip][ke] += dMob_dS[ip][ke] * potGrad; // dF/dS = { M [T * -dPc/dS1] + dM/dS1 * DPhi , M [T * dPc/dS2] + dM/dS2 * DPhi
            // }
          }


          uT += fluxVal[ip] / density2[ip];

          // add contribution from upstream cell mobility derivatives
          for( integer ke = 0; ke < 2; ++ke )
          {
//  duT_dP[ke] += dFlux_dP[ip][ke] - fluxVal[ip] * dDens_dP2[ip][ke] / density2[ip];
//  // duT_dP[ke] += dFlux_dP[ip][ke];

//  duT_dP[ke]  /= density2[ip];

            duT_dS[ke] += dFlux_dS[ip][ke] / density2[ip];
            // duT_dS[ke] += dFlux_dS[ip][ke];

          }

          potGrad_ip[ip] = potGrad;
          alpha_ip[ip] = alpha;

        } // loop over phases


        // this determines whether the local solver is needed becuase of heterogeneous capillary pressure regions
        // bool notOnInterface = std::fabs( jFMultiplier[0][0] - jFMultiplier[0][1] ) < 1  && std::fabs( jFMultiplier[1][0] - jFMultiplier[1][1] ) < 1;
       bool notOnInterface = !connectorHasInterfaceConditionQ;
        if( notOnInterface )
        {
          for( integer ip = 0; ip < 2; ++ip )
          {
            // populate local flux vector and derivatives
            stack.localFlux[k[0]*numEqn + ip] += m_dt * fluxVal[ip];
            stack.localFlux[k[1]*numEqn + ip] -= m_dt * fluxVal[ip];

            for( integer ke = 0; ke < 2; ++ke )
            {
              // pressure
              localIndex const localDofIndexPres = k[ke] * numDof;
              stack.localFluxJacobian[k[0]*numEqn + ip][localDofIndexPres] += m_dt * dFlux_dP[ip][ke];
              stack.localFluxJacobian[k[1]*numEqn + ip][localDofIndexPres] -= m_dt * dFlux_dP[ip][ke];

              // saturation
              localIndex const localDofIndexSat = k[ke] * numDof + 1;
              stack.localFluxJacobian[k[0]*numEqn + ip][localDofIndexSat] += m_dt * dFlux_dS[ip][ke];
              stack.localFluxJacobian[k[1]*numEqn + ip][localDofIndexSat] -= m_dt * dFlux_dS[ip][ke];
            }

            // Customize the kernel with this lambda
            kernelOp( k, seri, sesri, sei, connectionIndex, alpha_ip[ip], mobility, potGrad_ip[ip], fluxVal, dFlux_dP, dFlux_dS );

          }
        }
        else
        {


          bool converged = 0;


          // clear working arrays
          real64 halfFluxVal[numEqn][2]{};
          // real64 dhalfFlux1_dP[numEqn][2]{};
          // real64 dhalfFlux1_dS[numEqn][2]{};
          // real64 dhalfFlux2_dP[numEqn][2]{};
          // real64 dhalfFlux2_dS[numEqn][2]{};
          // real64 dhalfFlux_duT[numEqn][2]{};
          // real64 dhalfFlux_dpc[numEqn][2]{};


          // stdVector< real64 > JFMultipliers = {jFMultiplier[0][0], jFMultiplier[0][1]};
          stdVector< real64 > JFMultipliers = {0.0, 0.0};
          stdVector< real64 > trappedSats1 = {0.0, 0.0};
          stdVector< real64 > trappedSats2 = {0.0, 0.0};
          stdVector< real64 > transHats = {transHat[0], transHat[1]};
          stdVector< real64 > dTransHats_dP = {dTransHat_dP[0], dTransHat_dP[1]};
          stdVector< real64 > gravCoefHats = {gravCoefHat[0], gravCoefHat[1]};
          stdVector< real64 > gravCoefs = {gravCoef2[0], gravCoef2[1]};
          stdVector< real64 > cellCenterDuTdS = {duT_dP[0], duT_dP[1], duT_dS[0], duT_dS[1]};
          stdVector< real64 > cellCenterDens = {density2[0], density2[1]};
          stdVector< real64 > cellCenterDens_dP = {dDens_dP2[0][0], dDens_dP2[0][1], dDens_dP2[1][0], dDens_dP2[1][1]};
          // std::vector< RelativePermeabilityBase * > relPerms = {std::get< 0 >( m_interfaceConstitutivePairs_temp ), std::get< 0 >( m_interfaceConstitutivePairs_temp )};
          // std::vector< CapillaryPressureBase * > capPressures = {std::get< 1 >( m_interfaceConstitutivePairs_temp ), std::get< 1 >( m_interfaceConstitutivePairs_temp )};
          // std::vector< TwoPhaseImmiscibleFluid * > fluids = {std::get< 2 >( m_interfaceConstitutivePairs_temp ), std::get< 2 >( m_interfaceConstitutivePairs_temp )};

          // auto const & pairArray = m_interfaceConstitutivePairs[0];
          localIndex const surfaceRegionIndex = m_interfaceRegionByConnector.at(iconn);
auto const & pairArray = m_interfaceConstitutivePairs[surfaceRegionIndex];

std::vector< constitutive::RelativePermeabilityBase * > relPerms = {
  std::get<0>( pairArray[0] ),
  std::get<0>( pairArray[1] )
};

std::vector< constitutive::CapillaryPressureBase * > capPressures = {
  std::get<1>( pairArray[0] ),
  std::get<1>( pairArray[1] )
};

std::vector< constitutive::TwoPhaseImmiscibleFluid * > fluids = {
  std::get<2>( pairArray[0] ),
  std::get<2>( pairArray[1] )
};

          stdVector< real64 > phi = {halfFluxVal[0][0], halfFluxVal[0][1]};
          stdVector< real64 > grad_phi_P = {0.0, 0.0, 0.0, 0.0};
          stdVector< real64 > grad_phi_S = {0.0, 0.0, 0.0, 0.0};

          local_solver( uT, saturations, pressures, JFMultipliers, trappedSats1, trappedSats2, transHats, dTransHats_dP, gravCoefHats, gravCoefs, 
            cellCenterDuTdS, cellCenterDens, cellCenterDens_dP, relPerms, capPressures, fluids, phi, grad_phi_P, grad_phi_S, converged );
          

             fluxVal[0] = phi[0];
             fluxVal[1] = phi[1];
             dFlux_dP[0][0] = grad_phi_P[0];
             dFlux_dP[0][1] = grad_phi_P[1];
             dFlux_dP[1][0] = grad_phi_P[2];
             dFlux_dP[1][1] = grad_phi_P[3];
             dFlux_dS[0][0] = grad_phi_S[0];
             dFlux_dS[0][1] = grad_phi_S[1];
             dFlux_dS[1][0] = grad_phi_S[2];
             dFlux_dS[1][1] = grad_phi_S[3];

          // Global residual and jacobian update:
          for( integer ip = 0; ip < m_numPhases; ++ip )
          {
            // populate local flux vector and derivatives
            stack.localFlux[k[0]*numEqn + ip] += m_dt * fluxVal[ip];
            stack.localFlux[k[1]*numEqn + ip] -= m_dt * fluxVal[ip];

            for( integer ke = 0; ke < 2; ++ke )
            {
              // pressure
              localIndex const localDofIndexPres = k[ke] * numDof;
              stack.localFluxJacobian[k[0]*numEqn + ip][localDofIndexPres] += m_dt * dFlux_dP[ip][ke];
              stack.localFluxJacobian[k[1]*numEqn + ip][localDofIndexPres] -= m_dt * dFlux_dP[ip][ke];

              // saturation
              localIndex const localDofIndexSat = k[ke] * numDof + 1;
              stack.localFluxJacobian[k[0]*numEqn + ip][localDofIndexSat] += m_dt * dFlux_dS[ip][ke];
              stack.localFluxJacobian[k[1]*numEqn + ip][localDofIndexSat] -= m_dt * dFlux_dS[ip][ke];
            }

            // Customize the kernel with this lambda
            kernelOp( k, seri, sesri, sei, connectionIndex, alpha_ip[ip], mobility, potGrad_ip[ip], fluxVal, dFlux_dP, dFlux_dS );

          } // loop over phases
        } // end of else for interface conditions

        connectionIndex++;
      }
    } // loop over connection elements
  }

protected:

  /// Reference to the capillary pressure wrapper
  // CAPPRESWRAPPER const m_capPressureWrapper;
  // RELPERMWRAPPER const m_relPermWrapper;
  string_array const m_interfaceFaceSetNames;
  stdVector< std::array< std::tuple< constitutive::RelativePermeabilityBase *,
                                     constitutive::CapillaryPressureBase *,
                                     constitutive::TwoPhaseImmiscibleFluid * >, 2 > >  const m_interfaceConstitutivePairs;
  unordered_map< localIndex, localIndex > const m_interfaceRegionByConnector;
  // std::tuple< constitutive::RelativePermeabilityBase *,
  //             constitutive::CapillaryPressureBase *,
  //             constitutive::TwoPhaseImmiscibleFluid * > const m_interfaceConstitutivePairs_temp;

};



/****************************************** */

/**
 * @class FluxComputeKernelFactory
 */
class FluxComputeKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const hasCapPressure,
                   integer const useTotalMassEquation,
                   integer const checkPhasePresenceInGravity,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_EQN = 2;
    integer constexpr NUM_DOF = 2;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = FluxComputeKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
    typename kernelType::ImmiscibleMultiphaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::MultiphaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );

    kernelType kernel( numPhases, rankOffset, stencilWrapper, dofNumberAccessor,
                       flowAccessors, fluidAccessors, capPressureAccessors, permAccessors,
                       dt, localMatrix, localRhs, hasCapPressure, useTotalMassEquation,
                       checkPhasePresenceInGravity );
    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @tparam CAPPRESWRAPPER the type of the capillary pressure wrapper
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] capPresWrapper reference to the capillary pressure wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  // template< typename POLICY, typename STENCILWRAPPER, typename CAPPRESWRAPPER, typename RELPERMWRAPPER >
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const hasCapPressure,
                   integer const useTotalMassEquation,
                   integer const checkPhasePresenceInGravity,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                  //  CAPPRESWRAPPER const & capPresWrapper,
                  //  RELPERMWRAPPER const & relPermWrapper,
                   string_array const & interfaceFaceSetNames,
                   stdVector< std::array< std::tuple< constitutive::RelativePermeabilityBase *,
                                                      constitutive::CapillaryPressureBase *,
                                                      constitutive::TwoPhaseImmiscibleFluid * >, 2 > > const & interfaceConstitutivePairs,
                   unordered_map< localIndex, localIndex > const & interfaceRegionByConnector,
                  //  std::tuple< constitutive::RelativePermeabilityBase *,
                  //              constitutive::CapillaryPressureBase *,
                  //              constitutive::TwoPhaseImmiscibleFluid * > const & interfaceConstitutivePairs_temp,
                   ElementSubRegionBase const & subRegion,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_EQN = 2;
    integer constexpr NUM_DOF = 2;
    localIndex const domainSize = subRegion.size();
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    // using kernelType = FluxComputeInterfaceConditionKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER, CAPPRESWRAPPER, RELPERMWRAPPER >;
        using kernelType = FluxComputeInterfaceConditionKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
    typename kernelType::ImmiscibleMultiphaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::MultiphaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );

    // kernelType kernel( numPhases, rankOffset, stencilWrapper, capPresWrapper, relPermWrapper, dofNumberAccessor,
    //                    flowAccessors, fluidAccessors, capPressureAccessors, permAccessors,
    //                    dt, localMatrix, localRhs, hasCapPressure, useTotalMassEquation,
    //                    checkPhasePresenceInGravity, interfaceFaceSetNames, interfaceConstitutivePairs, interfaceRegionByConnector, interfaceConstitutivePairs_temp, domainSize );
    kernelType kernel( numPhases, rankOffset, stencilWrapper, dofNumberAccessor,
                       flowAccessors, fluidAccessors, capPressureAccessors, permAccessors,
                       dt, localMatrix, localRhs, hasCapPressure, useTotalMassEquation,
                       checkPhasePresenceInGravity, interfaceFaceSetNames, interfaceConstitutivePairs, interfaceRegionByConnector, domainSize );
    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};


/******************************** AccumulationKernel ********************************/

enum class KernelFlags
{
  TotalMassEquation = 1 << 0,  // 1

  /// Add more flags like that if needed:
  // Flag2 = 1 << 1, // 2
  // Flag3 = 1 << 2, // 4
  // Flag4 = 1 << 3, // 8
  // Flag5 = 1 << 4, // 16
  // Flag6 = 1 << 5, // 32
  // Flag7 = 1 << 6, // 64
  // Flag8 = 1 << 7  //128
};
/**
 * @class AccumulationKernel
 * @tparam NUM_EQN number of fluid phases
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */

template< integer NUM_EQN, integer NUM_DOF >
class AccumulationKernel
{
public:

  /// Number of fluid phases
  integer const m_numPhases;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_EQN;

  /**
   * @brief Constructor
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[inout] kernelFlags the kernel options
   */
  AccumulationKernel( localIndex const numPhases,
                      globalIndex const rankOffset,
                      string const dofKey,
                      ElementSubRegionBase const & subRegion,
                      constitutive::TwoPhaseImmiscibleFluid const & fluid,
                      constitutive::CoupledSolidBase const & solid,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs,
                      BitFlags< KernelFlags > const kernelFlags )
    : m_numPhases( numPhases ),
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity( solid.getPorosity() ),
    m_dPoro_dPres( solid.getDporosity_dPressure() ),
    m_phaseVolFrac( subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseMass_n( subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs ),
    m_kernelFlags( kernelFlags )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    // Pore volume information (used by both accumulation)

    /// Pore volume at time n+1
    real64 poreVolume = 0.0;

    /// Derivative of pore volume with respect to pressure
    real64 dPoreVolume_dPres = 0.0;

    // Residual information

    /// Index of the local row corresponding to this element
    localIndex localRow = -1;

    /// Indices of the matrix rows/columns corresponding to the dofs in this element
    globalIndex dofIndices[numDof]{};

    /// C-array storage for the element local residual vector (all equations except volume balance)
    real64 localResidual[numEqn]{};

    /// C-array storage for the element local Jacobian matrix (all equations except volume balance, all dofs)
    real64 localJacobian[numEqn][numDof]{};

  };

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    // initialize the pore volume
    stack.poreVolume = m_volume[ei] * m_porosity[ei][0];
    stack.dPoreVolume_dPres = m_volume[ei] * m_dPoro_dPres[ei][0];

    // set row index and degrees of freedom indices for this element
    stack.localRow = m_dofNumber[ei] - m_rankOffset;
    for( integer idof = 0; idof < numDof; ++idof )
    {
      stack.dofIndices[idof] = m_dofNumber[ei] + idof;
    }
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    constexpr int sign[2] = {1, -1};
    // ip - index of phase/component whose conservation equation is assembled
    // (i.e. row number in local matrix)
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      real64 const phaseMass = stack.poreVolume * m_phaseVolFrac[ei][ip] * m_phaseDens[ei][0][ip];
      real64 const phaseMass_n = m_phaseMass_n[ei][ip];

      stack.localResidual[ip] += phaseMass - phaseMass_n;

      real64 const dPhaseMass_dP = stack.dPoreVolume_dPres * m_phaseVolFrac[ei][ip] * m_phaseDens[ei][0][ip]
                                   + stack.poreVolume * m_phaseVolFrac[ei][ip] * m_dPhaseDens[ei][0][ip][0];
      stack.localJacobian[ip][0] += dPhaseMass_dP;

      real64 const dPhaseMass_dS = stack.poreVolume * m_phaseDens[ei][0][ip];
      stack.localJacobian[ip][1] += sign[ip] * dPhaseMass_dS;
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const GEOS_UNUSED_PARAM( ei ),
                 StackVariables & stack ) const
  {
    using namespace compositionalMultiphaseUtilities;

    if( m_kernelFlags.isSet( KernelFlags::TotalMassEquation ) )
    {
      // apply equation/variable change transformation to the component mass balance equations
      real64 work[numDof]{};
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numPhases, numDof, stack.localJacobian, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( m_numPhases, stack.localResidual );
    }

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numPhase - 1)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    integer const numRows = m_numPhases;
    for( integer i = 0; i < numRows; ++i )
    {
      m_localRhs[stack.localRow + i] += stack.localResidual[i];
      m_localMatrix.addToRow< serialAtomic >( stack.localRow + i,
                                              stack.dofIndices,
                                              stack.localJacobian[i],
                                              numDof );
    }
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.elemGhostRank( ei ) >= 0 )
      {
        return;
      }

      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, stack );
      kernelComponent.computeAccumulation( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const m_volume;

  /// Views on the porosity
  arrayView2d< real64 const > const m_porosity;
  arrayView2d< real64 const > const m_dPoro_dPres;

  /// Views on the phase volume fractions
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const m_phaseVolFrac;

  /// Views on the phase densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const m_dPhaseDens;

  // View on component amount (mass) from previous time step
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > m_phaseMass_n;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  BitFlags< KernelFlags > const m_kernelFlags;
};

/**
 * @class AccumulationKernelFactory
 */
class AccumulationKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[inout] useTotalMassEquation option for using total mass equation
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numPhases,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   constitutive::TwoPhaseImmiscibleFluid const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {

    geos::immiscibleMultiphaseKernels::kernelLaunchSelectorPhaseSwitch( numPhases, [&] ( auto NP )
    {
      integer constexpr NUM_EQN = NP();
      integer constexpr NUM_DOF = NP();

      BitFlags< KernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( KernelFlags::TotalMassEquation );

      AccumulationKernel< NUM_EQN, NUM_DOF > kernel( numPhases, rankOffset, dofKey, subRegion,
                                                     fluid, solid, localMatrix, localRhs, kernelFlags );
      AccumulationKernel< NUM_EQN, NUM_DOF >::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }

};



/******************************** PhaseMobilityKernel ********************************/

/**
 * @class PhaseMobilityKernel
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase mobilities
 */
template< integer NUM_PHASE >
class PhaseMobilityKernel
{
public:

  //using Base = MultiphaseFluidAccessors::PropertyKernelBase< NUM_COMP >;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  PhaseMobilityKernel( ObjectManagerBase & subRegion,
                       TwoPhaseImmiscibleFluid const & fluid,
                       RelativePermeabilityBase const & relperm )
    :
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseVisc( fluid.phaseViscosity() ),
    m_dPhaseVisc( fluid.dPhaseViscosity() ),
    m_phaseRelPerm( relperm.phaseRelPerm() ),
    m_dPhaseRelPerm_dPhaseVolFrac( relperm.dPhaseRelPerm_dPhaseVolFraction() ),
    m_phaseMob( subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMobility >() ),
    m_dPhaseMob( subRegion.getField< fields::immiscibleMultiphaseFlow::dPhaseMobility >() )
  {}

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      kernelComponent.compute( ei );
    } );
  }

  /**
   * @brief Compute the phase mobilities in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseMobilityKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                FUNC && phaseMobilityKernelOp = NoOpFunc{} ) const
  {
    using Deriv = immiscibleFlow::DerivativeOffset;

    arraySlice1d< real64, immiscibleFlow::USD_PHASE - 1 > const phaseMob = m_phaseMob[ei];
    arraySlice2d< real64, immiscibleFlow::USD_PHASE_DS - 1 > const dPhaseMob = m_dPhaseMob[ei];
    constexpr int sign[2] = {1, -1};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      real64 const density = m_phaseDens[ei][0][ip];
      real64 const dDens_dP = m_dPhaseDens[ei][0][ip][Deriv::dP];
      real64 const viscosity = m_phaseVisc[ei][0][ip];
      real64 const dVisc_dP = m_dPhaseVisc[ei][0][ip][Deriv::dP];

      real64 const relPerm = m_phaseRelPerm[ei][0][ip];

      real64 const mobility = relPerm * density / viscosity;

      phaseMob[ip] = mobility;
      dPhaseMob[ip][Deriv::dP] = mobility * (dDens_dP / density - dVisc_dP / viscosity);

      // for( integer jp = 0; jp < numPhase-1; ++jp )
      // {
      // Derivative matrix is currently diagonal. Implementation below handles missing off-diagonal entry.
      real64 const dRelPerm_dS = sign[ip] * m_dPhaseRelPerm_dPhaseVolFrac[ei][0][ip][ip];
      dPhaseMob[ip][Deriv::dS] = dRelPerm_dS * density / viscosity;
      // }


      // call the lambda in the phase loop to allow the reuse of the relperm, density, viscosity, and mobility
      // possible use: assemble the derivatives wrt temperature
      phaseMobilityKernelOp( ip, phaseMob[ip], dPhaseMob[ip] );
    }
  }

protected:

  // inputs

  /// Views on the phase densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseDens;

  //arrayView2d< real64 const, immiscibleFlow::USD_PHASE > m_phaseDens;
  //arrayView2d< real64 const, immiscibleFlow::USD_PHASE > m_dPhaseDens;

  /// Views on the phase viscosities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseVisc;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseVisc;
  //arrayView2d< real64 const, immiscibleFlow::USD_PHASE > m_phaseVisc;
  //arrayView2d< real64 const, immiscibleFlow::USD_PHASE > m_dPhaseVisc;

  /// Views on the phase relative permeabilities
  arrayView3d< real64 const, relperm::USD_RELPERM > m_phaseRelPerm;
  arrayView4d< real64 const, relperm::USD_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;

  // outputs

  /// Views on the phase mobilities
  arrayView2d< real64, immiscibleFlow::USD_PHASE > m_phaseMob;
  arrayView3d< real64, immiscibleFlow::USD_PHASE_DS > m_dPhaseMob;
};

/**
 * @class PhaseMobilityKernelFactory
 */
class PhaseMobilityKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] relperm the relperm model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numPhase,
                   ObjectManagerBase & subRegion,
                   TwoPhaseImmiscibleFluid const & fluid,
                   RelativePermeabilityBase const & relperm )
  {
    if( numPhase == 2 )
    {
      PhaseMobilityKernel< 2 > kernel( subRegion, fluid, relperm );
      PhaseMobilityKernel< 2 >::template launch< POLICY >( subRegion.size(), kernel );
    }
  }
};


struct FluidUpdateKernel
{
  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k] );
      }
    } );
  }
};

//******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 *
 * @brief Define the interface for the property kernel in charge of computing the residual norm
 */

/// Compile time value for the number of norms to compute
static constexpr integer numNorm = 1;

class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< numNorm >
{
public:


  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< numNorm >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      integer const numPhases,
                      ElementSubRegionBase const & subRegion,
                      constitutive::CoupledSolidBase const & solid,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_numPhases( numPhases ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity_n( solid.getPorosity_n() ),
    m_phaseMass_n( subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >() )
  {}

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 massNormalizer = 0;
    for( integer idof = 0; idof < m_numPhases; ++idof )
    {
      massNormalizer += LvArray::math::max( m_minNormalizer, m_phaseMass_n[ei][idof] );
    }

    for( integer idof = 0; idof < m_numPhases; ++idof )
    {
      real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / massNormalizer;
      if( valMass > stack.localValue[0] )
      {
        stack.localValue[0] = valMass;
      }
    }

    // for( integer idof = 0; idof < m_numPhases; ++idof )
    // {
    //   real64 const massNormalizer = LvArray::math::max( m_minNormalizer, m_phaseMass_n[ei][idof] );
    //   real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / massNormalizer;
    //   if( valMass > stack.localValue[0] )
    //   {
    //     stack.localValue[0] = valMass;
    //   }
    // }

  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    for( integer idof = 0; idof < m_numPhases; ++idof )
    {
      real64 const massNormalizer = LvArray::math::max( m_minNormalizer, m_phaseMass_n[ei][idof] );
      stack.localValue[0] += m_localResidual[stack.localRow + idof] * m_localResidual[stack.localRow + idof];
      stack.localNormalizer[0] += massNormalizer;
    }
  }


protected:

  /// Number of fluid phases
  integer const m_numPhases;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on porosity at the previous converged time step
  arrayView2d< real64 const > const m_porosity_n;
  /// View on mass at the previous converged time step
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const m_phaseMass_n;
};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[in] solid the solid model
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( physicsSolverBaseKernels::NormType const normType,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   constitutive::CoupledSolidBase const & solid,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[1],
                   real64 (& residualNormalizer)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, numPhases, subRegion, solid, minNormalizer );
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else  // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};



}  // namespace immiscible multiphasekernels


}  // namespace geos

 #endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_MULTIPHASEKERNELS_HPP
