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
 #include "constitutive/capillaryPressure/TableCapillaryPressureHysteresis.hpp"
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
constexpr bool ENABLE_LOCAL_SOLVER_DEBUG = false;

static void local_solver( real64 uT, stdVector< real64 > const & saturations, stdVector< real64 > const & pressures, stdVector< real64 > const & JFMultipliers,
                          stdVector< real64 > const & trappedSats1,
                          stdVector< real64 > const & trappedSats2, stdVector< fields::cappres::ModeIndexType > const & modes, stdVector< real64 > const & transHat,
                          stdVector< real64 > const & dTransHat_dP, stdVector< real64 > const & gravCoefHat, stdVector< real64 > const & gravCoef,
                          stdVector< real64 > const & cellCenterDuT, stdVector< real64 > const & cellCenterDens, stdVector< real64 > const & cellCenterDens_dP,
                          std::vector< RelativePermeabilityBase * > const & relPerms, std::vector< CapillaryPressureBase * > const & capPressures,
                          std::vector< TwoPhaseImmiscibleFluid * > const & fluids, stdVector< real64 > & phi, stdVector< real64 > & grad_phi_P, stdVector< real64 > & grad_phi_S, bool & converged,
                          stdVector< real64 > const & phaseMaxHistVolFrac1, stdVector< real64 > const & phaseMinHistVolFrac1,
                          stdVector< real64 > const & phaseMaxHistVolFrac2, stdVector< real64 > const & phaseMinHistVolFrac2,
                          real64 & warmStartPc )
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


          // std::ofstream outFile( "iterations2.csv" );


          // // Write data to the file
          // outFile << "Jacobian";
          // outFile << ",";
          // outFile << "residual";
          // outFile << ",";
          // outFile << "Fw_alpha";
          // outFile << ",";
          // outFile << "Fw_beta";
          // outFile << ",";
          // outFile << "Pc_int";
          // outFile << ",";
          // outFile << "Pc_int1";
          // outFile << ",";
          // outFile << "Pc_int2";
          // outFile << ",";
          // outFile << "Fn_alpha";
          // outFile << ",";
          // outFile << "Fn_beta";
          // outFile << ",";
          // outFile << "Vw_alpha";
          // outFile << ",";
          // outFile << "Vn_alpha";
          // outFile << ",";
          // outFile << "Vw_beta";
          // outFile << ",";
          // outFile << "Vn_beta";
          // outFile << ",";
          // outFile << "Gw_alpha";
          // outFile << ",";
          // outFile << "Gn_alpha";
          // outFile << ",";
          // outFile << "Gw_beta";
          // outFile << ",";
          // outFile << "Gn_beta";
          // outFile << ",";
          // outFile << "Cw_alpha";
          // outFile << ",";
          // outFile << "Cn_alpha";
          // outFile << ",";
          // outFile << "Cw_beta";
          // outFile << ",";
          // outFile << "Cn_beta";
          // outFile << ",";
          // outFile << "transHats0";
          // outFile << ",";
          // outFile << "transHats1";
          // outFile << ",";
          // outFile << "gravCoefHats";
          // outFile << ",";
          // outFile << "gravCoef0";
          // outFile << ",";
          // outFile << "gravCoef1";
          // outFile << ",";
          // outFile << "uT";
          // outFile << ",";
          // outFile << "duT_dP0";
          // outFile << ",";
          // outFile << "duT_dS0";
          // outFile << ",";
          // outFile << "duT_dP1";
          // outFile << ",";
          // outFile << "duT_dS1";
          // outFile << std::endl;



          // nonlinear solver's parameters
          real64 tol = 1.0e-10;
          int max_iter = 100;
          converged = 0;
          bool damping = true;

          // Local newton loop:

          // Use of the capillary pressure kernel wrapper

          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseVolFrac1( 1, 2 );
          StackArray< real64, 3, 2, constitutive::cappres::LAYOUT_CAPPRES > capPres1( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dPhaseVolFrac_dCapPres1( 1, 1, 2, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dCapPres1_dPhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > trappedVolFrac1( 1, 1, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > facePhaseVolFrac1( 1, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > faceTrappedVolFrac1( 1, 1, 2 );
          StackArray< real64, 3, 2, constitutive::cappres::LAYOUT_CAPPRES > faceCapPres1( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dfacePhaseVolFrac_dCapPres1( 1, 1, 2, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dCapPres1_dfacePhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseVolFrac2( 1, 2 );
          StackArray< real64, 3, 2, constitutive::cappres::LAYOUT_CAPPRES > capPres2( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dPhaseVolFrac_dCapPres2( 1, 1, 2, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dCapPres2_dPhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > trappedVolFrac2( 1, 1, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > facePhaseVolFrac2( 1, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > faceTrappedVolFrac2( 1, 1, 2 );
          StackArray< real64, 3, 2, constitutive::cappres::LAYOUT_CAPPRES > faceCapPres2( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dfacePhaseVolFrac_dCapPres2( 1, 1, 2, 2 );
          StackArray< real64, 4, 4, constitutive::cappres::LAYOUT_CAPPRES_DS > dCapPres2_dfacePhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 1, 2 > JFunc1( 2 );
          StackArray< real64, 1, 2 > JFunc2( 2 );


          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseMaxHistoricalVolFraction1( 1, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseMinHistoricalVolFraction1( 1, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseMaxHistoricalVolFraction2( 1, 2 );
          StackArray< real64, 2, 2, immiscibleFlow::LAYOUT_PHASE > phaseMinHistoricalVolFraction2( 1, 2 );
          StackArray< real64, 2, 2, compflow::LAYOUT_PHASE > phaseMode2PeakVolFraction1( 1, 2 );
          StackArray< real64, 2, 2, compflow::LAYOUT_PHASE > phaseMode2PeakVolFraction2( 1, 2 );
          phaseMode2PeakVolFraction1[0][0] = 0.0;
          phaseMode2PeakVolFraction1[0][1] = 0.0;
          phaseMode2PeakVolFraction2[0][0] = 0.0;
          phaseMode2PeakVolFraction2[0][1] = 0.0;
          arraySlice1d< real64, compflow::USD_PHASE - 1 > phaseMode2PeakSlice1 = phaseMode2PeakVolFraction1[0];
          arraySlice1d< real64, compflow::USD_PHASE - 1 > phaseMode2PeakSlice2 = phaseMode2PeakVolFraction2[0];

          if( phaseMaxHistVolFrac1.size() >= 2 && phaseMaxHistVolFrac1[0] >= 0.0 )
          {

            phaseMaxHistoricalVolFraction1[0][0] = phaseMaxHistVolFrac1[0];
            phaseMaxHistoricalVolFraction1[0][1] = phaseMaxHistVolFrac1[1];
            phaseMinHistoricalVolFraction1[0][0] = phaseMinHistVolFrac1[0];
            phaseMinHistoricalVolFraction1[0][1] = phaseMinHistVolFrac1[1];
          }
          else
          {

            phaseMaxHistoricalVolFraction1[0][0] = saturations[0];
            phaseMaxHistoricalVolFraction1[0][1] = 1.0 - saturations[0];
            phaseMinHistoricalVolFraction1[0][0] = saturations[0];
            phaseMinHistoricalVolFraction1[0][1] = 1.0 - saturations[0];
          }

          if( phaseMaxHistVolFrac2.size() >= 2 && phaseMaxHistVolFrac2[0] >= 0.0 )
          {

            phaseMaxHistoricalVolFraction2[0][0] = phaseMaxHistVolFrac2[0];
            phaseMaxHistoricalVolFraction2[0][1] = phaseMaxHistVolFrac2[1];
            phaseMinHistoricalVolFraction2[0][0] = phaseMinHistVolFrac2[0];
            phaseMinHistoricalVolFraction2[0][1] = phaseMinHistVolFrac2[1];
          }
          else
          {
            phaseMaxHistoricalVolFraction2[0][0] = saturations[1];
            phaseMaxHistoricalVolFraction2[0][1] = 1.0 - saturations[1];
            phaseMinHistoricalVolFraction2[0][0] = saturations[1];
            phaseMinHistoricalVolFraction2[0][1] = 1.0 - saturations[1];
          }

          // compute relative permeability for both cell centers:

          trappedVolFrac1[0][0][0] = trappedSats1[0];
          trappedVolFrac1[0][0][1] = trappedSats1[1];

          trappedVolFrac2[0][0][0] = trappedSats2[0];
          trappedVolFrac2[0][0][1] = trappedSats2[1];

          faceTrappedVolFrac1[0][0][0] = trappedSats1[0];
          faceTrappedVolFrac1[0][0][1] = trappedSats1[1];

          faceTrappedVolFrac2[0][0][0] = trappedSats2[0];
          faceTrappedVolFrac2[0][0][1] = trappedSats2[1];

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
          else if constexpr (std::is_same_v< T1, TableCapillaryPressureHysteresis >) {
            auto preComputeMode0 = modes[0];
            capPresWrapper1.compute( phaseVolFrac1[0],
                                     phaseMaxHistoricalVolFraction1[0],
                                     phaseMinHistoricalVolFraction1[0],
                                     trappedVolFrac1[0][0],
                                     capPres1[0][0],
                                     dCapPres1_dPhaseVolFrac[0][0],
                                     preComputeMode0,
                                     phaseMode2PeakSlice1 );

            if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
              std::cout << "  [BOUNDS_DEBUG] Cell0: T1=TableCapillaryPressureHysteresis"
                        << ", preComputeMode=" << static_cast< int >(preComputeMode0)
                        << ", postComputeMode=" << static_cast< int >(modes[0])
                        << ", S=" << phaseVolFrac1[0][0]
                        << ", Pc=" << capPres1[0][0][0] << std::endl;
            }

            if( modes[0] == fields::cappres::ModeIndexType::DRAINAGE ||
                modes[0] == fields::cappres::ModeIndexType::IMBIBITION )
            {
              if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) { std::cout << "  [BOUNDS_DEBUG] Cell0: DRAINAGE/IMBIBITION path" << std::endl; }
              integer const tableIdx1 = static_cast< integer >(modes[0]);
              real64 Pc1_at_S1, dPc1_at_S1;
              real64 Pc1_at_S0, dPc1_at_S0;
              real64 const S_one = 1.0;
              real64 const S_zero = 0.0;
              capPresWrapper1.computeRawTablePc( tableIdx1, S_one, Pc1_at_S1, dPc1_at_S1 );
              capPresWrapper1.computeRawTablePc( tableIdx1, S_zero, Pc1_at_S0, dPc1_at_S0 );
              Pc1_min = LvArray::math::min( Pc1_at_S1, Pc1_at_S0 );
              Pc1_max = LvArray::math::max( Pc1_at_S1, Pc1_at_S0 );

              if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) { std::cout << "  [BOUNDS_DEBUG] Cell0: Pc1_min=" << Pc1_min
                                                                   << " (table@S=1:" << Pc1_at_S1 << ", table@S=0:" << Pc1_at_S0 << ")" << std::endl; }
            }
            else
            {
              if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) { std::cout << "  [BOUNDS_DEBUG] Cell0: SCANNING path, modes[0]=" << static_cast< int >(modes[0]) << std::endl; }
              real64 const S_min_hist1 = phaseMinHistoricalVolFraction1[0][0];
              real64 const S_max_hist1 = phaseMaxHistoricalVolFraction1[0][0];

              auto scanBoundsMode1 = modes[0];

              if( modes[0] == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION ||
                  modes[0] == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE ||
                  modes[0] == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
              {
                capPresWrapper1.computeScanningCurvePcRange(
                  S_min_hist1, S_max_hist1, modes[0], Pc1_min, Pc1_max );
              }
              else
              {
                facePhaseVolFrac1[0][1] = 0.0;
                facePhaseVolFrac1[0][0] = 1.0;
                scanBoundsMode1 = modes[0];
                capPresWrapper1.compute( facePhaseVolFrac1[0],
                                         phaseMaxHistoricalVolFraction1[0],
                                         phaseMinHistoricalVolFraction1[0],
                                         faceTrappedVolFrac1[0][0],
                                         faceCapPres1[0][0],
                                         dCapPres1_dfacePhaseVolFrac[0][0],
                                         scanBoundsMode1,
                                         phaseMode2PeakSlice1 );
                Pc1_min = faceCapPres1[0][0][0];

                facePhaseVolFrac1[0][1] = 1.0;
                facePhaseVolFrac1[0][0] = 0.0;
                scanBoundsMode1 = modes[0];
                capPresWrapper1.compute( facePhaseVolFrac1[0],
                                         phaseMaxHistoricalVolFraction1[0],
                                         phaseMinHistoricalVolFraction1[0],
                                         faceTrappedVolFrac1[0][0],
                                         faceCapPres1[0][0],
                                         dCapPres1_dfacePhaseVolFrac[0][0],
                                         scanBoundsMode1,
                                         phaseMode2PeakSlice1 );
                Pc1_max = faceCapPres1[0][0][0];
              }
            }
          }
          else
          {
            if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) { std::cout << "  [BOUNDS_DEBUG] Cell0: GENERIC (non-hysteresis) path" << std::endl; }
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
            if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) { std::cout << "  [BOUNDS_DEBUG] Cell0 GENERIC: Pc1_min=" << Pc1_min << ", Pc1_max=" << Pc1_max << std::endl; }
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
          else if constexpr (std::is_same_v< T2, TableCapillaryPressureHysteresis >) {
            auto preComputeMode1 = modes[1];
            capPresWrapper2.compute( phaseVolFrac2[0],
                                     phaseMaxHistoricalVolFraction2[0],
                                     phaseMinHistoricalVolFraction2[0],
                                     trappedVolFrac2[0][0],
                                     capPres2[0][0],
                                     dCapPres2_dPhaseVolFrac[0][0],
                                     preComputeMode1,
                                     phaseMode2PeakSlice2 );

            if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
              std::cout << "  [BOUNDS_DEBUG] Cell1: T2=TableCapillaryPressureHysteresis"
                        << ", preComputeMode=" << static_cast< int >(preComputeMode1)
                        << ", postComputeMode=" << static_cast< int >(modes[1])
                        << ", S=" << phaseVolFrac2[0][0]
                        << ", Pc=" << capPres2[0][0][0] << std::endl;
            }

            if( modes[1] == fields::cappres::ModeIndexType::DRAINAGE ||
                modes[1] == fields::cappres::ModeIndexType::IMBIBITION )
            {
              if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) { std::cout << "  [BOUNDS_DEBUG] Cell1: DRAINAGE/IMBIBITION path" << std::endl; }
              integer const tableIdx2 = static_cast< integer >(modes[1]);
              real64 Pc2_at_S1, dPc2_at_S1;
              real64 Pc2_at_S0, dPc2_at_S0;
              real64 const S_one_2 = 1.0;
              real64 const S_zero_2 = 0.0;
              capPresWrapper2.computeRawTablePc( tableIdx2, S_one_2, Pc2_at_S1, dPc2_at_S1 );
              capPresWrapper2.computeRawTablePc( tableIdx2, S_zero_2, Pc2_at_S0, dPc2_at_S0 );
              Pc2_min = LvArray::math::min( Pc2_at_S1, Pc2_at_S0 );
              Pc2_max = LvArray::math::max( Pc2_at_S1, Pc2_at_S0 );

              if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) { std::cout << "  [BOUNDS_DEBUG] Cell1: Pc2_min=" << Pc2_min
                                                                   << " (table@S=1:" << Pc2_at_S1 << ", table@S=0:" << Pc2_at_S0 << ")" << std::endl; }
            }
            else
            {
              if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) { std::cout << "  [BOUNDS_DEBUG] Cell1: SCANNING path, modes[1]=" << static_cast< int >(modes[1]) << std::endl; }
              real64 const S_min_hist2 = phaseMinHistoricalVolFraction2[0][0];
              real64 const S_max_hist2 = phaseMaxHistoricalVolFraction2[0][0];
              auto scanBoundsMode2 = modes[1];

              if( modes[1] == fields::cappres::ModeIndexType::DRAINAGE_TO_IMBIBITION ||
                  modes[1] == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE ||
                  modes[1] == fields::cappres::ModeIndexType::IMBIBITION_TO_DRAINAGE_FROM_SCANNING )
              {
                capPresWrapper2.computeScanningCurvePcRange(
                  S_min_hist2, S_max_hist2, modes[1], Pc2_min, Pc2_max );
              }
              else
              {
                facePhaseVolFrac2[0][1] = 0.0;
                facePhaseVolFrac2[0][0] = 1.0;
                scanBoundsMode2 = modes[1];
                capPresWrapper2.compute( facePhaseVolFrac2[0],
                                         phaseMaxHistoricalVolFraction2[0],
                                         phaseMinHistoricalVolFraction2[0],
                                         faceTrappedVolFrac2[0][0],
                                         faceCapPres2[0][0],
                                         dCapPres2_dfacePhaseVolFrac[0][0],
                                         scanBoundsMode2,
                                         phaseMode2PeakSlice2 );
                Pc2_min = faceCapPres2[0][0][0];

                facePhaseVolFrac2[0][1] = 1.0;
                facePhaseVolFrac2[0][0] = 0.0;
                scanBoundsMode2 = modes[1];
                capPresWrapper2.compute( facePhaseVolFrac2[0],
                                         phaseMaxHistoricalVolFraction2[0],
                                         phaseMinHistoricalVolFraction2[0],
                                         faceTrappedVolFrac2[0][0],
                                         faceCapPres2[0][0],
                                         dCapPres2_dfacePhaseVolFrac[0][0],
                                         scanBoundsMode2,
                                         phaseMode2PeakSlice2 );
                Pc2_max = faceCapPres2[0][0][0];
              }
            }
          }
          else
          {
            if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) { std::cout << "  [BOUNDS_DEBUG] Cell1: GENERIC (non-hysteresis) path" << std::endl; }
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
            if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) { std::cout << "  [BOUNDS_DEBUG] Cell1 GENERIC: Pc2_min=" << Pc2_min << ", Pc2_max=" << Pc2_max << std::endl; }
          }

          if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
            std::cout << "  [BOUNDS_DEBUG] FINAL: modes[0]=" << static_cast< int >(modes[0])
                      << ", modes[1]=" << static_cast< int >(modes[1])
                      << ", Pc1_min=" << Pc1_min << ", Pc1_max=" << Pc1_max
                      << ", Pc2_min=" << Pc2_min << ", Pc2_max=" << Pc2_max << std::endl;
          }

          // Use of the relative permeability kernel wrapper


          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > faceRelPerm1( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::relperm::LAYOUT_RELPERM_DS > dfacePhaseRelPerm1_dPhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > relPerm1( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::relperm::LAYOUT_RELPERM_DS > dPhaseRelPerm1_dPhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > faceRelPerm2( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::relperm::LAYOUT_RELPERM_DS > dfacePhaseRelPerm2_dPhaseVolFrac( 1, 1, 2, 2 );
          StackArray< real64, 3, 2, constitutive::relperm::LAYOUT_RELPERM > relPerm2( 1, 1, 2 );
          StackArray< real64, 4, 4, constitutive::relperm::LAYOUT_RELPERM_DS > dPhaseRelPerm2_dPhaseVolFrac( 1, 1, 2, 2 );


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

          real64 Pc_union_min = fmin( Pc1_min, Pc2_min );
          real64 Pc_union_max = fmax( Pc1_max, Pc2_max );

          real64 Pc_min_all = fmax( Pc1_min, Pc2_min );
          real64 Pc_max_all = fmin( Pc1_max, Pc2_max );
          bool rangesOverlap = Pc_min_all < Pc_max_all;

          real64 Pc_bracket_lo = Pc_union_min;
          real64 Pc_bracket_hi = Pc_union_max;
          real64 res_bracket_lo = 0.0;
          real64 res_bracket_hi = 0.0;
          bool bracket_lo_set = false;
          bool bracket_hi_set = false;

          constexpr bool ENABLE_WARM_START = true;
          real64 Pc_int;
          bool usedWarmStart = false;
          if( ENABLE_WARM_START && warmStartPc > 0.0 )
          {
            Pc_int = warmStartPc;
            usedWarmStart = true;
          }
          else
          {
            Pc_int = rangesOverlap ? ( Pc_min_all + Pc_max_all ) / 2.0
                                   : LvArray::math::max( Pc1, Pc2 );
          }
          Pc_int = fmin( Pc_union_max, fmax( Pc_int, Pc_union_min ) );

          real64 Pc_int_iterate = Pc_int;

          int iter = 0;
          int div = 0;

          int n_same_pc = 0;
          real64 prev_Pc_evaluated = -1.0e30;
          int stuck_probe_level = 0;

          constexpr bool ENABLE_BEST_SOLUTION_FALLBACK = false;
          constexpr real64 FALLBACK_TOL_FACTOR = 100.0;
          real64 best_Pc = 0.0;
          real64 best_absRes = 1.0e30;
          bool fallback_used = false;

          constexpr bool ENABLE_SWEEP_FALLBACK = true;
          constexpr real64 SWEEP_DS = 0.005;
          constexpr integer SWEEP_MIN_POINTS = 100;
          constexpr real64 SWEEP_ACCEPT_TOL_FACTOR = 100.0;
          bool sweep_active = false;
          real64 sweep_best_Pc = 0.0;
          real64 sweep_best_absRes = 1.0e30;
          int sweep_newton_iter_end = 0;

          real64 last_local_residual = 0.0;
          real64 last_local_jacobian = 0.0;
          real64 last_facePhaseVolFrac1_0 = 0.0;
          real64 last_facePhaseVolFrac2_0 = 0.0;
          real64 last_faceCapPres1 = 0.0;
          real64 last_faceCapPres2 = 0.0;
          integer last_iter = 0;

          while( iter < max_iter )
          {

            Pc_int_iterate = Pc_int;

            constexpr bool ENABLE_INITIAL_SAT_DEBUG = false;
            if constexpr (ENABLE_INITIAL_SAT_DEBUG) {
              if( iter == 0 )
              {
                std::cout << "[INITIAL_SAT_DEBUG] Local solver iter=" << iter << ", Initial S1=" << phaseVolFrac1[0][0]
                          << ", Initial S2=" << phaseVolFrac2[0][0] << std::endl;
                std::cout << "[INITIAL_SAT_DEBUG] Initial Pc1=" << capPres1[0][0][0]
                          << ", Initial Pc2=" << capPres2[0][0][0] << std::endl;
              }
            }

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

            {
              bool bracketReady = bracket_lo_set && bracket_hi_set;
              real64 clamp_lo = bracketReady ? Pc_bracket_lo : Pc_union_min;
              real64 clamp_hi = bracketReady ? Pc_bracket_hi : Pc_union_max;
              Pc_int = fmin( clamp_hi, fmax( Pc_int, clamp_lo ) );
            }

            faceCapPres1[0][0][0] = fmin( Pc1_max, Pc_int );
            faceCapPres2[0][0][0] = fmin( Pc2_max, Pc_int );
            bool const faceCapPres1_clamped = (faceCapPres1[0][0][0] < Pc_int - 1e-10);
            bool const faceCapPres2_clamped = (faceCapPres2[0][0][0] < Pc_int - 1e-10);

            JFunc1[0] = JFMultipliers[0];
            JFunc2[0] = JFMultipliers[1];

            constexpr bool ENABLE_TCH_COMPUTEINV_DEBUG = false;
            constexpr bool ENABLE_TCH_COMPUTE_DEBUG = false;

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
            else if constexpr (std::is_same_v< T3, TableCapillaryPressureHysteresis >) {
              if constexpr (ENABLE_TCH_COMPUTEINV_DEBUG) {
                if( iter < 5 || iter % 10 == 0 || iter >= max_iter - 5 )
                {
                  std::cout << "[TCH_COMPUTEINV_DEBUG] Cell1, iter=" << iter << std::endl;
                  std::cout << "  BEFORE computeInv: Pc=" << faceCapPres1[0][0][0]
                            << ", S=" << facePhaseVolFrac1[0][0]
                            << ", mode=" << modes[0] << std::endl;
                  std::cout << "  phaseMaxHistoricalVolFraction1[0][0]=" << phaseMaxHistoricalVolFraction1[0][0]
                            << ", phaseMinHistoricalVolFraction1[0][0]=" << phaseMinHistoricalVolFraction1[0][0] << std::endl;
                }
              }

              auto tempMode0 = modes[0];
              capPresWrapper1.computeInv( facePhaseVolFrac1[0],
                                          phaseMaxHistoricalVolFraction1[0],
                                          phaseMinHistoricalVolFraction1[0],
                                          phaseMode2PeakSlice1,
                                          faceTrappedVolFrac1[0][0],
                                          faceCapPres1[0][0],
                                          dfacePhaseVolFrac_dCapPres1[0][0],
                                          tempMode0 );
              real64 dS_dPc_after_computeInv_TCH = dfacePhaseVolFrac_dCapPres1[0][0][0][0];
              facePhaseVolFrac1[0][0] = fmin( 1.0, fmax( facePhaseVolFrac1[0][0], 0.0 ));
              facePhaseVolFrac1[0][1] = fmin( 1.0, fmax( facePhaseVolFrac1[0][1], 0.0 ));

              if constexpr (ENABLE_TCH_COMPUTEINV_DEBUG) {
                if( iter < 5 || iter % 10 == 0 || iter >= max_iter - 5 )
                {
                  std::cout << "  AFTER computeInv: Pc=" << faceCapPres1[0][0][0]
                            << ", S=" << facePhaseVolFrac1[0][0]
                            << ", dS_dPc=" << dS_dPc_after_computeInv_TCH << std::endl;
                }
              }

              real64 Pc_before_compute_TCH = faceCapPres1[0][0][0];
              if( modes[0] == fields::cappres::ModeIndexType::DRAINAGE ||
                  modes[0] == fields::cappres::ModeIndexType::IMBIBITION )
              {
                real64 rawPc = 0.0, rawDPcDS = 0.0;
                capPresWrapper1.computeRawTablePc( static_cast< integer >( modes[0] ),
                                                   facePhaseVolFrac1[0][0],
                                                   rawPc, rawDPcDS );
                faceCapPres1[0][0][0] = rawPc;
                dCapPres1_dfacePhaseVolFrac[0][0][0][0] = rawDPcDS;
              }
              else
              {
                auto tempComputeMode0 = modes[0];
                capPresWrapper1.compute( facePhaseVolFrac1[0],
                                         phaseMaxHistoricalVolFraction1[0],
                                         phaseMinHistoricalVolFraction1[0],
                                         faceTrappedVolFrac1[0][0],
                                         faceCapPres1[0][0],
                                         dCapPres1_dfacePhaseVolFrac[0][0],
                                         tempComputeMode0,
                                         phaseMode2PeakSlice1 );
              }
              real64 dPc_dS_after_compute_TCH = dCapPres1_dfacePhaseVolFrac[0][0][0][0];

              if constexpr (ENABLE_TCH_COMPUTE_DEBUG) {
                if( iter < 5 || iter % 10 == 0 || iter >= max_iter - 5 )
                {
                  std::cout << "[TCH_COMPUTE_DEBUG] Cell1, iter=" << iter << std::endl;
                  std::cout << "  BEFORE compute: Pc=" << Pc_before_compute_TCH
                            << ", S=" << facePhaseVolFrac1[0][0] << std::endl;
                  std::cout << "  AFTER compute: Pc=" << faceCapPres1[0][0][0]
                            << ", dPc_dS=" << dPc_dS_after_compute_TCH << std::endl;
                  std::cout << "  Pc change=" << (faceCapPres1[0][0][0] - Pc_before_compute_TCH) << std::endl;
                  real64 product = dS_dPc_after_computeInv_TCH * dPc_dS_after_compute_TCH;
                  std::cout << "  dS_dPc * dPc_dS = " << product
                            << " (should be ~1.0, error=" << std::abs( product - 1.0 ) << ")" << std::endl;
                }
              }

            }
            else
            {
              capPresWrapper1.computeInv( facePhaseVolFrac1[0],
                                          faceCapPres1[0][0],
                                          dfacePhaseVolFrac_dCapPres1[0][0] );
              real64 dS_dPc_after_computeInv_TCP = dfacePhaseVolFrac_dCapPres1[0][0][0][0];
              facePhaseVolFrac1[0][0] = fmin( 1.0, fmax( facePhaseVolFrac1[0][0], 0.0 ));
              facePhaseVolFrac1[0][1] = fmin( 1.0, fmax( facePhaseVolFrac1[0][1], 0.0 ));
              //get derivatives:
              real64 Pc_before_compute_TCP = faceCapPres1[0][0][0];
              real64 dPc_dS_before_TCP = dCapPres1_dfacePhaseVolFrac[0][0][0][0];
              capPresWrapper1.compute( facePhaseVolFrac1[0],
                                       faceCapPres1[0][0],
                                       dCapPres1_dfacePhaseVolFrac[0][0] );
              real64 dPc_dS_after_compute_TCP = dCapPres1_dfacePhaseVolFrac[0][0][0][0];
              constexpr bool ENABLE_COMPUTE_DEBUG = false;
              if constexpr (ENABLE_COMPUTE_DEBUG) {
                std::cout << "[COMPUTE_DEBUG_TCP] After computeInv(): dS_dPc=" << dS_dPc_after_computeInv_TCP
                          << ", S=" << facePhaseVolFrac1[0][0] << ", Pc=" << faceCapPres1[0][0][0] << std::endl;
                std::cout << "[COMPUTE_DEBUG_TCP] After compute(): Pc_before=" << Pc_before_compute_TCP
                          << ", Pc_after=" << faceCapPres1[0][0][0]
                          << ", dPc_dS_before=" << dPc_dS_before_TCP
                          << ", dPc_dS_after=" << dPc_dS_after_compute_TCP
                          << ", S=" << facePhaseVolFrac1[0][0] << std::endl;
                real64 product = dS_dPc_after_computeInv_TCP * dPc_dS_after_compute_TCP;
                std::cout << "[COMPUTE_DEBUG_TCP] dS_dPc * dPc_dS = " << product
                          << " (should be ~1.0, error=" << std::abs( product - 1.0 ) << ")" << std::endl;
                real64 Pc_diff = faceCapPres1[0][0][0] - Pc_before_compute_TCP;
                std::cout << "[COMPUTE_DEBUG_TCP] Pc change after compute() = " << Pc_diff
                          << " (should be ~0.0)" << std::endl;
              }
            }

            if( faceCapPres1_clamped )
            {
              dfacePhaseVolFrac_dCapPres1[0][0][0][0] = 0.0;
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
            else if constexpr (std::is_same_v< T4, TableCapillaryPressureHysteresis >) {
              if constexpr (ENABLE_TCH_COMPUTEINV_DEBUG) {
                if( iter < 5 || iter % 10 == 0 || iter >= max_iter - 5 )
                {
                  std::cout << "[TCH_COMPUTEINV_DEBUG] Cell2, iter=" << iter << std::endl;
                  std::cout << "  BEFORE computeInv: Pc=" << faceCapPres2[0][0][0]
                            << ", S=" << facePhaseVolFrac2[0][0]
                            << ", mode=" << modes[1] << std::endl;
                  std::cout << "  phaseMaxHistoricalVolFraction2[0][0]=" << phaseMaxHistoricalVolFraction2[0][0]
                            << ", phaseMinHistoricalVolFraction2[0][0]=" << phaseMinHistoricalVolFraction2[0][0] << std::endl;
                }
              }

              auto tempMode1 = modes[1];
              capPresWrapper2.computeInv( facePhaseVolFrac2[0],
                                          phaseMaxHistoricalVolFraction2[0],
                                          phaseMinHistoricalVolFraction2[0],
                                          phaseMode2PeakSlice2,
                                          faceTrappedVolFrac2[0][0],
                                          faceCapPres2[0][0],
                                          dfacePhaseVolFrac_dCapPres2[0][0],
                                          tempMode1 );
              facePhaseVolFrac2[0][0] = fmin( 1.0, fmax( facePhaseVolFrac2[0][0], 0.0 ));
              facePhaseVolFrac2[0][1] = fmin( 1.0, fmax( facePhaseVolFrac2[0][1], 0.0 ));
              real64 dS_dPc_after_computeInv_TCH2 = dfacePhaseVolFrac_dCapPres2[0][0][0][0];

              if constexpr (ENABLE_TCH_COMPUTEINV_DEBUG) {
                if( iter < 5 || iter % 10 == 0 || iter >= max_iter - 5 )
                {
                  std::cout << "  AFTER computeInv: Pc=" << faceCapPres2[0][0][0]
                            << ", S=" << facePhaseVolFrac2[0][0]
                            << ", dS_dPc=" << dS_dPc_after_computeInv_TCH2 << std::endl;
                }
              }

              real64 Pc_before_compute_TCH2 = faceCapPres2[0][0][0];

              if( modes[1] == fields::cappres::ModeIndexType::DRAINAGE ||
                  modes[1] == fields::cappres::ModeIndexType::IMBIBITION )
              {
                real64 rawPc2 = 0.0, rawDPcDS2 = 0.0;
                capPresWrapper2.computeRawTablePc( static_cast< integer >( modes[1] ),
                                                   facePhaseVolFrac2[0][0],
                                                   rawPc2, rawDPcDS2 );
                faceCapPres2[0][0][0] = rawPc2;
                dCapPres2_dfacePhaseVolFrac[0][0][0][0] = rawDPcDS2;
              }
              else
              {
                auto tempComputeMode1 = modes[1];
                capPresWrapper2.compute( facePhaseVolFrac2[0],
                                         phaseMaxHistoricalVolFraction2[0],
                                         phaseMinHistoricalVolFraction2[0],
                                         faceTrappedVolFrac2[0][0],
                                         faceCapPres2[0][0],
                                         dCapPres2_dfacePhaseVolFrac[0][0],
                                         tempComputeMode1,
                                         phaseMode2PeakSlice2 );
              }

              if constexpr (ENABLE_TCH_COMPUTE_DEBUG) {
                if( iter < 5 || iter % 10 == 0 || iter >= max_iter - 5 )
                {
                  std::cout << "[TCH_COMPUTE_DEBUG] Cell2, iter=" << iter << std::endl;
                  std::cout << "  BEFORE compute: Pc=" << Pc_before_compute_TCH2
                            << ", S=" << facePhaseVolFrac2[0][0] << std::endl;
                  std::cout << "  AFTER compute: Pc=" << faceCapPres2[0][0][0]
                            << ", dPc_dS=" << dCapPres2_dfacePhaseVolFrac[0][0][0] << std::endl;
                  std::cout << "  Pc change=" << (faceCapPres2[0][0][0] - Pc_before_compute_TCH2) << std::endl;
                  real64 dPc_dS_after_compute_TCH2 = dCapPres2_dfacePhaseVolFrac[0][0][0][0];
                  real64 product = dS_dPc_after_computeInv_TCH2 * dPc_dS_after_compute_TCH2;
                  std::cout << "  dS_dPc * dPc_dS = " << product
                            << " (should be ~1.0, error=" << std::abs( product - 1.0 ) << ")" << std::endl;
                }
              }

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

              real64 Pc_before_compute_TCP2 = faceCapPres2[0][0][0];
              capPresWrapper2.compute( facePhaseVolFrac2[0],
                                       faceCapPres2[0][0],
                                       dCapPres2_dfacePhaseVolFrac[0][0] );
              constexpr bool ENABLE_COMPUTE_DEBUG = false;
              if constexpr (ENABLE_COMPUTE_DEBUG) {
                std::cout << "[COMPUTE_DEBUG_TCP] After compute() cell2: Pc_before=" << Pc_before_compute_TCP2
                          << ", Pc_after=" << faceCapPres2[0][0][0]
                          << ", dPc_dS=" << dCapPres2_dfacePhaseVolFrac[0][0][0]
                          << ", S=" << facePhaseVolFrac2[0][0] << std::endl;
              }
            }

            if( faceCapPres2_clamped )
            {
              dfacePhaseVolFrac_dCapPres2[0][0][0][0] = 0.0;
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

            constexpr bool ENABLE_RELPERM_DEBUG = false;
            if constexpr (ENABLE_RELPERM_DEBUG) {
              if( iter == 0 )
              {
                std::cout << "[RELPERM_DEBUG] After relPerm1 compute: S1_face=" << facePhaseVolFrac1[0][0]
                          << ", faceRelPerm1[0][0][0]=" << faceRelPerm1[0][0][0]
                          << ", faceRelPerm1[0][0][1]=" << faceRelPerm1[0][0][1] << std::endl;
                std::cout << "[RELPERM_DEBUG] dfacePhaseRelPerm1_dPhaseVolFrac[0][0][0][0]="
                          << dfacePhaseRelPerm1_dPhaseVolFrac[0][0][0][0]
                          << ", dfacePhaseRelPerm1_dPhaseVolFrac[0][0][1][1]="
                          << dfacePhaseRelPerm1_dPhaseVolFrac[0][0][1][1] << std::endl;
              }
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

            if constexpr (ENABLE_RELPERM_DEBUG) {
              if( iter == 0 )
              {
                std::cout << "[RELPERM_DEBUG] After relPerm2 compute: S2_face=" << facePhaseVolFrac2[0][0]
                          << ", faceRelPerm2[0][0][0]=" << faceRelPerm2[0][0][0]
                          << ", faceRelPerm2[0][0][1]=" << faceRelPerm2[0][0][1] << std::endl;
                std::cout << "[RELPERM_DEBUG] dfacePhaseRelPerm2_dPhaseVolFrac[0][0][0][0]="
                          << dfacePhaseRelPerm2_dPhaseVolFrac[0][0][0][0]
                          << ", dfacePhaseRelPerm2_dPhaseVolFrac[0][0][1][1]="
                          << dfacePhaseRelPerm2_dPhaseVolFrac[0][0][1][1] << std::endl;
              }
            }


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

                    constexpr bool ENABLE_DCAPPRES_DS_DEBUG = false;
                    if constexpr (ENABLE_DCAPPRES_DS_DEBUG) {
                      if( iter == 0 && ix == 0 && ip == 0 && ke == 0 )
                      {
                        std::cout << "[DCAPPRES_DS_DEBUG] iter=" << iter << ", ix=" << ix << ", ip=" << ip << ", ke=" << ke << std::endl;
                        std::cout << "[DCAPPRES_DS_DEBUG] dCapPres_dS=" << dCapPres_dS << " (from cell " << (ke+1) << ")" << std::endl;
                        std::cout << "[DCAPPRES_DS_DEBUG] dCapPres1_dfacePhaseVolFrac[0][0][0][0]=" << dCapPres1_dfacePhaseVolFrac[0][0][0][0] << std::endl;
                        std::cout << "[DCAPPRES_DS_DEBUG] dCapPres2_dfacePhaseVolFrac[0][0][0][0]=" << dCapPres2_dfacePhaseVolFrac[0][0][0][0] << std::endl;
                      }
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

                  constexpr bool ENABLE_UPWIND_DEBUG = false;
                  if constexpr (ENABLE_UPWIND_DEBUG) {
                    if( iter == 0 && ip == 0 )
                    {
                      std::cout << "[UPWIND_DEBUG] iter=" << iter << ", ip=" << ip << ", ix=" << ix << ", k_up_0[ip]=" << k_up_0[ip] << ", k_up_1[ip]=" << k_up_1[ip] << std::endl;
                      std::cout << "[UPWIND_DEBUG] potGrad=" << potGrad << ", capGrad[ip]=" << capGrad[ip] << std::endl;
                      std::cout << "[UPWIND_DEBUG] dfacePhaseRelPerm1_dPhaseVolFrac[0][0][ip][ip]=" << dfacePhaseRelPerm1_dPhaseVolFrac[0][0][ip][ip] << std::endl;
                      std::cout << "[UPWIND_DEBUG] dPhaseRelPerm1_dPhaseVolFrac[0][0][ip][ip]=" << dPhaseRelPerm1_dPhaseVolFrac[0][0][ip][ip] << std::endl;
                    }
                  }

                  if( k_up_0[ip] == 1 && ix == 0 )
                  {
                    mobility[ip] = faceRelPerm1[0][0][ip] / viscosity[ip];
                    dMob_dP[ip][k_up_0[ip]] = mobility[ip] * (-dVisc_dP[ip][ix] / viscosity[ip]);

                    dMob_dS[ip][k_up_0[ip]] = sign[ip] * dfacePhaseRelPerm1_dPhaseVolFrac[0][0][ip][ip] / viscosity[ip];
                    if constexpr (ENABLE_UPWIND_DEBUG) {
                      if( iter == 0 && ip == 0 )
                      {
                        std::cout << "[UPWIND_DEBUG] Branch: k_up_0[ip]==1 && ix==0" << std::endl;
                        std::cout << "[UPWIND_DEBUG] dMob_dS[ip][k_up_0[ip]]=" << dMob_dS[ip][k_up_0[ip]] << ", ix=" << ix << std::endl;
                      }
                    }
                  }
                  else if( k_up_1[ip] == 0 && ix == 1 )
                  {
                    mobility[ip] = relPerm2[0][0][ip] / viscosity[ip];
                    dMob_dP[ip][0] = mobility[ip] * (-dVisc_dP[ip][ix] / viscosity[ip]);
                    dMob_dS[ip][0] = sign[ip] * dPhaseRelPerm2_dPhaseVolFrac[0][0][ip][ip] / viscosity[ip];
                    if constexpr (ENABLE_UPWIND_DEBUG) {
                      if( iter == 0 && ip == 0 )
                      {
                        std::cout << "[UPWIND_DEBUG] Branch: k_up_1[ip]==0 && ix==1" << std::endl;
                        std::cout << "[UPWIND_DEBUG] dMob_dS[ip][0]=" << dMob_dS[ip][0] << ", ix=" << ix << std::endl;
                      }
                    }
                  }
                  else if( k_up_1[ip] == 1 && ix == 1 )
                  {
                    mobility[ip] = faceRelPerm2[0][0][ip] / viscosity[ip];
                    dMob_dP[ip][1] = mobility[ip] * (-dVisc_dP[ip][ix] / viscosity[ip]);
                    dMob_dS[ip][1] = sign[ip] * dfacePhaseRelPerm2_dPhaseVolFrac[0][0][ip][ip] / viscosity[ip];
                    if constexpr (ENABLE_UPWIND_DEBUG) {
                      if( iter == 0 && ip == 0 )
                      {
                        std::cout << "[UPWIND_DEBUG] Branch: k_up_1[ip]==1 && ix==1" << std::endl;
                        std::cout << "[UPWIND_DEBUG] dMob_dS[ip][1]=" << dMob_dS[ip][1] << ", ix=" << ix << std::endl;
                      }
                    }
                  }
                  else
                  {
                    mobility[ip] = relPerm1[0][0][ip] / viscosity[ip];
                    dMob_dP[ip][ix] = mobility[ip] * (-dVisc_dP[ip][ix] / viscosity[ip]);
                    dMob_dS[ip][ix] = sign[ip] * dPhaseRelPerm1_dPhaseVolFrac[0][0][ip][ip] / viscosity[ip];
                    if constexpr (ENABLE_UPWIND_DEBUG) {
                      if( iter == 0 && ip == 0 )
                      {
                        std::cout << "[UPWIND_DEBUG] Branch: else (default)" << std::endl;
                        std::cout << "[UPWIND_DEBUG] dMob_dS[ip][ix]=" << dMob_dS[ip][ix] << ", ix=" << ix << std::endl;
                      }
                    }
                  }
                  real64 constexpr eps = 0.0;
                  total_mobility += mobility[ip] + eps;
                } // loop over phases

                constexpr bool ENABLE_MOBILITY_DEBUG = false;
                if constexpr (ENABLE_MOBILITY_DEBUG) {
                  if( iter == 0 && ix == 0 )
                  {
                    std::cout << "[MOBILITY_DEBUG] ix=" << ix << ", S1_face=" << facePhaseVolFrac1[0][0]
                              << ", S2_face=" << facePhaseVolFrac2[0][0] << std::endl;
                    std::cout << "[MOBILITY_DEBUG] mobility[0]=" << mobility[0] << ", mobility[1]=" << mobility[1]
                              << ", total_mobility=" << total_mobility << std::endl;
                    std::cout << "[MOBILITY_DEBUG] dMob_dS[0][0]=" << dMob_dS[0][0] << ", dMob_dS[0][1]=" << dMob_dS[0][1] << std::endl;
                    std::cout << "[MOBILITY_DEBUG] dMob_dS[1][0]=" << dMob_dS[1][0] << ", dMob_dS[1][1]=" << dMob_dS[1][1] << std::endl;
                    if( ix == 0 )
                    {
                      std::cout << "[MOBILITY_DEBUG] faceRelPerm1[0][0][0]=" << faceRelPerm1[0][0][0]
                                << ", faceRelPerm1[0][0][1]=" << faceRelPerm1[0][0][1] << std::endl;
                      std::cout << "[MOBILITY_DEBUG] dfacePhaseRelPerm1_dPhaseVolFrac[0][0][0][0]="
                                << dfacePhaseRelPerm1_dPhaseVolFrac[0][0][0][0]
                                << ", dfacePhaseRelPerm1_dPhaseVolFrac[0][0][1][1]="
                                << dfacePhaseRelPerm1_dPhaseVolFrac[0][0][1][1] << std::endl;
                    }
                  }
                }

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

                      constexpr bool ENABLE_DV_DEBUG = false;
                      if constexpr (ENABLE_DV_DEBUG) {
                        if( iter == 0 && ip == 0 && ke == 0 )
                        {
                          std::cout << "[DV_DEBUG] ix=" << ix << ", ip=" << ip << ", ke=" << ke
                                    << ", dV_dS=" << dV_dS << ", dV1_dS[ip][ke]=" << dV1_dS[ip][ke] << std::endl;
                          std::cout << "[DV_DEBUG] dMob_dS[0][ke]=" << dMob_dS[0][ke] << ", dMob_dS[1][ke]=" << dMob_dS[1][ke] << std::endl;
                          std::cout << "[DV_DEBUG] mobility[0]=" << mobility[0] << ", mobility[1]=" << mobility[1] << std::endl;
                          std::cout << "[DV_DEBUG] uT=" << uT << ", total_mobility=" << total_mobility << std::endl;
                        }
                      }
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
                      if( std::fabs( facePhaseVolFrac1[0][0] - 1.0 ) > 1e-8 )
                      {
                        dC1_dpc[ip][ke] =  dC1_dS[ip][ke] * dfacePhaseVolFrac_dCapPres1[0][0][0][0];
                      }
                      else
                      {
                        dC1_dpc[ip][ke] = dC_dS_term1 * dfacePhaseVolFrac_dCapPres1[0][0][0][0];
                      }

                      // GEOS_UNUSED_VAR( dC1_dpc[ip][ke] );
                    }
                    else
                    {
                      dC2_dS[ip][ke] = dC_dS_term1 + dC_dS_term2 * dC_dS_term3;
                      dhalfFlux2_dP[ip][ke] -= dC_dP;
                      dhalfFlux2_dS[ip][ke] -= dC2_dS[ip][ke];
                      if( std::fabs( facePhaseVolFrac2[0][0] - 1.0 ) > 1e-8 )
                      {
                        dC2_dpc[ip][ke] =  dC2_dS[ip][ke] * dfacePhaseVolFrac_dCapPres2[0][0][0][0];
                      }
                      else
                      {
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

            real64 constexpr eps2 = 1e-12;
            // real64 constexpr eps2 = 0.0;
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



            real64 const denom = dhalfFlux_dpc[0][0] - dhalfFlux_dpc[0][1];
            real64 const denomAbs = LvArray::math::abs( denom );
            local_jacobian = (denomAbs < eps2) ? (denom >= 0.0 ? eps2 : -eps2) : denom;
            local_residual = halfFluxVal[0][0] - halfFluxVal[0][1];



            last_local_residual = local_residual;
            last_local_jacobian = local_jacobian;
            last_facePhaseVolFrac1_0 = facePhaseVolFrac1[0][0];
            last_facePhaseVolFrac2_0 = facePhaseVolFrac2[0][0];
            last_faceCapPres1 = faceCapPres1[0][0][0];
            last_faceCapPres2 = faceCapPres2[0][0][0];
            last_iter = iter;

            constexpr bool ENABLE_NEWTON_ITER_DEBUG = false;
            if constexpr (ENABLE_NEWTON_ITER_DEBUG) {
              if( iter < 10 || iter % 10 == 0 || iter >= max_iter - 3 )
              {
                std::cout << "  [NI] it=" << iter
                          << " Pc=" << Pc_int
                          << " S1=" << facePhaseVolFrac1[0][0]
                          << " S2=" << facePhaseVolFrac2[0][0]
                          << " Pc1f=" << faceCapPres1[0][0][0]
                          << " Pc2f=" << faceCapPres2[0][0][0]
                          << " dSdPc1=" << dfacePhaseVolFrac_dCapPres1[0][0][0][0]
                          << " dSdPc2=" << dfacePhaseVolFrac_dCapPres2[0][0][0][0]
                          << " res=" << local_residual
                          << " jac=" << local_jacobian
                          << " dPc=" << (std::abs( local_jacobian ) > 1e-30 ? local_residual/local_jacobian : 0.0)
                          << std::endl;
              }
            }

            if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
              if( false )
              {
                std::cout << "[LOCAL_SOLVER_DEBUG_ITER] iter=" << iter << ", Pc_int=" << Pc_int
                          << ", Pc1_face=" << faceCapPres1[0][0][0] << ", Pc2_face=" << faceCapPres2[0][0][0] << std::endl;
                std::cout << "[LOCAL_SOLVER_DEBUG_ITER] S1_face=" << facePhaseVolFrac1[0][0]
                          << ", S2_face=" << facePhaseVolFrac2[0][0] << std::endl;
                std::cout << "[LOCAL_SOLVER_DEBUG_ITER] halfFluxVal[0][0]=" << halfFluxVal[0][0]
                          << ", halfFluxVal[0][1]=" << halfFluxVal[0][1] << std::endl;
                std::cout << "[LOCAL_SOLVER_DEBUG_ITER] dhalfFlux_dpc[0][0]=" << dhalfFlux_dpc[0][0]
                          << ", dhalfFlux_dpc[0][1]=" << dhalfFlux_dpc[0][1] << std::endl;
                std::cout << "[LOCAL_SOLVER_DEBUG_ITER] local_residual=" << local_residual
                          << ", local_jacobian=" << local_jacobian << std::endl;
              }
            }

            {
              real64 absRes = std::fabs( local_residual );
              if( absRes < best_absRes )
              {
                best_absRes = absRes;
                best_Pc = Pc_int;
              }
            }

            bool bracketEstablished = bracket_lo_set && bracket_hi_set;
            if( !bracketEstablished )
            {
              if( !bracket_lo_set && !bracket_hi_set )
              {
                Pc_bracket_lo = Pc_int;
                res_bracket_lo = local_residual;
                bracket_lo_set = true;
              }
              else
              {
                if( res_bracket_lo * local_residual < 0.0 )
                {
                  Pc_bracket_hi = Pc_int;
                  res_bracket_hi = local_residual;
                  bracket_hi_set = true;
                  if( Pc_bracket_lo > Pc_bracket_hi )
                  {
                    real64 tmpPc = Pc_bracket_lo;
                    real64 tmpRes = res_bracket_lo;
                    Pc_bracket_lo = Pc_bracket_hi;
                    res_bracket_lo = res_bracket_hi;
                    Pc_bracket_hi = tmpPc;
                    res_bracket_hi = tmpRes;
                  }
                }
                else
                {
                  Pc_bracket_lo = Pc_int;
                  res_bracket_lo = local_residual;
                }
              }
            }
            else
            {
              if( local_residual * res_bracket_lo <= 0.0 )
              {
                Pc_bracket_hi = Pc_int;
                res_bracket_hi = local_residual;
              }
              else
              {
                Pc_bracket_lo = Pc_int;
                res_bracket_lo = local_residual;
              }
            }
            bracketEstablished = bracket_lo_set && bracket_hi_set;

            constexpr bool ENABLE_BRACKET_DEBUG = false;
            if constexpr (ENABLE_BRACKET_DEBUG) {
              if( iter < 5 || iter % 10 == 0 )
              {
                std::cout << "  [BRACKET] iter=" << iter
                          << " established=" << bracketEstablished
                          << " lo=" << Pc_bracket_lo << " (res=" << res_bracket_lo << ")"
                          << " hi=" << Pc_bracket_hi << " (res=" << res_bracket_hi << ")"
                          << " width=" << (Pc_bracket_hi - Pc_bracket_lo) << std::endl;
              }
            }

            real64 const convergeTol = ( ENABLE_BEST_SOLUTION_FALLBACK && fallback_used ) ? ( FALLBACK_TOL_FACTOR * tol ) : tol;
            bool const eitherClamped = faceCapPres1_clamped || faceCapPres2_clamped;
            bool const rejectClampedConvergence = eitherClamped && !bracketEstablished && iter < max_iter - 5;
            if( std::fabs( local_residual ) < convergeTol && !rejectClampedConvergence )
            {
              converged = 1;
              // outFile << GEOS_FMT( "{:10.10e}", local_jacobian );
              // outFile << GEOS_FMT( ",{:10.10e}", local_residual );
              // outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[0][0] );
              // outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[0][1] );
              // outFile << GEOS_FMT( ",{:10.10e}", Pc_int_iterate );
              // outFile << GEOS_FMT( ",{:10.10e}", faceCapPres1[0][0][0] );
              // outFile << GEOS_FMT( ",{:10.10e}", faceCapPres2[0][0][0] );
              // outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[1][0] );
              // outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[1][1] );
              // outFile << GEOS_FMT( ",{:10.10e}", viscous[0][0] );
              // outFile << GEOS_FMT( ",{:10.10e}", viscous[1][0] );
              // outFile << GEOS_FMT( ",{:10.10e}", viscous[0][1] );
              // outFile << GEOS_FMT( ",{:10.10e}", viscous[1][1] );
              // outFile << GEOS_FMT( ",{:10.10e}", bouyancy[0][0] );
              // outFile << GEOS_FMT( ",{:10.10e}", bouyancy[1][0] );
              // outFile << GEOS_FMT( ",{:10.10e}", bouyancy[0][1] );
              // outFile << GEOS_FMT( ",{:10.10e}", bouyancy[1][1] );
              // outFile << GEOS_FMT( ",{:10.10e}", capillarity[0][0] );
              // outFile << GEOS_FMT( ",{:10.10e}", capillarity[1][0] );
              // outFile << GEOS_FMT( ",{:10.10e}", capillarity[0][1] );
              // outFile << GEOS_FMT( ",{:10.10e}", capillarity[1][1] );
              // outFile << GEOS_FMT( ",{:10.10e}", transHat[0] );
              // outFile << GEOS_FMT( ",{:10.10e}", transHat[1] );
              // outFile << GEOS_FMT( ",{:10.10e}", gravCoefHat[0] );
              // outFile << GEOS_FMT( ",{:10.10e}", gravCoef[0] );
              // outFile << GEOS_FMT( ",{:10.10e}", gravCoef[1] );
              // outFile << GEOS_FMT( ",{:10.10e}", uT );
              // outFile << GEOS_FMT( ",{:10.10e}", duT_dP[0] );
              // outFile << GEOS_FMT( ",{:10.10e}", duT_dS[0] );
              // outFile << GEOS_FMT( ",{:10.10e}", duT_dP[1] );
              // outFile << GEOS_FMT( ",{:10.10e}", duT_dS[1] );
              // outFile << std::endl;
              break; // Converged
            }

            {
              bool atLoBound = (std::fabs( Pc_int - Pc_union_min ) < 1e-10 * (std::fabs( Pc_union_min ) + 1.0));
              bool atHiBound = (std::fabs( Pc_int - Pc_union_max ) < 1e-10 * (std::fabs( Pc_union_max ) + 1.0));
              if( (atLoBound || atHiBound) && std::fabs( local_residual ) < 10.0 * tol && !rejectClampedConvergence )
              {
                converged = 1;
                break;
              }
            }

            if( iter >= max_iter - 1 )
            {
              if constexpr ( ENABLE_SWEEP_FALLBACK )
              {
                if( !sweep_active )
                {
                  sweep_active = true;
                  sweep_newton_iter_end = iter;
                  sweep_best_Pc = best_Pc;
                  sweep_best_absRes = best_absRes;
                  bracket_lo_set = false;
                  bracket_hi_set = false;
                  Pc_int = Pc_union_min;
                  max_iter += 500;
                  if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
                    std::cout << "[SWEEP] Activated at iter=" << iter
                              << ", sweeping [" << Pc_union_min << ", " << Pc_union_max << "]"
                              << ", dS=" << SWEEP_DS << std::endl;
                  }
                  iter++;
                  continue;  // re-enter loop to evaluate at Pc_union_min
                }
                else
                {
                  real64 const sweepAcceptTol = SWEEP_ACCEPT_TOL_FACTOR * tol;
                  if( sweep_best_absRes < sweepAcceptTol )
                  {
                    Pc_int = sweep_best_Pc;
                    converged = 1;
                    if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
                      std::cout << "[SWEEP] Accepted best: Pc=" << sweep_best_Pc
                                << ", |res|=" << sweep_best_absRes
                                << " < sweepTol=" << sweepAcceptTol << std::endl;
                    }
                    break;
                  }
                  else
                  {
                    if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
                      std::cout << "[SWEEP] DIVERGED: sweep_best_Pc=" << sweep_best_Pc
                                << ", sweep_best_absRes=" << sweep_best_absRes
                                << ", sweepAcceptTol=" << sweepAcceptTol << std::endl;
                    }
                    div = 1;
                    local_jacobian = 0.0;
                    break;
                  }
                }
              }
              else
              {
                if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
                  std::cout << "[LOCAL_SOLVER_DEBUG] DIVERGED: best_Pc=" << best_Pc
                            << ", best_absRes=" << best_absRes
                            << ", fallbackTol=" << FALLBACK_TOL_FACTOR * tol << std::endl;
                }
                div = 1;
                local_jacobian = 0.0;
                break;
              }
            }

            if( sweep_active && !bracketEstablished )
            {
              real64 const absRes_sweep = std::fabs( local_residual );
              if( absRes_sweep < sweep_best_absRes )
              {
                sweep_best_absRes = absRes_sweep;
                sweep_best_Pc = Pc_int;
              }

              real64 const absDSdPc1 = std::fabs( dfacePhaseVolFrac_dCapPres1[0][0][0][0] );
              real64 const absDSdPc2 = std::fabs( dfacePhaseVolFrac_dCapPres2[0][0][0][0] );
              real64 const max_absDSdPc = fmax( absDSdPc1, absDSdPc2 );
              real64 const dPc_max_cap = ( Pc_union_max - Pc_union_min ) / static_cast< real64 >( SWEEP_MIN_POINTS );
              real64 dPc_step;
              if( max_absDSdPc > 1e-20 )
              {
                dPc_step = SWEEP_DS / max_absDSdPc;
              }
              else
              {
                dPc_step = dPc_max_cap;
              }
              dPc_step = fmin( dPc_step, dPc_max_cap );

              real64 const Pc_next = Pc_int + dPc_step;

              if( Pc_next >= Pc_union_max )
              {
                Pc_int = sweep_best_Pc;
                max_iter = iter + 2;
                if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
                  std::cout << "[SWEEP] Completed full sweep, no sign change found."
                            << " sweep_best_Pc=" << sweep_best_Pc
                            << ", sweep_best_absRes=" << sweep_best_absRes << std::endl;
                }
              }
              else
              {
                Pc_int = Pc_next;
              }

              if constexpr (ENABLE_BRACKET_DEBUG) {
                if( (iter - sweep_newton_iter_end) % 20 == 0 || Pc_next >= Pc_union_max )
                {
                  std::cout << "  [SWEEP] Pc=" << Pc_int
                            << ", dPc_step=" << dPc_step
                            << ", best_Pc=" << sweep_best_Pc
                            << ", best_absRes=" << sweep_best_absRes << std::endl;
                }
              }
            }
            else if( sweep_active && bracketEstablished )
            {
              if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
                std::cout << "[SWEEP] Bracket found! Switching to bisection."
                          << " bracket=[" << Pc_bracket_lo << ", " << Pc_bracket_hi << "]" << std::endl;
              }
              Pc_int = ( Pc_bracket_lo + Pc_bracket_hi ) / 2.0;
            }
            else
            {
              real64 deltaPc = local_residual / local_jacobian;

              if( damping )
              {
                real64 max_dpc = fmax( fabs( dCapPres1_dfacePhaseVolFrac[0][0][0][0] ), fabs( dCapPres2_dfacePhaseVolFrac[0][0][0][0] ));
                real64 sign = std::copysign( 1.0, deltaPc );
                deltaPc = fmin( fabs( deltaPc ), max_dpc * 0.1 );
                deltaPc *= sign;
              }

              real64 Pc_newton = Pc_int - deltaPc;

              bool useBisection = false;
              if( bracketEstablished )
              {
                real64 brak_width = Pc_bracket_hi - Pc_bracket_lo;
                if( Pc_newton <= Pc_bracket_lo || Pc_newton >= Pc_bracket_hi ||
                    std::fabs( deltaPc ) > 0.75 * brak_width )
                {
                  useBisection = true;
                }
              }

              if( useBisection )
              {
                Pc_int = ( Pc_bracket_lo + Pc_bracket_hi ) / 2.0;
              }
              else
              {
                Pc_int = fmin( Pc_union_max, fmax( Pc_newton, Pc_union_min ) );
              }

              if( std::fabs( Pc_int - prev_Pc_evaluated ) < 1e-10 * std::fabs( Pc_int ) + 1e-30 )
              {
                n_same_pc++;
              }
              else
              {
                n_same_pc = 0;
              }
              prev_Pc_evaluated = Pc_int_iterate;

              if( n_same_pc >= 3 && !bracketEstablished )
              {
                real64 Pc_probe = 0.0;
                real64 const Pc_inner_max = fmin( Pc1_max, Pc2_max );
                if( stuck_probe_level == 0 )
                {
                  Pc_probe = ( Pc1 + Pc2 ) / 2.0;
                }
                else if( stuck_probe_level == 1 )
                {
                  Pc_probe = Pc_inner_max - 1.0;
                }
                else if( stuck_probe_level == 2 )
                {
                  Pc_probe = Pc_inner_max;
                }
                else
                {
                  if( std::fabs( Pc_int - Pc_union_min ) < std::fabs( Pc_int - Pc_union_max ) )
                    Pc_probe = Pc_union_max;
                  else
                    Pc_probe = Pc_union_min;
                }

                if( std::fabs( Pc_probe - Pc_int ) > 1e-10 * ( std::fabs( Pc_int ) + 1.0 ) )
                {
                  Pc_int = Pc_probe;
                }
                else
                {
                  stuck_probe_level++;
                  continue;
                }
                stuck_probe_level++;
                n_same_pc = 0;

                if constexpr (ENABLE_BRACKET_DEBUG) {
                  std::cout << "  [STUCK] Forced probe (level " << (stuck_probe_level - 1)
                            << ") to Pc_int=" << Pc_int << std::endl;
                }
              }

              if constexpr (ENABLE_BRACKET_DEBUG) {
                if( iter < 5 || iter % 10 == 0 )
                {
                  char const * stepName = useBisection ? "BISECTION" : "NEWTON";
                  std::cout << "  [STEP] " << stepName
                            << " -> Pc_int=" << Pc_int << std::endl;
                }
              }
            } // end step selection

            // truncate the updated capillary pressure (extended capillary pressure condition) for reporting/plotting:

            real64 faceCapPres1_plot = fmin( Pc1_max, fmax( Pc_int, Pc1_min ));
            real64 faceCapPres2_plot = fmin( Pc2_max, fmax( Pc_int, Pc2_min ));
            faceCapPres1_plot = fmin( Pc2_max, fmax( faceCapPres1_plot, Pc2_min ));
            faceCapPres2_plot = fmin( Pc1_max, fmax( faceCapPres2_plot, Pc1_min ));

            // Write data to the file
            // outFile << GEOS_FMT( "{:10.10e}", local_jacobian );
            // outFile << GEOS_FMT( ",{:10.10e}", local_residual );
            // outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[0][0] );
            // outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[0][1] );
            // outFile << GEOS_FMT( ",{:10.10e}", Pc_int_iterate );
            // outFile << GEOS_FMT( ",{:10.10e}", faceCapPres1[0][0][0] );
            // outFile << GEOS_FMT( ",{:10.10e}", faceCapPres2[0][0][0] );
            // outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[1][0] );
            // outFile << GEOS_FMT( ",{:10.10e}", halfFluxVal[1][1] );
            // outFile << GEOS_FMT( ",{:10.10e}", viscous[0][0] );
            // outFile << GEOS_FMT( ",{:10.10e}", viscous[1][0] );
            // outFile << GEOS_FMT( ",{:10.10e}", viscous[0][1] );
            // outFile << GEOS_FMT( ",{:10.10e}", viscous[1][1] );
            // outFile << GEOS_FMT( ",{:10.10e}", bouyancy[0][0] );
            // outFile << GEOS_FMT( ",{:10.10e}", bouyancy[1][0] );
            // outFile << GEOS_FMT( ",{:10.10e}", bouyancy[0][1] );
            // outFile << GEOS_FMT( ",{:10.10e}", bouyancy[1][1] );
            // outFile << GEOS_FMT( ",{:10.10e}", capillarity[0][0] );
            // outFile << GEOS_FMT( ",{:10.10e}", capillarity[1][0] );
            // outFile << GEOS_FMT( ",{:10.10e}", capillarity[0][1] );
            // outFile << GEOS_FMT( ",{:10.10e}", capillarity[1][1] );
            // outFile << GEOS_FMT( ",{:10.10e}", transHat[0] );
            // outFile << GEOS_FMT( ",{:10.10e}", transHat[1] );
            // outFile << GEOS_FMT( ",{:10.10e}", gravCoefHat[0] );
            // outFile << GEOS_FMT( ",{:10.10e}", gravCoef[0] );
            // outFile << GEOS_FMT( ",{:10.10e}", gravCoef[1] );
            // outFile << GEOS_FMT( ",{:10.10e}", uT );
            // outFile << GEOS_FMT( ",{:10.10e}", duT_dP[0] );
            // outFile << GEOS_FMT( ",{:10.10e}", duT_dS[0] );
            // outFile << GEOS_FMT( ",{:10.10e}", duT_dP[1] );
            // outFile << GEOS_FMT( ",{:10.10e}", duT_dS[1] );
            // outFile << std::endl;

            iter++;

          } // while loop

          if( converged )
          {
            warmStartPc = Pc_int;
          }

          if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
            const char * modeNames[] = {
              "DRAINAGE",
              "IMBIBITION",
              "DRAINAGE_TO_IMBIBITION",
              "IMBIBITION_TO_DRAINAGE",
              "IMBIBITION_TO_DRAINAGE_FROM_SCANNING"
            };
            integer mode0 = static_cast< integer >(modes[0]);
            integer mode1 = static_cast< integer >(modes[1]);
            const char * name0 = (mode0 >= 0 && mode0 <= 4) ? modeNames[mode0] : "UNKNOWN";
            const char * name1 = (mode1 >= 0 && mode1 <= 4) ? modeNames[mode1] : "UNKNOWN";

            std::cout << "[LOCAL_SOLVER_DEBUG] converged=" << converged << ", iter=" << last_iter << ", div=" << div << std::endl;
            std::cout << "[LOCAL_SOLVER_DEBUG] Modes: cell0=" << name0 << " (" << mode0 << "), cell1=" << name1 << " (" << mode1 << ")" << std::endl;
            std::cout << "[LOCAL_SOLVER_DEBUG] Capillary Pressures: Pc1=" << Pc1 << ", Pc2=" << Pc2 << ", Pc_int=" << Pc_int << std::endl;
            std::cout << "[LOCAL_SOLVER_DEBUG] Pc bounds: Pc1_min=" << Pc1_min << ", Pc1_max=" << Pc1_max << ", Pc2_min=" << Pc2_min << ", Pc2_max=" << Pc2_max << std::endl;
            std::cout << "[LOCAL_SOLVER_DEBUG] Initial saturations: S1=" << saturations[0] << ", S2=" << saturations[1] << std::endl;
            std::cout << "[LOCAL_SOLVER_DEBUG] Last face saturations: S1_face=" << last_facePhaseVolFrac1_0 << ", S2_face=" << last_facePhaseVolFrac2_0 << std::endl;
            std::cout << "[LOCAL_SOLVER_DEBUG] Last face capillary pressures: Pc1_face=" << last_faceCapPres1 << ", Pc2_face=" << last_faceCapPres2 << std::endl;
            std::cout << "[LOCAL_SOLVER_DEBUG] Last residual: local_residual=" << last_local_residual << ", local_jacobian=" << last_local_jacobian << std::endl;
            std::cout << "[LOCAL_SOLVER_DEBUG] Best solution: best_Pc=" << best_Pc << ", best_absRes=" << best_absRes
                      << ", fallbackTol=" << FALLBACK_TOL_FACTOR * tol << ", fallback_used=" << fallback_used << std::endl;
            std::cout << "[LOCAL_SOLVER_DEBUG] Warm-start: used=" << usedWarmStart << ", warmStartPc=" << warmStartPc << std::endl;
            if( sweep_active )
            {
              std::cout << "[LOCAL_SOLVER_DEBUG] Sweep: active, sweep_best_Pc=" << sweep_best_Pc
                        << ", sweep_best_absRes=" << sweep_best_absRes
                        << ", sweep_iter=" << (last_iter - sweep_newton_iter_end) << std::endl;
            }
          }

          if( converged )
          {

            real64 constexpr eps3 = 1e-12;
            real64 const denom_deriv = dhalfFlux_dpc[0][0] - dhalfFlux_dpc[0][1];
            real64 const denomAbs_deriv = LvArray::math::abs( denom_deriv );
            real64 const denom_protected = (denomAbs_deriv < eps3) ? (denom_deriv >= 0.0 ? eps3 : -eps3) : denom_deriv;

            real64 const dPc_int_dS1 =(-1.0) * (dhalfFlux1_dS[0][0] +  dhalfFlux_duT[0][0] * duT_dS[0] - dhalfFlux_duT[0][1] * duT_dS[0]) / denom_protected;
            real64 const dPc_int_dS2 =(-1.0) * (dhalfFlux_duT[0][0] * duT_dS[1] - dhalfFlux2_dS[0][0] - dhalfFlux_duT[0][1] * duT_dS[1]) / denom_protected;
            real64 const dPc_int_du =(-1.0) * (dhalfFlux_duT[0][0]  - dhalfFlux_duT[0][1]) / denom_protected;

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

            // std::cout << "dhalfFlux1_dS[0][0]=" << dhalfFlux1_dS[0][0] << ", duT_dS[0]=" << duT_dS[0] << std::endl;
            // std::cout << "dhalfFlux_dpc[0][1]=" << dhalfFlux_dpc[0][1] << ", dPc_int_dS1=" << dPc_int_dS1 << std::endl;
            // std::cout << "dhalfFlux1_dS[1][0]=" << dhalfFlux1_dS[1][0] << ", duT_dS[1]=" << duT_dS[1] << std::endl;
            // std::cout << "dhalfFlux2_dS[0][0]=" << dhalfFlux2_dS[0][0] << ", dPc_int_dS2=" << dPc_int_dS2 << std::endl;

            constexpr bool ENABLE_GLOBAL_DERIVATIVES_DEBUG = false;
            if constexpr (ENABLE_GLOBAL_DERIVATIVES_DEBUG) {
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG] fluxVal[0]=" << fluxVal[0] << ", fluxVal[1]=" << fluxVal[1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG] dFlux_dP[0][0]=" << dFlux_dP[0][0] << ", dFlux_dP[0][1]=" << dFlux_dP[0][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG] dFlux_dP[1][0]=" << dFlux_dP[1][0] << ", dFlux_dP[1][1]=" << dFlux_dP[1][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG] dFlux_dS[0][0]=" << dFlux_dS[0][0] << ", dFlux_dS[0][1]=" << dFlux_dS[0][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG] dFlux_dS[1][0]=" << dFlux_dS[1][0] << ", dFlux_dS[1][1]=" << dFlux_dS[1][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG] Intermediate values:" << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   dPc_int_dS1=" << dPc_int_dS1 << ", dPc_int_dS2=" << dPc_int_dS2 << ", dPc_int_du=" << dPc_int_du << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   halfFluxVal[0][0]=" << halfFluxVal[0][0] << ", halfFluxVal[0][1]=" << halfFluxVal[0][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   halfFluxVal[1][0]=" << halfFluxVal[1][0] << ", halfFluxVal[1][1]=" << halfFluxVal[1][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   dhalfFlux_dpc[0][0]=" << dhalfFlux_dpc[0][0] << ", dhalfFlux_dpc[0][1]=" << dhalfFlux_dpc[0][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   dhalfFlux_dpc[1][0]=" << dhalfFlux_dpc[1][0] << ", dhalfFlux_dpc[1][1]=" << dhalfFlux_dpc[1][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   dhalfFlux1_dS[0][0]=" << dhalfFlux1_dS[0][0] << ", dhalfFlux1_dS[0][1]=" << dhalfFlux1_dS[0][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   dhalfFlux1_dS[1][0]=" << dhalfFlux1_dS[1][0] << ", dhalfFlux1_dS[1][1]=" << dhalfFlux1_dS[1][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   dhalfFlux2_dS[0][0]=" << dhalfFlux2_dS[0][0] << ", dhalfFlux2_dS[0][1]=" << dhalfFlux2_dS[0][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   dhalfFlux2_dS[1][0]=" << dhalfFlux2_dS[1][0] << ", dhalfFlux2_dS[1][1]=" << dhalfFlux2_dS[1][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   dhalfFlux_duT[0][0]=" << dhalfFlux_duT[0][0] << ", dhalfFlux_duT[0][1]=" << dhalfFlux_duT[0][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   dhalfFlux_duT[1][0]=" << dhalfFlux_duT[1][0] << ", dhalfFlux_duT[1][1]=" << dhalfFlux_duT[1][1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   duT_dS[0]=" << duT_dS[0] << ", duT_dS[1]=" << duT_dS[1] << std::endl;
              std::cout << "[GLOBAL_DERIVATIVES_DEBUG]   density2[0]=" << density2[0] << ", density2[1]=" << density2[1] << std::endl;
            }

          }
          else
          {

            if constexpr (ENABLE_LOCAL_SOLVER_DEBUG) {
              std::cout << "**********************Diverged*******************" << std::endl;

              const char * modeNames[] = {
                "DRAINAGE",
                "IMBIBITION",
                "DRAINAGE_TO_IMBIBITION",
                "IMBIBITION_TO_DRAINAGE",
                "IMBIBITION_TO_DRAINAGE_FROM_SCANNING"
              };
              integer mode0 = static_cast< integer >(modes[0]);
              integer mode1 = static_cast< integer >(modes[1]);
              const char * name0 = (mode0 >= 0 && mode0 <= 4) ? modeNames[mode0] : "UNKNOWN";
              const char * name1 = (mode1 >= 0 && mode1 <= 4) ? modeNames[mode1] : "UNKNOWN";
              std::cout << "  [DEBUG] Modes: cell0=" << name0 << " (" << mode0 << "), cell1=" << name1 << " (" << mode1 << ")" << std::endl;
              std::cout << "  [DEBUG] Capillary Pressures: Pc1=" << Pc1 << ", Pc2=" << Pc2 << ", Pc_int=" << Pc_int << std::endl;
              std::cout << "  [DEBUG] Pc bounds: Pc1_min=" << Pc1_min << ", Pc1_max=" << Pc1_max << ", Pc2_min=" << Pc2_min << ", Pc2_max=" << Pc2_max << std::endl;
              std::cout << "  [DEBUG] Initial saturations: S1=" << saturations[0] << ", S2=" << saturations[1] << std::endl;
              std::cout << "  [DEBUG] Last iteration values: iter=" << last_iter << ", div=" << div << std::endl;
              std::cout << "  [DEBUG] Last face saturations: S1_face=" << last_facePhaseVolFrac1_0 << ", S2_face=" << last_facePhaseVolFrac2_0 << std::endl;
              std::cout << "  [DEBUG] Last face capillary pressures: Pc1_face=" << last_faceCapPres1 << ", Pc2_face=" << last_faceCapPres2 << std::endl;
              std::cout << "  [DEBUG] Last residual: local_residual=" << last_local_residual << ", local_jacobian=" << last_local_jacobian << std::endl;
            }

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


// Close the file after writing
          // outFile.close();

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
    StencilMaterialAccessors< CapillaryPressureBase,
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
                                       stdVector< std::array< localIndex, 2 > > const & interfacePairRegionIndices,
                                       arrayView1d< localIndex const > const & interfaceRegionLookup,
                                       arrayView1d< real64 > const & convergedPcInt,
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
    m_interfacePairRegionIndices( interfacePairRegionIndices ),
    m_interfaceRegionLookup( interfaceRegionLookup ),
    m_convergedPcInt( convergedPcInt )
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
    if( anyInterfaceConditionsQ )
    {
      connectorHasInterfaceConditionQ =
        ( iconn < m_interfaceRegionLookup.size() && m_interfaceRegionLookup[iconn] >= 0 );
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
        // bool notOnInterface = std::fabs( jFMultiplier[0][0] - jFMultiplier[0][1] ) < 1  && std::fabs( jFMultiplier[1][0] -
        // jFMultiplier[1][1] ) < 1;
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
          // trappedSats will be extracted from models below and passed to local_solver
          stdVector< real64 > trappedSats1 = {0.0, 0.0};
          stdVector< real64 > trappedSats2 = {0.0, 0.0};
          stdVector< fields::cappres::ModeIndexType > modes = {static_cast< fields::cappres::ModeIndexType >(0), static_cast< fields::cappres::ModeIndexType >(0)};
          stdVector< real64 > transHats = {transHat[0], transHat[1]};
          stdVector< real64 > dTransHats_dP = {dTransHat_dP[0], dTransHat_dP[1]};
          stdVector< real64 > gravCoefHats = {gravCoefHat[0], gravCoefHat[1]};
          stdVector< real64 > gravCoefs = {gravCoef2[0], gravCoef2[1]};
          stdVector< real64 > cellCenterDuTdS = {duT_dP[0], duT_dP[1], duT_dS[0], duT_dS[1]};
          stdVector< real64 > cellCenterDens = {density2[0], density2[1]};
          stdVector< real64 > cellCenterDens_dP = {dDens_dP2[0][0], dDens_dP2[0][1], dDens_dP2[1][0], dDens_dP2[1][1]};
          // std::vector< RelativePermeabilityBase * > relPerms = {std::get< 0 >( m_interfaceConstitutivePairs_temp ), std::get< 0 >(
          // m_interfaceConstitutivePairs_temp )};
          // std::vector< CapillaryPressureBase * > capPressures = {std::get< 1 >( m_interfaceConstitutivePairs_temp ), std::get< 1 >(
          // m_interfaceConstitutivePairs_temp )};
          // std::vector< TwoPhaseImmiscibleFluid * > fluids = {std::get< 2 >( m_interfaceConstitutivePairs_temp ), std::get< 2 >(
          // m_interfaceConstitutivePairs_temp )};

          // auto const & pairArray = m_interfaceConstitutivePairs[0];
          localIndex const surfaceRegionIndex = m_interfaceRegionLookup[iconn];
          auto const & pairArray = m_interfaceConstitutivePairs[surfaceRegionIndex];
          auto const & pairReg = m_interfacePairRegionIndices[surfaceRegionIndex];

// Determine canonical ordering based on global DOF numbers.
// Always put the cell with the SMALLER global DOF as local_solver side 0.
// This guarantees both ranks process the same interface connection identically,
// since local_solver is NOT symmetric w.r.t. swapping its two sides.
          globalIndex const dof0 = m_dofNumber[seri[0]][sesri[0]][sei[0]];
          globalIndex const dof1 = m_dofNumber[seri[1]][sesri[1]][sei[1]];

// p[i] maps canonical side i to stencil index (0 or 1)
          localIndex p[2] = {0, 1}; // default: stencil cell 0 has smaller DOF
          if( dof0 > dof1 )
          {
            // Stencil cell 1 has smaller DOF — swap to make it canonical side 0
            p[0] = 1;
            p[1] = 0;
          }

// Build local_solver inputs using mapped indices so that
// canonical side 0 always has the smaller global DOF

// uT is defined as flow from stencil cell 0 to stencil cell 1.
// If we swapped the cells, the direction reverses.
// duT derivatives must also be negated and cell-index-swapped.
          real64 const uT_sign = (p[0] == 0) ? 1.0 : -1.0;
          real64 const uT_ls = uT_sign * uT;

          stdVector< real64 > saturations_ls = {saturations[p[0]], saturations[p[1]]};
          stdVector< real64 > pressures_ls = {pressures[p[0]], pressures[p[1]]};
// transHat sign convention: transHat[0]=+T_0, transHat[1]=-T_1.
// When stencil is reversed, swap and negate to maintain: ls[0]=+T_canon0, ls[1]=-T_canon1
          stdVector< real64 > transHats_ls = {uT_sign * transHats[p[0]], uT_sign * transHats[p[1]]};
          stdVector< real64 > dTransHats_dP_ls = {uT_sign * dTransHats_dP[p[0]], uT_sign * dTransHats_dP[p[1]]};
          stdVector< real64 > gravCoefHats_ls = {gravCoefHats[p[0]], gravCoefHats[p[1]]};
          stdVector< real64 > gravCoefs_ls = {gravCoefs[p[0]], gravCoefs[p[1]]};
          stdVector< real64 > cellCenterDuTdS_ls = {uT_sign * cellCenterDuTdS[p[0]], uT_sign * cellCenterDuTdS[p[1]],
                                                    uT_sign * cellCenterDuTdS[2 + p[0]], uT_sign * cellCenterDuTdS[2 + p[1]]};
          stdVector< real64 > cellCenterDens_ls = {cellCenterDens[p[0]], cellCenterDens[p[1]]};
          stdVector< real64 > cellCenterDens_dP_ls = {cellCenterDens_dP[p[0]], cellCenterDens_dP[p[1]],
                                                      cellCenterDens_dP[2 + p[0]], cellCenterDens_dP[2 + p[1]]};

// Map sei indices for extracting per-element constitutive data.
// sei_ls[ke] is the element index for canonical side ke.
          localIndex const sei_ls[2] = {sei[p[0]], sei[p[1]]};

// Match constitutive models from the pair to the canonical ordering.
// pairArray[0] corresponds to pairReg[0], pairArray[1] to pairReg[1].
// Canonical side 0 is stencil cell p[0] with region seri[p[0]].
// We need to find which pair side matches each canonical side's region.
          localIndex cp[2] = {0, 1}; // cp[i] = which pair side goes with canonical side i
          if( pairReg[0] >= 0 && pairReg[1] >= 0 )
          {
            localIndex const canonReg0 = seri[p[0]]; // region of canonical side 0
            if( canonReg0 == pairReg[1] )
            {
              // canonical side 0's region matches pair side 1
              cp[0] = 1;
              cp[1] = 0;
            }
          }

          std::vector< constitutive::RelativePermeabilityBase * > relPerms = {
            std::get< 0 >( pairArray[cp[0]] ),
            std::get< 0 >( pairArray[cp[1]] )
          };

          std::vector< constitutive::CapillaryPressureBase * > capPressures = {
            std::get< 1 >( pairArray[cp[0]] ),
            std::get< 1 >( pairArray[cp[1]] )
          };

          std::vector< constitutive::TwoPhaseImmiscibleFluid * > fluids = {
            std::get< 2 >( pairArray[cp[0]] ),
            std::get< 2 >( pairArray[cp[1]] )
          };

          // Extract historical volume fractions, trapped saturations, and mode from TableCapillaryPressureHysteresis models
          // Initialize with sentinel values to indicate they haven't been set
          // Use -1.0 as a sentinel value (invalid for saturations which are 0-1)
          stdVector< real64 > phaseMaxHistoricalVolFraction1 = {-1.0, -1.0};
          stdVector< real64 > phaseMinHistoricalVolFraction1 = {-1.0, -1.0};
          stdVector< real64 > phaseMaxHistoricalVolFraction2 = {-1.0, -1.0};
          stdVector< real64 > phaseMinHistoricalVolFraction2 = {-1.0, -1.0};

          // Trapped saturations are available in all capillary pressure models (via base class),
          // but we only extract them when using hysteresis models since they're most meaningful there
          // Initialize with sentinel values to indicate they haven't been set
          stdVector< real64 > trappedSats1_extracted = {-1.0, -1.0};
          stdVector< real64 > trappedSats2_extracted = {-1.0, -1.0};

          for( integer ke = 0; ke < 2; ++ke )
          {
            constitutive::TableCapillaryPressureHysteresis * capPresHyst =
              dynamic_cast< constitutive::TableCapillaryPressureHysteresis * >(capPressures[ke]);

            if( capPresHyst != nullptr )
            {
              // Get the mode for this cell
              auto const & modeArray = capPresHyst->getField< fields::cappres::mode >().reference();
              if( sei_ls[ke] < static_cast< localIndex >(modeArray.size()) )
              {
                modes[ke] = static_cast< fields::cappres::ModeIndexType >(modeArray[sei_ls[ke]]);
              }

              // Get historical volume fractions for this cell
              auto const & maxHistArray = capPresHyst->getField< fields::cappres::phaseMaxHistoricalVolFraction >().reference();
              auto const & minHistArray = capPresHyst->getField< fields::cappres::phaseMinHistoricalVolFraction >().reference();

              if( sei_ls[ke] < static_cast< localIndex >(maxHistArray.size( 0 )) && m_numPhases > 0 )
              {
                for( integer ip = 0; ip < m_numPhases; ++ip )
                {
                  if( ke == 0 )
                  {
                    phaseMaxHistoricalVolFraction1[ip] = maxHistArray[sei_ls[ke]][ip];
                    phaseMinHistoricalVolFraction1[ip] = minHistArray[sei_ls[ke]][ip];
                  }
                  else
                  {
                    phaseMaxHistoricalVolFraction2[ip] = maxHistArray[sei_ls[ke]][ip];
                    phaseMinHistoricalVolFraction2[ip] = minHistArray[sei_ls[ke]][ip];
                  }
                }
              }

              // Get trapped volume fractions for this cell
              // Note: phaseTrappedVolFraction is available in all capillary pressure models (via base class)
              auto const & trappedArray = capPresHyst->getField< fields::cappres::phaseTrappedVolFraction >().reference();

              if( sei_ls[ke] < static_cast< localIndex >(trappedArray.size( 0 )) && m_numPhases > 0 )
              {
                for( integer ip = 0; ip < m_numPhases; ++ip )
                {
                  if( ke == 0 )
                  {
                    trappedSats1_extracted[ip] = trappedArray[sei_ls[ke]][0][ip]; // [element][subregion][phase]
                  }
                  else
                  {
                    trappedSats2_extracted[ip] = trappedArray[sei_ls[ke]][0][ip];
                  }
                }
              }
            }
            else
            {
              // For non-hysteresis models, we can still try to get trapped saturations if available
              // but they're typically not meaningful for non-hysteresis models
              // We'll extract them anyway in case they were set (e.g., for consistency)
              auto const & trappedArray = capPressures[ke]->getField< fields::cappres::phaseTrappedVolFraction >().reference();

              if( sei_ls[ke] < static_cast< localIndex >(trappedArray.size( 0 )) && m_numPhases > 0 )
              {
                for( integer ip = 0; ip < m_numPhases; ++ip )
                {
                  if( ke == 0 )
                  {
                    trappedSats1_extracted[ip] = trappedArray[sei_ls[ke]][0][ip];
                  }
                  else
                  {
                    trappedSats2_extracted[ip] = trappedArray[sei_ls[ke]][0][ip];
                  }
                }
              }
            }
          }

          stdVector< real64 > phi = {halfFluxVal[0][0], halfFluxVal[0][1]};
          stdVector< real64 > grad_phi_P = {0.0, 0.0, 0.0, 0.0};
          stdVector< real64 > grad_phi_S = {0.0, 0.0, 0.0, 0.0};


          // Use extracted trapped saturations if available, otherwise use defaults
          if( trappedSats1_extracted[0] >= 0.0 )
          {
            trappedSats1 = trappedSats1_extracted;
          }
          if( trappedSats2_extracted[0] >= 0.0 )
          {
            trappedSats2 = trappedSats2_extracted;
          }

          real64 warmStartPc = ( iconn < m_convergedPcInt.size() )
                               ? m_convergedPcInt[iconn] : -1.0;

          local_solver( uT_ls, saturations_ls, pressures_ls, JFMultipliers, trappedSats1, trappedSats2, modes, transHats_ls, dTransHats_dP_ls, gravCoefHats_ls, gravCoefs_ls,
                        cellCenterDuTdS_ls, cellCenterDens_ls, cellCenterDens_dP_ls, relPerms, capPressures, fluids, phi, grad_phi_P, grad_phi_S, converged,
                        phaseMaxHistoricalVolFraction1, phaseMinHistoricalVolFraction1, phaseMaxHistoricalVolFraction2, phaseMinHistoricalVolFraction2,
                        warmStartPc );

          if( iconn < m_convergedPcInt.size() )
          {
            m_convergedPcInt[iconn] = warmStartPc;
          }

          // Remap local_solver outputs back to stencil ordering for global assembly.
          // local_solver returns phi[0]=flux(phase0), phi[1]=flux(phase1) (phase-indexed).
          // grad_phi_P layout: [0]=dFlux0/dP_canon0, [1]=dFlux0/dP_canon1,
          //                    [2]=dFlux1/dP_canon0, [3]=dFlux1/dP_canon1.
          // We need to remap canonical cell indices back to stencil cell indices.
          // If canonical ordering differs from stencil ordering (p[0]==1), negate fluxes
          // (direction reversal) and swap cell derivative indices.

          if( p[0] == 0 )
          {
            // No swap needed - stencil ordering matches canonical ordering
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
          }
          else
          {
            // Stencil ordering is reversed relative to canonical.
            // Negate fluxes (direction reversal) and swap cell derivative indices.
            // Canon side 0 -> stencil cell 1, canon side 1 -> stencil cell 0
            fluxVal[0] = -phi[0];
            fluxVal[1] = -phi[1];
            // dFlux/dP: swap cell indices and negate
            dFlux_dP[0][0] = -grad_phi_P[1];     // dFlux0/dP_stencil0 = -dFlux0/dP_pair1
            dFlux_dP[0][1] = -grad_phi_P[0];     // dFlux0/dP_stencil1 = -dFlux0/dP_pair0
            dFlux_dP[1][0] = -grad_phi_P[3];     // dFlux1/dP_stencil0 = -dFlux1/dP_pair1
            dFlux_dP[1][1] = -grad_phi_P[2];     // dFlux1/dP_stencil1 = -dFlux1/dP_pair0
            dFlux_dS[0][0] = -grad_phi_S[1];
            dFlux_dS[0][1] = -grad_phi_S[0];
            dFlux_dS[1][0] = -grad_phi_S[3];
            dFlux_dS[1][1] = -grad_phi_S[2];
          }

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
  /// Region indices for each side of the interface pair (used for ordering consistency)
  stdVector< std::array< localIndex, 2 > > const m_interfacePairRegionIndices;
  arrayView1d< localIndex const > const m_interfaceRegionLookup;

  /// Per-instance warm-start capillary pressure array (one entry per connection)
  arrayView1d< real64 > const m_convergedPcInt;

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
   * @param[inout] convergedPcInt warm-start capillary pressure array (one entry per connection),
   *               populated by the caller using global face indices for MPI invariance
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
                   stdVector< std::array< localIndex, 2 > > const & interfacePairRegionIndices,
                   unordered_map< localIndex, localIndex > const & interfaceRegionByConnector,
                   ElementSubRegionBase const & subRegion,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   arrayView1d< real64 > const & convergedPcInt )
  {
    integer constexpr NUM_EQN = 2;
    integer constexpr NUM_DOF = 2;
    localIndex const domainSize = subRegion.size();
    localIndex const numConn = stencilWrapper.size();

    // Pre-compute the interface region lookup as a flat array for parallel-safe access.
    // This replaces unordered_map::find/at inside the parallel kernel, which is not thread-safe.
    array1d< localIndex > interfaceRegionLookup( numConn );
    interfaceRegionLookup.setValues< serialPolicy >( -1 );
    for( auto const & entry : interfaceRegionByConnector )
    {
      if( entry.first < numConn )
      {
        interfaceRegionLookup[entry.first] = entry.second;
      }
    }

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = FluxComputeInterfaceConditionKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
    typename kernelType::ImmiscibleMultiphaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::MultiphaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );

    kernelType kernel( numPhases, rankOffset, stencilWrapper, dofNumberAccessor,
                       flowAccessors, fluidAccessors, capPressureAccessors, permAccessors,
                       dt, localMatrix, localRhs, hasCapPressure, useTotalMassEquation,
                       checkPhasePresenceInGravity, interfaceFaceSetNames, interfaceConstitutivePairs,
                       interfacePairRegionIndices,
                       interfaceRegionLookup.toViewConst(), convergedPcInt, domainSize );
    kernelType::template launch< POLICY >( numConn, kernel );
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
