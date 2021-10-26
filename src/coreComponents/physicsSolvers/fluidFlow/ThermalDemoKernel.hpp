/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ThermalDemoKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALDEMOKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALDEMOKERNELS_HPP

#include "DemoKernel.hpp"

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "CompositionalMultiphaseUtilities.hpp"

namespace geosx
{

namespace ThermalDemoKernels
{

using namespace constitutive;

/**
 * @class ThermalAccumulationDemoKernel
 * @brief Define the interface for the thermal accumulation kernel
 */
template< localIndex NUM_COMP, localIndex NUM_DOF >
class ThermalAccumulationDemoKernel : public DemoKernels::AccumulationDemoKernel< NUM_COMP, NUM_DOF >
{
public:
  /// Alias for the base class;
  using Base = DemoKernels::AccumulationDemoKernel< NUM_COMP, NUM_DOF >;

  using Base::NC;
  using Base::NDOF;
  using Base::m_numPhases;
  using Base::m_rankOffset;
  using Base::m_dofNumber;
  using Base::m_elemGhostRank;
  using Base::m_volume;
  using Base::m_porosityOld;
  using Base::m_porosityNew;
  using Base::m_dPoro_dPres;
  using Base::m_dCompFrac_dCompDens;
  using Base::m_phaseVolFracOld;
  using Base::m_phaseVolFrac;
  using Base::m_dPhaseVolFrac_dPres;
  using Base::m_dPhaseVolFrac_dCompDens;
  using Base::m_phaseDensOld;
  using Base::m_phaseDens;
  using Base::m_dPhaseDens_dPres;
  using Base::m_dPhaseDens_dComp;
  using Base::m_phaseCompFracOld;
  using Base::m_phaseCompFrac;
  using Base::m_dPhaseCompFrac_dPres;
  using Base::m_dPhaseCompFrac_dComp;
  using Base::m_localMatrix;
  using Base::m_localRhs;

  /**
   * @brief Constructor
   */
  ThermalAccumulationDemoKernel( localIndex const numPhases,
                                 globalIndex const rankOffset,
                                 arrayView1d< globalIndex const > const & dofNumber,
                                 arrayView1d< integer const > const & elemGhostRank,
                                 arrayView1d< real64 const > const & volume,
                                 arrayView2d< real64 const > const & porosityOld,
                                 arrayView2d< real64 const > const & porosityNew,
                                 arrayView2d< real64 const > const & dPoro_dPres,
                                 arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
                                 arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFracOld,
                                 arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
                                 arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
                                 arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dTemp,
                                 arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dCompDens,
                                 arrayView2d< real64 const, compflow::USD_PHASE > const & phaseDensOld,
                                 arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
                                 arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
                                 arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dTemp,
                                 arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
                                 arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & phaseCompFracOld,
                                 arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & phaseCompFrac,
                                 arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & dPhaseCompFrac_dPres,
                                 arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & dPhaseCompFrac_dTemp,
                                 arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dPhaseCompFrac_dComp,
                                 arrayView2d< real64 const, compflow::USD_PHASE > const & phaseInternalEnergyOld,
                                 arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseInternalEnergy,
                                 arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseInternalEnergy_dPres,
                                 arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseInternalEnergy_dTemp,
                                 arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseInternalEnergy_dComp,
                                 arrayView1d< real64 const > const & rockInternalEnergyOld,
                                 arrayView2d< real64 const > const & rockInternalEnergy,
                                 arrayView2d< real64 const > const & dRockInternalEnergy_dTemp,
                                 arrayView2d< real64 const > const & rockDensity,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
    : Base( numPhases,
            rankOffset,
            dofNumber,
            elemGhostRank,
            volume,
            porosityOld,
            porosityNew,
            dPoro_dPres,
            dCompFrac_dCompDens,
            phaseVolFracOld,
            phaseVolFrac,
            dPhaseVolFrac_dPres,
            dPhaseVolFrac_dCompDens,
            phaseDensOld,
            phaseDens,
            dPhaseDens_dPres,
            dPhaseDens_dComp,
            phaseCompFracOld,
            phaseCompFrac,
            dPhaseCompFrac_dPres,
            dPhaseCompFrac_dComp,
            localMatrix,
            localRhs ),
    m_dPhaseVolFrac_dTemp( dPhaseVolFrac_dTemp ),
    m_dPhaseDens_dTemp( dPhaseDens_dTemp ),
    m_dPhaseCompFrac_dTemp( dPhaseCompFrac_dTemp ),
    m_phaseInternalEnergyOld( phaseInternalEnergyOld ),
    m_phaseInternalEnergy( phaseInternalEnergy ),
    m_dPhaseInternalEnergy_dPres( dPhaseInternalEnergy_dPres ),
    m_dPhaseInternalEnergy_dTemp( dPhaseInternalEnergy_dTemp ),
    m_dPhaseInternalEnergy_dComp( dPhaseInternalEnergy_dComp ),
    m_rockInternalEnergyOld( rockInternalEnergyOld ),
    m_rockInternalEnergy( rockInternalEnergy ),
    m_dRockInternalEnergy_dTemp( dRockInternalEnergy_dTemp ),
    m_rockDensity( rockDensity )
  {
    // or maybe now we want to get the views from here by passing the subRegion
  }

  // This seems useless, but if I remove it, it does not compile 
  struct StackVariables : public Base::StackVariables
  {
public:

    GEOSX_HOST_DEVICE
    StackVariables()
      : Base::StackVariables()
    {}

    using Base::StackVariables::localRow;
    using Base::StackVariables::dofIndices;
    using Base::StackVariables::localResidual;
    using Base::StackVariables::localJacobian;
  };


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void compute( localIndex const ei,
                StackVariables & stack ) const
  {
    // Here I need help: I don't want the [&], I would like [=] instead
    // But, if I do that I cannot compile... (read-only variable)
    // So this will have to be changed for the GPU
    // Maybe I can manually pass the local residual/jacobian "from the other way" with the poreVolume and phaseAmount
    Base::compute( ei, stack, [&] GEOSX_HOST_DEVICE ( localIndex const ip,
                                                      real64 const & poreVolumeNew,
                                                      real64 const &,
                                                      real64 const &,
                                                      real64 const & phaseAmountNew,
                                                      real64 const & phaseAmountOld,
                                                      real64 const & dPhaseAmount_dP,
                                                      real64 const (&dPhaseAmount_dC)[ NC ] )
    {
      // We are in the loop over phases, ip provides the current phase index.
      // We have to do two things:
      //   1- Assemble the derivatives of the component mass balance equations with respect to temperature
      //   2- Assemble the accumulation term of the energy equation

      // compute solid volumes (maybe move to stackVariables to avoid recomputing)
      real64 const solidVolumeNew = m_volume[ei] * ( 1.0 - m_porosityNew[ei][0] );
      real64 const solidVolumeOld = m_volume[ei] * ( 1.0 - m_porosityOld[ei][0] );
      real64 const dSolidVolume_dPres = -m_volume[ei] * m_dPoro_dPres[ei][0];

      /// rock internal energy (same, move) 
      real64 const rockInternalEnergyOld = m_rockInternalEnergyOld[ei];
      real64 const rockInternalEnergy = m_rockInternalEnergy[ei][0];
      real64 const dRockInternalEnergy_dTemp = m_dRockInternalEnergy_dTemp[ei][0];
      real64 const rockDensity = m_rockDensity[ei][0];

      // construct the slices (find a better way)
      arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];

      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > dPhaseVolFrac_dTemp = m_dPhaseVolFrac_dTemp[ei];

      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens = m_phaseDens[ei][0];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > dPhaseDens_dTemp = m_dPhaseDens_dTemp[ei][0];

      arraySlice2d< real64 const, multifluid::USD_PHASE_COMP-2 > phaseCompFrac = m_phaseCompFrac[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_COMP-2 > dPhaseCompFrac_dTemp = m_dPhaseCompFrac_dTemp[ei][0];

      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseInternalEnergyOld = m_phaseInternalEnergyOld[ei];
      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > phaseInternalEnergy = m_phaseInternalEnergy[ei][0];
      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > dPhaseInternalEnergy_dPres = m_dPhaseInternalEnergy_dPres[ei][0];
      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > dPhaseInternalEnergy_dTemp = m_dPhaseInternalEnergy_dTemp[ei][0];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > dPhaseInternalEnergy_dComp = m_dPhaseInternalEnergy_dComp[ei][0];


      // Step 1: assemble the derivatives of the component mass balance equations with respect to temperature

      real64 const dPoreVolume_dTemp = 0.0;
      real64 const dPhaseAmount_dT = dPoreVolume_dTemp * phaseVolFrac[ip] * phaseDens[ip]
                                     + poreVolumeNew * (dPhaseVolFrac_dTemp[ip] * phaseDens[ip]
                                                        + phaseVolFrac[ip] * dPhaseDens_dTemp[ip]);

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        stack.localJacobian[ic][NC+1] += dPhaseAmount_dT * phaseCompFrac[ip][ic]
                                         + phaseAmountNew * dPhaseCompFrac_dTemp[ip][ic];
      }


      // Step 2: assemble the accumulation term of the energy equation

      real64 dPhaseInternalEnergy_dC[NC]{};

      real64 const phaseEnergyNew = phaseAmountNew * phaseInternalEnergy[ip];
      real64 const phaseEnergyOld = phaseAmountOld * phaseInternalEnergyOld[ip];

      real64 const solidEnergyNew = solidVolumeNew * rockInternalEnergy * rockDensity;
      real64 const solidEnergyOld = solidVolumeOld * rockInternalEnergyOld * rockDensity;

      // local accumulation
      stack.localResidual[NC] = phaseEnergyNew - phaseEnergyOld
                                + solidEnergyNew - solidEnergyOld;

      real64 const dPhaseEnergy_dP = dPhaseAmount_dP * phaseInternalEnergy[ip] + phaseAmountNew * dPhaseInternalEnergy_dPres[ip];
      real64 const dPhaseEnergy_dT = dPhaseAmount_dT * phaseInternalEnergy[ip] + phaseAmountNew * dPhaseInternalEnergy_dTemp[ip];

      real64 const dSolidInternalEnergy_dP = dSolidVolume_dPres * rockInternalEnergy * rockDensity;
      real64 const dSolidInternalEnergy_dT = solidVolumeNew * dRockInternalEnergy_dTemp;

      // derivatives w.r.t. pressure and temperature
      stack.localJacobian[NC][0]      += dPhaseEnergy_dP + dSolidInternalEnergy_dP;
      stack.localJacobian[NC][NDOF-1] += dPhaseEnergy_dT + dSolidInternalEnergy_dT;

      // derivatives w.r.t. component densities
      applyChainRule( NC, dCompFrac_dCompDens, dPhaseInternalEnergy_dComp[ip], dPhaseInternalEnergy_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        real64 const dPhaseEnergy_dC = phaseInternalEnergy[ip] * dPhaseAmount_dC[jc]
                                       + dPhaseInternalEnergy_dC[jc] * phaseAmountNew;

        stack.localJacobian[NC][jc + 1] += dPhaseEnergy_dC;
      }
    } );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // Step 1: assemble the component mass balance equations
    Base::complete( ei, stack );

    // Step 2: assemble the energy equation
    m_localRhs[stack.localRow + NDOF-1] += stack.localResidual[NDOF-2];
    m_localMatrix.template addToRow< serialAtomic >( stack.localRow + NDOF-1,
                                                     stack.dofIndices,
                                                     stack.localJacobian[NDOF-2],
                                                     NDOF );

  }

private:

  /// Views on derivatives wrt to temperature for phase volume fraction, density, and phase comp fraction
  arrayView2d< real64 const, compflow::USD_PHASE > const & m_dPhaseVolFrac_dTemp;
  arrayView3d< real64 const, multifluid::USD_PHASE > const & m_dPhaseDens_dTemp;
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & m_dPhaseCompFrac_dTemp;

  /// Views on phase internal energy
  arrayView2d< real64 const, compflow::USD_PHASE > const & m_phaseInternalEnergyOld;
  arrayView3d< real64 const, multifluid::USD_PHASE > const & m_phaseInternalEnergy;
  arrayView3d< real64 const, multifluid::USD_PHASE > const & m_dPhaseInternalEnergy_dPres;
  arrayView3d< real64 const, multifluid::USD_PHASE > const & m_dPhaseInternalEnergy_dTemp;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & m_dPhaseInternalEnergy_dComp;

  /// Views on rock internal energy
  arrayView1d< real64 const > const & m_rockInternalEnergyOld;
  arrayView2d< real64 const > const & m_rockInternalEnergy;
  arrayView2d< real64 const > const & m_dRockInternalEnergy_dTemp;

  /// View on rock density
  arrayView2d< real64 const > const & m_rockDensity;

};

} // namespace ThermalDemoKernels

} // namespace geosx

#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_THERMALDEMOKERNELS_HPP
