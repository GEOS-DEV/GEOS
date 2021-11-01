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
 * @file CompositionalMultiphaseBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geosx
{

namespace CompositionalMultiphaseBaseKernels
{

using namespace constitutive;

static constexpr real64 minDensForDivision = 1e-10;

/******************************** ComponentFractionKernel ********************************/

/**
 * @brief Functions to compute component fractions from global component densities (mass or molar)
 */
struct ComponentFractionKernel
{
  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
  compute( arraySlice1d< real64 const, compflow::USD_COMP - 1 > compDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > dCompDens,
           arraySlice1d< real64, compflow::USD_COMP - 1 > compFrac,
           arraySlice2d< real64, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens );

  template< localIndex NC >
  static void
  launch( localIndex const size,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView2d< real64, compflow::USD_COMP > const & compFrac,
          arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens );

  template< localIndex NC >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView2d< real64, compflow::USD_COMP > const & compFrac,
          arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens );
};

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @brief Functions to compute phase volume fractions (saturations) and derivatives
 */
struct PhaseVolumeFractionKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  compute( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & compDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & dCompDens,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dComp,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseFrac,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseFrac_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseFrac_dComp,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseFrac_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseFrac_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp );
};


/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          real64 const temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp, compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          real64 const temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k] + dPres[k], temp, compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          real64 const temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k] + dPres[k], temp, compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          real64 const temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp, compFrac[k] );
      }
    } );
  }
};

/******************************** RelativePermeabilityUpdateKernel ********************************/

struct RelativePermeabilityUpdateKernel
{
  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( localIndex const size,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

/******************************** CapillaryPressureUpdateKernel ********************************/

struct CapillaryPressureUpdateKernel
{
  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  launch( localIndex const size,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of accumulation and volume balance
 */
template< localIndex NUM_COMP, localIndex NUM_DOF >
class ElementBasedAssemblyKernel
{
public:

  /// Compile time value for the number of components
  static constexpr localIndex NC = NUM_COMP;

  /// Compute time value for the number of degrees of freedom
  static constexpr localIndex NDOF = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr localIndex NEQ = NUM_DOF;

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
   */
  ElementBasedAssemblyKernel( localIndex const numPhases,
                              globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              MultiFluidBase const & fluid,
                              CoupledSolidBase const & solid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    : m_numPhases( numPhases ),
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_volume( subRegion.getElementVolume() ),
    m_porosityOld( solid.getOldPorosity() ),
    m_porosityNew( solid.getPorosity() ),
    m_dPoro_dPres( solid.getDporosity_dPressure() ),
    m_dCompFrac_dCompDens( subRegion.getReference< array3d< real64, compflow::LAYOUT_COMP_DC > >( CompositionalMultiphaseBase::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() ) ),
    m_phaseVolFracOld( subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( CompositionalMultiphaseBase::viewKeyStruct::phaseVolumeFractionOldString() ) ),
    m_phaseVolFrac( subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( CompositionalMultiphaseBase::viewKeyStruct::phaseVolumeFractionString() ) ),
    m_dPhaseVolFrac_dPres( subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( CompositionalMultiphaseBase::viewKeyStruct::dPhaseVolumeFraction_dPressureString() ) ),
    m_dPhaseVolFrac_dCompDens( subRegion.getReference< array3d< real64, compflow::LAYOUT_PHASE_DC > >( CompositionalMultiphaseBase::viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() ) ),
    m_phaseDensOld( subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( CompositionalMultiphaseBase::viewKeyStruct::phaseDensityOldString() ) ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens_dPres( fluid.dPhaseDensity_dPressure() ),
    m_dPhaseDens_dComp( fluid.dPhaseDensity_dGlobalCompFraction() ),
    m_phaseCompFracOld( subRegion.getReference< array3d< real64, compflow::LAYOUT_PHASE_COMP > >( CompositionalMultiphaseBase::viewKeyStruct::phaseComponentFractionOldString() ) ),
    m_phaseCompFrac( fluid.phaseCompFraction() ),
    m_dPhaseCompFrac_dPres( fluid.dPhaseCompFraction_dPressure() ),
    m_dPhaseCompFrac_dComp( fluid.dPhaseCompFraction_dGlobalCompFraction() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    GEOSX_HOST_DEVICE
    StackVariables() {}

    // Pore volume information (used by both accumulation and volume balance)

    /// Pore volume at time n+1
    real64 poreVolumeNew = 0.0;

    /// Pore volume at the previous converged time step
    real64 poreVolumeOld = 0.0;

    /// Derivative of pore volume with respect to pressure
    real64 dPoreVolume_dPres = 0.0;

    // Residual information

    /// Index of the local row corresponding to this element
    localIndex localRow = -1;

    /// Indices of the matrix rows/columns corresponding to the dofs in this element
    globalIndex dofIndices[NDOF]{};

    /// C-array storage for the element local residual vector (all equations except volume balance)
    real64 localResidual[NEQ]{};

    /// C-array storage for the element local Jacobian matrix (all equations except volume balance, all dogs)
    real64 localJacobian[NEQ][NDOF]{};

  };

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 getElemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    // initialize the pore volume
    stack.poreVolumeNew = m_volume[ei] * m_porosityNew[ei][0];
    stack.poreVolumeOld = m_volume[ei] * m_porosityOld[ei][0];
    stack.dPoreVolume_dPres = m_volume[ei] * m_dPoro_dPres[ei][0];

    // set row index and degrees of freedom indices for this element
    stack.localRow = m_dofNumber[ei] - m_rankOffset;
    for( integer idof = 0; idof < NDOF; ++idof )
    {
      stack.dofIndices[idof] = m_dofNumber[ei] + idof;
    }
  }

  /**
   * @brief Internal struct to provide no-op defaults used in the inclusion
   *   of lambda functions into kernel component functions.
   * @struct NoOpFunctors
   */
  struct NoOpFunctors
  {

    /**
     * @brief operator() no-op used for additional terms in the residual / Jacobian in the accumulation term
     * @tparam NC the number of components
     * @param[in] ip the phase index
     * @param[in] phaseAmountNew the volume of phase ip at time n+1
     * @param[in] phaseAmountOld the volume of phase ip at the previous converged time step
     * @param[in] dPhaseAmount_dPres the derivative of the volume of phase ip with respect to pressure
     * @param[in] dPhaseAmount_dCompDens the derivatives of the volume of phase ip with respect to component densities
     */
    template< localIndex NC >
    GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE constexpr
    void operator() ( localIndex const ip,
                      real64 const & phaseAmountNew,
                      real64 const & phaseAmountOld,
                      real64 const & dPhaseAmount_dPres,
                      real64 const (&dPhaseAmount_dCompDens )[ NC ] )
    {
      GEOSX_UNUSED_VAR( ip, phaseAmountNew, phaseAmountOld, dPhaseAmount_dPres, dPhaseAmount_dCompDens );
    }

    /**
     * @brief operator() no-op used for additional terms in the residual / Jacobian in the volume balance term
     * @param[in] oneMinusPhaseVolumeFracSum one minus the sum of the phase volume fractions in the element
     */
    GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE constexpr
    void operator() ( real64 const & oneMinusPhaseVolFracSum )
    {
      GEOSX_UNUSED_VAR( oneMinusPhaseVolFracSum );
    }

  };

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam LAMBDA the type of the functor that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] lambda the functor used to customize the kernel
   */
  template< typename LAMBDA = NoOpFunctors >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack,
                            LAMBDA && lambda = NoOpFunctors{} ) const
  {
    // construct the slices for variables accessed multiple times
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFracOld = m_phaseVolFracOld[ei];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[ei];

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseDensOld = m_phaseDensOld[ei];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens = m_phaseDens[ei][0];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > dPhaseDens_dPres = m_dPhaseDens_dPres[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens_dComp = m_dPhaseDens_dComp[ei][0];

    arraySlice2d< real64 const, compflow::USD_PHASE_COMP-1 > phaseCompFracOld = m_phaseCompFracOld[ei];
    arraySlice2d< real64 const, multifluid::USD_PHASE_COMP-2 > phaseCompFrac = m_phaseCompFrac[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_COMP-2 > dPhaseCompFrac_dPres = m_dPhaseCompFrac_dPres[ei][0];
    arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC-2 > dPhaseCompFrac_dComp = m_dPhaseCompFrac_dComp[ei][0];

    // temporary work arrays
    real64 dPhaseAmount_dC[NC]{};
    real64 dPhaseCompFrac_dC[NC]{};

    // sum contributions to component accumulation from each phase
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      real64 const phaseAmountNew = stack.poreVolumeNew * phaseVolFrac[ip] * phaseDens[ip];
      real64 const phaseAmountOld = stack.poreVolumeOld * phaseVolFracOld[ip] * phaseDensOld[ip];

      real64 const dPhaseAmount_dP = stack.dPoreVolume_dPres * phaseVolFrac[ip] * phaseDens[ip]
                                     + stack.poreVolumeNew * (dPhaseVolFrac_dPres[ip] * phaseDens[ip]
                                                              + phaseVolFrac[ip] * dPhaseDens_dPres[ip]);

      // assemble density dependence
      applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dPhaseAmount_dC );
      for( integer jc = 0; jc < NC; ++jc )
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                              + phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
        dPhaseAmount_dC[jc] *= stack.poreVolumeNew;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for( integer ic = 0; ic < NC; ++ic )
      {
        real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ip][ic];
        real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOld[ip][ic];

        real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic]
                                           + phaseAmountNew * dPhaseCompFrac_dPres[ip][ic];

        stack.localResidual[ic] += phaseCompAmountNew - phaseCompAmountOld;
        stack.localJacobian[ic][0] += dPhaseCompAmount_dP;

        // jc - index of component w.r.t. whose compositional var the derivative is being taken
        // (i.e. col number in local matrix)

        // assemble phase composition dependence
        applyChainRule( NC, dCompFrac_dCompDens, dPhaseCompFrac_dComp[ip][ic], dPhaseCompFrac_dC );
        for( integer jc = 0; jc < NC; ++jc )
        {
          real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                             + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
          stack.localJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
        }
      }

      // call the lambda in the phase loop to allow the reuse of the phase amounts and their derivatives
      // possible use: assemble the derivatives wrt temperature, and the accumulation term of the energy equation for this phase
      lambda( ip, phaseAmountNew, phaseAmountOld, dPhaseAmount_dP, dPhaseAmount_dC );

    }
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam LAMBDA the type of the functor that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] lambda the functor used to customize the kernel
   */
  template< typename LAMBDA = NoOpFunctors >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeVolumeBalance( localIndex const ei,
                             StackVariables & stack,
                             LAMBDA && lambda = NoOpFunctors{} ) const
  {
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[ei];

    real64 oneMinusPhaseVolFracSum = 1.0;

    // sum contributions to component accumulation from each phase
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      oneMinusPhaseVolFracSum -= phaseVolFrac[ip];
      stack.localJacobian[NC][0] -= dPhaseVolFrac_dPres[ip];

      for( integer jc = 0; jc < NC; ++jc )
      {
        stack.localJacobian[NC][jc+1] -= dPhaseVolFrac_dCompDens[ip][jc];
      }
    }

    // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
    stack.localResidual[NC] = stack.poreVolumeNew * oneMinusPhaseVolFracSum;
    for( integer idof = 0; idof < NC+1; ++idof )
    {
      stack.localJacobian[NC][idof] *= stack.poreVolumeNew;
    }
    stack.localJacobian[NC][0] += stack.dPoreVolume_dPres * oneMinusPhaseVolFracSum;

    // call the lambda in the phase loop to allow the reuse of the phase amounts and their derivatives
    // possible use: assemble the derivatives wrt temperature, and use oneMinusPhaseVolFracSum if poreVolumeNew depends on temperature
    lambda( oneMinusPhaseVolFracSum );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void complete( localIndex const GEOSX_UNUSED_PARAM( ei ),
                 StackVariables & stack ) const
  {
    using namespace CompositionalMultiphaseUtilities;

    // apply equation/variable change transformation to the component mass balance equations
    real64 work[NDOF]{};
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NDOF, stack.localJacobian, work );
    shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, stack.localResidual );

    // add contribution to residual and jacobian into:
    // - the component mass balance equations
    // - the volume balance equations
    for( integer i = 0; i < NC+1; ++i )
    {
      m_localRhs[stack.localRow + i] += stack.localResidual[i];
      m_localMatrix.addToRow< serialAtomic >( stack.localRow + i,
                                              stack.dofIndices,
                                              stack.localJacobian[i],
                                              NDOF );
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
    GEOSX_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.getElemGhostRank( ei ) >= 0 )
      {
        return;
      }

      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, stack );
      kernelComponent.computeAccumulation( ei, stack );
      kernelComponent.computeVolumeBalance( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// Number of fluid phases
  localIndex const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const m_volume;

  /// Views on the porosity
  arrayView2d< real64 const > const m_porosityOld;
  arrayView2d< real64 const > const m_porosityNew;
  arrayView2d< real64 const > const m_dPoro_dPres;

  /// Views on the derivatives of comp fractions wrt component density
  arrayView3d< real64 const, compflow::USD_COMP_DC > const m_dCompFrac_dCompDens;

  /// Views on the phase volume fractions (excluding derivative wrt temperature)
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFracOld;
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFrac;
  arrayView2d< real64 const, compflow::USD_PHASE > const m_dPhaseVolFrac_dPres;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > const m_dPhaseVolFrac_dCompDens;

  /// Views on the phase densities (excluding derivative wrt temperature)
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseDensOld;
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseDens;
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_dPhaseDens_dPres;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const m_dPhaseDens_dComp;

  /// Views on the phase component fraction (excluding derivative wrt temperature)
  arrayView3d< real64 const, compflow::USD_PHASE_COMP > const m_phaseCompFracOld;
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const m_phaseCompFrac;
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const m_dPhaseCompFrac_dPres;
  arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const m_dPhaseCompFrac_dComp;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

};

namespace internal
{

template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: type should be integral" );

  switch( value )
  {
    case 1:
    { lambda( std::integral_constant< T, 1 >() ); return; }
    case 2:
    { lambda( std::integral_constant< T, 2 >() ); return; }
    case 3:
    { lambda( std::integral_constant< T, 3 >() ); return; }
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return; }
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return; }
    default:
    { GEOSX_ERROR( "Unsupported number of components: " << value ); }
  }
}

} // namespace internal


/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] isIsothermal flag specifying whether the assembly is isothermal or non-isothermal
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( bool const isIsothermal,
                   localIndex const numComps,
                   localIndex const numPhases,
                   globalIndex const rankOffset,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   MultiFluidBase const & fluid,
                   CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
    {
      if( isIsothermal )
      {
        localIndex constexpr NUM_COMP = NC();
        localIndex constexpr NUM_DOF = NC()+1;
        ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
        kernel( numPhases, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
        ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::template
        launch< POLICY, ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF > >( subRegion.size(), kernel );
      }
      else
      {
        GEOSX_ERROR( "CompositionalMultiphaseBase: Thermal simulation is not supported yet: " );
        /*
           TODO: uncomment and move when the thermal kernel is ready
           localIndex constexpr NUM_COMP = NC();
           localIndex constexpr NUM_DOF = NC()+2;
           ThermalElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
           kernel( numPhases, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
           ThermalElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::template
           launch< POLICY, ThermalElementBasedAssemblyKernel< NUM_COMP, NUM_DOF > >( subRegion.size(), kernel );
         */
      }
    } );
  }

};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{

  template< typename POLICY, typename REDUCE_POLICY >
  static void launch( arrayView1d< real64 const > const & localResidual,
                      globalIndex const rankOffset,
                      localIndex const numComponents,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< integer const > const & ghostRank,
                      arrayView1d< real64 const > const & refPoro,
                      arrayView1d< real64 const > const & volume,
                      arrayView1d< real64 const > const & totalDensOld,
                      real64 & localResidualNorm )
  {
    RAJA::ReduceSum< REDUCE_POLICY, real64 > localSum( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;
        real64 const normalizer = totalDensOld[ei] * refPoro[ei] * volume[ei];

        for( localIndex idof = 0; idof < numComponents + 1; ++idof )
        {
          real64 const val = localResidual[localRow + idof] / normalizer;
          localSum += val * val;
        }
      }
    } );
    localResidualNorm += localSum.get();
  }

};


/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY, typename REDUCE_POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          localIndex const numComponents,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          integer const allowCompDensChopping,
          real64 const scalingFactor )
  {
    real64 constexpr eps = minDensForDivision;

    RAJA::ReduceMin< REDUCE_POLICY, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;
        {
          real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[localRow];
          check.min( newPres >= 0.0 );
        }

        // if component density chopping is not allowed, the time step fails if a component density is negative
        // otherwise, we just check that the total density is positive, and negative component densities
        // will be chopped (i.e., set to zero) in ApplySystemSolution)
        if( !allowCompDensChopping )
        {
          for( localIndex ic = 0; ic < numComponents; ++ic )
          {
            real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            check.min( newDens >= 0.0 );
          }
        }
        else
        {
          real64 totalDens = 0.0;
          for( localIndex ic = 0; ic < numComponents; ++ic )
          {
            real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            totalDens += (newDens > 0.0) ? newDens : 0.0;
          }
          check.min( totalDens >= eps );
        }
      }
    } );
    return check.get();
  }

};


/******************************** Kernel launch machinery ********************************/

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector1( localIndex const numComp, ARGS && ... args )
{
  internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
  {
    KERNELWRAPPER::template launch< NC() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector2( localIndex const numComp, localIndex const numPhase, ARGS && ... args )
{
  // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
  if( numPhase == 2 )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      KERNELWRAPPER::template launch< NC(), 2 >( std::forward< ARGS >( args ) ... );
    } );
  }
  else if( numPhase == 3 )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      KERNELWRAPPER::template launch< NC(), 3 >( std::forward< ARGS >( args ) ... );
    } );
  }
  else
  {
    GEOSX_ERROR( "Unsupported number of phases: " << numPhase );
  }
}

} // namespace CompositionalMultiphaseBaseKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP
