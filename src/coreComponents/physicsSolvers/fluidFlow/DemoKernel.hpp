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
 * @file DemoKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_DEMOKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_DEMOKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "CompositionalMultiphaseUtilities.hpp"

namespace geosx
{

namespace DemoKernels
{

using namespace constitutive;


/**
 * @class AccumulationDemoKernel
 * @brief Define the interface for the isothermal accumulation kernel
 *
 * Will ultimately also inherit from a FV cell-based base kernel (reused for VolBal, etc)
 */
template< localIndex NUM_COMP, localIndex NUM_DOF >
class AccumulationDemoKernel
{
public:

  /// Compile time value for the number of components
  static constexpr localIndex NC = NUM_COMP;

  /// Compute time value for the number of degrees of freedom
  static constexpr localIndex NDOF = NUM_DOF;

  /**
   * @brief Constructor
   */
  AccumulationDemoKernel( localIndex const numPhases,
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
                          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dCompDens,
                          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseDensOld,
                          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
                          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
                          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
                          arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & phaseCompFracOld,
                          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & phaseCompFrac,
                          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & dPhaseCompFrac_dPres,
                          arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dPhaseCompFrac_dComp,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
    : m_numPhases( numPhases ),
    m_rankOffset( rankOffset ),
    m_dofNumber( dofNumber ),
    m_elemGhostRank( elemGhostRank ),
    m_volume( volume ),
    m_porosityOld( porosityOld ),
    m_porosityNew( porosityNew ),
    m_dPoro_dPres( dPoro_dPres ),
    m_dCompFrac_dCompDens( dCompFrac_dCompDens ),
    m_phaseVolFracOld( phaseVolFracOld ),
    m_phaseVolFrac( phaseVolFrac ),
    m_dPhaseVolFrac_dPres( dPhaseVolFrac_dPres ),
    m_dPhaseVolFrac_dCompDens( dPhaseVolFrac_dCompDens ),
    m_phaseDensOld( phaseDensOld ),
    m_phaseDens( phaseDens ),
    m_dPhaseDens_dPres( dPhaseDens_dPres ),
    m_dPhaseDens_dComp( dPhaseDens_dComp ),
    m_phaseCompFracOld( phaseCompFracOld ),
    m_phaseCompFrac( phaseCompFrac ),
    m_dPhaseCompFrac_dPres( dPhaseCompFrac_dPres ),
    m_dPhaseCompFrac_dComp( dPhaseCompFrac_dComp ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {
    // or maybe now we want to get the views from here by passing the subRegion
  }


  /**
   * @struct StackVariables
   * @brief Kernel variables allocated on the stack
   */
  struct StackVariables
  {
public:

    GEOSX_HOST_DEVICE
    StackVariables() {}

    /// local Jacobian information

    /// Index of the local row corresponding to this element
    localIndex localRow = -1;

    /// Indices of the matrix rows/columns corresponding to the dofs in this element
    globalIndex dofIndices[NDOF]{};

    /// C-array storage for the element local residual vector (all equations except volume balance)
    real64 localResidual[NDOF-1]{};

    /// C-array storage for the element local Jacobian matrix (all equations except volume balance, all dogs)
    real64 localJacobian[NDOF-1][NDOF]{};

  };


  /**
   * @brief Performs the setup phase for the kernel.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    // set DOF indices for this block
    stack.localRow = m_dofNumber[ei] - m_rankOffset;
    for( localIndex idof = 0; idof < NDOF; ++idof )
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
     * @brief operator() no-op used for additional terms in the residual / Jacobian
     */
    template< localIndex NC >
    GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE constexpr
    void operator() ( localIndex const ip,
                      real64 const & poreVolumeNew,
                      real64 const & poreVolumeOld,
                      real64 const & dPoreVolume_dPres,
                      real64 const & phaseAmountNew,
                      real64 const & phaseAmountOld,
                      real64 const & dPhaseAmount_dP,
                      real64 const (&dPhaseAmount_dC )[ NC ] )
    {
      GEOSX_UNUSED_VAR( ip,
                        poreVolumeNew, poreVolumeOld, dPoreVolume_dPres,
                        phaseAmountNew, phaseAmountOld, dPhaseAmount_dP, dPhaseAmount_dC );
    }

  };

  /**
   * @brief Compute the local contributions to the residual and Jacobian
   */
  template< typename LAMBDA = NoOpFunctors >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void compute( localIndex const ei,
                StackVariables & stack,
                LAMBDA && lambda = NoOpFunctors{} ) const
  {
    // compute pore volumes
    real64 const poreVolumeNew = m_volume[ei] * m_porosityNew[ei][0];
    real64 const poreVolumeOld = m_volume[ei] * m_porosityOld[ei][0];
    real64 const dPoreVolume_dPres = m_volume[ei] * m_dPoro_dPres[ei][0];

    // construct the slices (did not find a good way to create them in StackVariables)
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
    for( localIndex ip = 0; ip < m_numPhases; ++ip )
    {
      real64 const phaseAmountNew = poreVolumeNew * phaseVolFrac[ip] * phaseDens[ip];
      real64 const phaseAmountOld = poreVolumeOld * phaseVolFracOld[ip] * phaseDensOld[ip];

      real64 const dPhaseAmount_dP = dPoreVolume_dPres * phaseVolFrac[ip] * phaseDens[ip]
                                     + poreVolumeNew * (dPhaseVolFrac_dPres[ip] * phaseDens[ip]
                                                        + phaseVolFrac[ip] * dPhaseDens_dPres[ip]);

      // assemble density dependence
      applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dPhaseAmount_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                              + phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
        dPhaseAmount_dC[jc] *= poreVolumeNew;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for( localIndex ic = 0; ic < NC; ++ic )
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
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                             + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
          stack.localJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
        }
      }

      // call the lambda in the middle of the phase loop to reuse all the phase amounts and their derivatives
      // inside the lambda, assemble the derivatives wrt temperature, and the accumulation term of the energy equation for this phase
      lambda( ip,
              poreVolumeNew, poreVolumeOld, dPoreVolume_dPres,
              phaseAmountNew, phaseAmountOld, dPhaseAmount_dP, dPhaseAmount_dC );

    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( ei );
    
    using namespace CompositionalMultiphaseUtilities;

    // apply equation/variable change transformation(s) for the component mass balance equations
    // this should work whether NDOF = NC + 1 or NDOF = NC + 2
    real64 work[NDOF];
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NDOF, stack.localJacobian, work );
    shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, stack.localResidual );

    // add contribution to residual and jacobian into mass balance equations
    // when called by the derived class, this will assemble the derivatives wrt temp in the mass balance equations
    // but this will not touch the energy equations
    for( localIndex i = 0; i < NC; ++i )
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
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// Number of fluid phases
  localIndex const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const & m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const & m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const & m_volume;

  /// Views on the porosity
  arrayView2d< real64 const > const & m_porosityOld;
  arrayView2d< real64 const > const & m_porosityNew;
  arrayView2d< real64 const > const & m_dPoro_dPres;

  /// Views on the derivatives of comp fractions wrt component density
  arrayView3d< real64 const, compflow::USD_COMP_DC > const & m_dCompFrac_dCompDens;

  /// Views on the phase volume fractions (excluding derivative wrt temperature)
  arrayView2d< real64 const, compflow::USD_PHASE > const & m_phaseVolFracOld;
  arrayView2d< real64 const, compflow::USD_PHASE > const & m_phaseVolFrac;
  arrayView2d< real64 const, compflow::USD_PHASE > const & m_dPhaseVolFrac_dPres;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > const & m_dPhaseVolFrac_dCompDens;

  /// Views on the phase densities (excluding derivative wrt temperature)
  arrayView2d< real64 const, compflow::USD_PHASE > const & m_phaseDensOld;
  arrayView3d< real64 const, multifluid::USD_PHASE > const & m_phaseDens;
  arrayView3d< real64 const, multifluid::USD_PHASE > const & m_dPhaseDens_dPres;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & m_dPhaseDens_dComp;

  /// Views on the phase component fraction (excluding derivative wrt temperature)
  arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & m_phaseCompFracOld;
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & m_phaseCompFrac;
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & m_dPhaseCompFrac_dPres;
  arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & m_dPhaseCompFrac_dComp;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const & m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const & m_localRhs;

};

} // namespace DemoKernels

} // namespace geosx

#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_DEMOKERNELS_HPP
