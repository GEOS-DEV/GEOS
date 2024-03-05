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
 * @file IsothermalCompositionalMultiphaseHybridFVMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/HybridFVMHelperKernels.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseHybridFVMKernels
{

using namespace constitutive;


/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_FACE number of faces per element
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @tparam IP the type of inner product
 * @brief Define the interface for the assembly kernel in charge of flux terms for the hybrid FVM scheme
 */
template< integer NUM_FACE, integer NUM_COMP, integer NUM_PHASE, typename IP >
class ElementBasedAssemblyKernel
{
public:

  /// Compile time value for the number of faces
  static constexpr integer numFace = NUM_FACE;

  /// Compile time value for the number of components
  static constexpr integer numComp = NUM_COMP;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /// Compile time value for the number of dofs
  static constexpr integer numDof = NUM_COMP+1;

  /// Min allowed value of the total mobility
  static constexpr real64 minTotalMobility = 1e-12;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using CompFlowAccessors =
    StencilAccessors< fields::flow::phaseMobility,
                      fields::flow::dPhaseMobility,
                      fields::flow::dGlobalCompFraction_dGlobalCompDensity >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::dPhaseDensity,
                              fields::multifluid::phaseMassDensity,
                              fields::multifluid::dPhaseMassDensity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;

  /**
   * @brief Constructor
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] er the index of the region
   * @param[in] esr the index of the subregion
   * @param[in] lengthTolerance tolerance used in the transmissibility computations
   * @param[in] faceDofKey the string key to retrieve the face degrees of freedom numbers
   * @param[in] nodeManager the node manager
   * @param[in] faceManager the face manager
   * @param[in] subRegion the element subregion
   * @param[in] dofNumberAccessor accessors for the element dof numbers
   * @param[in] compFlowAccessors accessors for the flow variables
   * @param[in] multiFluidAccessors accessors for the multifluid variables
   * @param[in] permeability the permeability model
   * @param[in] regionFilter the region filter
   * @param[in] dt the time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              localIndex const er,
                              localIndex const esr,
                              real64 const & lengthTolerance,
                              string const faceDofKey,
                              NodeManager const & nodeManager,
                              FaceManager const & faceManager,
                              CellElementSubRegion const & subRegion,
                              DofNumberAccessor const & dofNumberAccessor,
                              CompFlowAccessors const & compFlowAccessors,
                              MultiFluidAccessors const & multiFluidAccessors,
                              constitutive::PermeabilityBase const & permeability,
                              SortedArrayView< localIndex const > const & regionFilter,
                              real64 const & dt,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    :
    m_rankOffset( rankOffset ),
    m_er( er ),
    m_esr( esr ),
    m_lengthTolerance( lengthTolerance ),
    m_dt( dt ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_elemDofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_faceGhostRank( faceManager.ghostRank() ),
    m_faceDofNumber( faceManager.getReference< array1d< globalIndex > >( faceDofKey ) ),
    m_elemToFaces( subRegion.faceList().toViewConst() ),
    m_elemCenter( subRegion.getElementCenter() ),
    m_elemVolume( subRegion.getElementVolume() ),
    m_elemGravCoef( subRegion.getField< fields::flow::gravityCoefficient >() ),
    m_faceToNodes( faceManager.nodeList().toViewConst() ),
    m_faceGravCoef( faceManager.getField< fields::flow::gravityCoefficient >() ),
    m_mimFaceGravCoef( faceManager.getField< fields::flow::mimGravityCoefficient >() ),
    m_regionFilter( regionFilter ),
    m_nodePosition( nodeManager.referencePosition() ),
    m_elemRegionList( faceManager.elementRegionList() ),
    m_elemSubRegionList( faceManager.elementSubRegionList() ),
    m_elemList( faceManager.elementList() ),
    m_elemPerm( permeability.permeability() ),
    m_transMultiplier( faceManager.getField< fields::flow::transMultiplier >() ),
    m_elemPres( subRegion.getField< fields::flow::pressure >() ),
    m_facePres( faceManager.getField< fields::flow::facePressure >() ),
    m_phaseMob( compFlowAccessors.get( fields::flow::phaseMobility {} ) ),
    m_dPhaseMob( compFlowAccessors.get( fields::flow::dPhaseMobility {} ) ),
    m_dCompFrac_dCompDens( compFlowAccessors.get( fields::flow::dGlobalCompFraction_dGlobalCompDensity {} ) ),
    m_phaseDens( multiFluidAccessors.get( fields::multifluid::phaseDensity {} ) ),
    m_dPhaseDens( multiFluidAccessors.get( fields::multifluid::dPhaseDensity {} ) ),
    m_phaseMassDens( multiFluidAccessors.get( fields::multifluid::phaseMassDensity {} ) ),
    m_dPhaseMassDens( multiFluidAccessors.get( fields::multifluid::dPhaseMassDensity {} ) ),
    m_phaseCompFrac( multiFluidAccessors.get( fields::multifluid::phaseCompFraction {} ) ),
    m_dPhaseCompFrac( multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
    GEOS_HOST_DEVICE
    StackVariables()
      : transMatrix( numFace, numFace ),
      transMatrixGrav( numFace, numFace )
    {}

    // transmissibility matrices
    stackArray2d< real64, numFace *numFace > transMatrix;
    stackArray2d< real64, numFace *numFace > transMatrixGrav;

    // one sided volumetric fluxes
    real64 oneSidedVolFlux[numFace]{};
    real64 dOneSidedVolFlux_dPres[numFace]{};
    real64 dOneSidedVolFlux_dFacePres[numFace][numFace]{};
    real64 dOneSidedVolFlux_dCompDens[numFace][numComp]{};

    // divergence of fluxes
    real64 divMassFluxes[numComp]{};
    real64 dDivMassFluxes_dElemVars[numComp][numDof*(numFace+1)]{};
    real64 dDivMassFluxes_dFaceVars[numComp][numFace]{};

    // auxiliary variables for upwinding

    // upwinding phase buoyancy transport coefficients
    real64 upwPhaseViscCoef[numPhase][numComp]{};
    real64 dUpwPhaseViscCoef_dPres[numPhase][numComp]{};
    real64 dUpwPhaseViscCoef_dCompDens[numPhase][numComp][numComp]{};

    // gravity term: ( \rho_l - \rho_m ) g \Delta z
    real64 phaseGravTerm[numPhase][numPhase-1]{};
    real64 dPhaseGravTerm_dPres[numPhase][numPhase-1][2]{};
    real64 dPhaseGravTerm_dCompDens[numPhase][numPhase-1][2][numComp]{};

    // upwinding phase buoyancy transport coefficients
    real64 upwPhaseGravCoef[numPhase][numPhase-1][numComp]{};
    real64 dUpwPhaseGravCoef_dPres[numPhase][numPhase-1][numComp][2]{};
    real64 dUpwPhaseGravCoef_dCompDens[numPhase][numPhase-1][numComp][2][numComp]{};

    // dof numbers
    localIndex cellCenteredEqnRowIndex{};
    localIndex faceCenteredEqnRowIndex[numFace]{};
    globalIndex upwViscDofNumber{};
    globalIndex elemDofColIndices[ numDof*(numFace+1) ]{};
    globalIndex faceDofColIndices[ numFace ]{};
  };


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    stack.cellCenteredEqnRowIndex = m_elemDofNumber[m_er][m_esr][ei] - m_rankOffset;
    for( integer idof = 0; idof < numDof; ++idof )
    {
      // note that the other elemDofColIndices will be collected in the compute
      stack.elemDofColIndices[idof] = m_elemDofNumber[m_er][m_esr][ei] + idof;
    }
    for( integer iFaceLoc = 0; iFaceLoc < numFace; ++iFaceLoc )
    {
      stack.faceCenteredEqnRowIndex[iFaceLoc] = m_faceDofNumber[m_elemToFaces[ei][iFaceLoc]] - m_rankOffset;
      stack.faceDofColIndices[iFaceLoc] = m_faceDofNumber[m_elemToFaces[ei][iFaceLoc]];
    }
  }

  /**
   * @brief In a given element, compute the one-sided volumetric fluxes at this element's faces
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeGradient( localIndex const ei,
                        StackVariables & stack ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    real64 dPhaseMassDens_dC[numPhase][numComp]{};
    real64 dPresDif_dCompDens[numComp]{};
    real64 dPhaseGravDif_dCompDens[numComp]{};
    real64 dPhaseMobPotDif_dCompDens[numComp]{};

    // 0) precompute dPhaseDens_dC since it is always computed at the element center
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[m_er][m_esr][ei],
                      m_dPhaseMassDens[m_er][m_esr][ei][0][ip],
                      dPhaseMassDens_dC[ip],
                      Deriv::dC );
    }

    for( integer iFaceLoc = 0; iFaceLoc < numFace; ++iFaceLoc )
    {

      // now in the following nested loop,
      // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
      for( integer jFaceLoc = 0; jFaceLoc < numFace; ++jFaceLoc )
      {

        // depth difference between element center and face center
        real64 const ccGravCoef = m_elemGravCoef[ei];
        real64 const fGravCoef = m_faceGravCoef[m_elemToFaces[ei][jFaceLoc]];
        real64 const gravCoefDif = ccGravCoef - fGravCoef;

        for( integer ip = 0; ip < numPhase; ++ip )
        {

          // 1) compute the potential diff between the cell center and the face center
          real64 const ccPres = m_elemPres[ei];
          real64 const fPres  = m_facePres[m_elemToFaces[ei][jFaceLoc]];

          // pressure difference
          real64 const presDif = ccPres - fPres;
          real64 const dPresDif_dPres = 1;
          real64 const dPresDif_dFacePres = -1;
          for( integer ic = 0; ic < numComp; ++ic )
          {
            dPresDif_dCompDens[ic] = 0.0; // no capillary pressure
          }

          // gravity term
          real64 const phaseGravDif = m_phaseMassDens[m_er][m_esr][ei][0][ip] * gravCoefDif;
          real64 const dPhaseGravDif_dPres = m_dPhaseMassDens[m_er][m_esr][ei][0][ip][Deriv::dP] * gravCoefDif;
          for( localIndex ic = 0; ic < numComp; ++ic )
          {
            dPhaseGravDif_dCompDens[ic] = dPhaseMassDens_dC[ip][ic] * gravCoefDif;
          }
          // no density evaluated at the face center

          // potential difference
          real64 const phasePotDif = presDif - phaseGravDif;
          real64 const phaseMobPotDif = m_phaseMob[m_er][m_esr][ei][ip] * phasePotDif;
          real64 const dPhaseMobPotDif_dPres = m_dPhaseMob[m_er][m_esr][ei][ip][Deriv::dP] * phasePotDif
                                               + m_phaseMob[m_er][m_esr][ei][ip] * (dPresDif_dPres - dPhaseGravDif_dPres);
          real64 const dPhaseMobPotDif_dFacePres = m_phaseMob[m_er][m_esr][ei][ip] * dPresDif_dFacePres;
          for( integer ic = 0; ic < numComp; ++ic )
          {
            dPhaseMobPotDif_dCompDens[ic] = m_dPhaseMob[m_er][m_esr][ei][ip][Deriv::dC+ic] * phasePotDif
                                            + m_phaseMob[m_er][m_esr][ei][ip] * (dPresDif_dCompDens[ic] - dPhaseGravDif_dCompDens[ic]);
          }

          // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
          stack.oneSidedVolFlux[iFaceLoc] += stack.transMatrix[iFaceLoc][jFaceLoc] * phaseMobPotDif;
          stack.dOneSidedVolFlux_dPres[iFaceLoc] += stack.transMatrix[iFaceLoc][jFaceLoc] * dPhaseMobPotDif_dPres;
          stack.dOneSidedVolFlux_dFacePres[iFaceLoc][jFaceLoc] += stack.transMatrix[iFaceLoc][jFaceLoc] * dPhaseMobPotDif_dFacePres;
          for( integer ic = 0; ic < numComp; ++ic )
          {
            stack.dOneSidedVolFlux_dCompDens[iFaceLoc][ic] += stack.transMatrix[iFaceLoc][jFaceLoc] * dPhaseMobPotDif_dCompDens[ic];
          }
        }
      }
    }
  }

  /**
   * @brief In a given element, compute the upwinded term \lambda_p / \lambda_T at the face
   * @param[in] local the indices triplet of the local cell
   * @param[in] neighbor the indices triplet of the neighbor cell
   * @param[in] iFaceLoc the index of the face
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void upwindViscousCoefficient( localIndex const (&local)[3],
                                 localIndex const (&neighbor)[3],
                                 localIndex const iFaceLoc,
                                 StackVariables & stack ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    real64 dUpwMobRatio_dCompDens[numComp]{};
    real64 dUpwDensMobRatio_dCompDens[numComp]{};
    real64 dPhaseDens_dC[numComp]{};
    real64 dPhaseCompFrac_dC[numComp]{};

    // 1) Upwind
    localIndex const er  = ( stack.oneSidedVolFlux[iFaceLoc] > 0 ) ? local[0] : neighbor[0];
    localIndex const esr = ( stack.oneSidedVolFlux[iFaceLoc] > 0 ) ? local[1] : neighbor[1];
    localIndex const ei  = ( stack.oneSidedVolFlux[iFaceLoc] > 0 ) ? local[2] : neighbor[2];

    // 2) Compute total mobility: \lambda_T = \sum_{\ell} \lambda_{\ell}
    real64 totalMob = 0;
    real64 dTotalMob_dPres = 0;
    real64 dTotalMob_dCompDens[numComp]{};
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      totalMob += m_phaseMob[er][esr][ei][ip];
      dTotalMob_dPres += m_dPhaseMob[er][esr][ei][ip][Deriv::dP];
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dTotalMob_dCompDens[ic] += m_dPhaseMob[er][esr][ei][ip][Deriv::dC+ic];
      }
    }
    real64 const totalMobInv = 1.0 / totalMob;

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      // 3) Compute viscous mobility ratio: \frac{\lambda_{\ell}}{\lambda_T}
      real64 const upwMobRatio = m_phaseMob[er][esr][ei][ip] * totalMobInv;
      real64 const dUpwMobRatio_dPres = ( m_dPhaseMob[er][esr][ei][ip][Deriv::dP] - upwMobRatio * dTotalMob_dPres )
                                        * totalMobInv;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dUpwMobRatio_dCompDens[ic] = ( m_dPhaseMob[er][esr][ei][ip][Deriv::dC+ic] - upwMobRatio * dTotalMob_dCompDens[ic] )
                                     * totalMobInv;
      }

      // 4) Multiply mobility ratio by phase density: \rho^{up}_{\ell} \frac{\lambda_{\ell}}{\lambda_T}
      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[er][esr][ei],
                      m_dPhaseDens[er][esr][ei][0][ip],
                      dPhaseDens_dC,
                      Deriv::dC );
      real64 const upwDensMobRatio = m_phaseDens[er][esr][ei][0][ip] * upwMobRatio;
      real64 const dUpwDensMobRatio_dPres = m_dPhaseDens[er][esr][ei][0][ip][Deriv::dP] * upwMobRatio
                                            + m_phaseDens[er][esr][ei][0][ip] * dUpwMobRatio_dPres;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dUpwDensMobRatio_dCompDens[ic] = dPhaseDens_dC[ic] * upwMobRatio
                                         + m_phaseDens[er][esr][ei][0][ip] * dUpwMobRatio_dCompDens[ic];
      }

      // 5) Multiply density mobility ratio by phase comp fraction: x_{c,\ell} \rho^{up}_{\ell} \frac{\lambda_{\ell}}{\lambda_T}
      for( integer ic = 0; ic < numComp; ++ic )
      {
        applyChainRule( numComp,
                        m_dCompFrac_dCompDens[er][esr][ei],
                        m_dPhaseCompFrac[er][esr][ei][0][ip][ic],
                        dPhaseCompFrac_dC,
                        Deriv::dC );
        stack.upwPhaseViscCoef[ip][ic] = m_phaseCompFrac[er][esr][ei][0][ip][ic] * upwDensMobRatio;
        stack.dUpwPhaseViscCoef_dPres[ip][ic] = m_dPhaseCompFrac[er][esr][ei][0][ip][ic][Deriv::dP] * upwDensMobRatio
                                                + m_phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dPres;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          stack.dUpwPhaseViscCoef_dCompDens[ip][ic][jc] = dPhaseCompFrac_dC[jc] * upwDensMobRatio
                                                          + m_phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dCompDens[jc];
        }
      }
    }

    // 6) Save the dof number of the upwind cell
    stack.upwViscDofNumber = m_elemDofNumber[er][esr][ei];

  }

  /**
   * @brief In a given element, compute the term \lambda_p / \lambda_T u_T at the face
   * @param[in] local the indices triplet of the local cell
   * @param[in] neighbor the indices triplet of the neighbor cell
   * @param[in] iFaceLoc the index of the face
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void assembleViscousFlux( localIndex const (&local)[3],
                            localIndex const (&neighbor)[3],
                            localIndex const iFaceLoc,
                            StackVariables & stack ) const
  {
    localIndex const elemVarsOffset = numDof*(iFaceLoc+1);
    globalIndex const elemDofNumber = m_elemDofNumber[local[0]][local[1]][local[2]];
    globalIndex const neighborDofNumber = m_elemDofNumber[neighbor[0]][neighbor[1]][neighbor[2]];

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      for( integer ic = 0; ic < numComp; ++ic )
      {
        // compute the mass flux at the one-sided face plus its derivatives
        // add the newly computed flux to the sum

        real64 const dt_upwPhaseViscCoef = m_dt * stack.upwPhaseViscCoef[ip][ic];

        // residual
        stack.divMassFluxes[ic] += dt_upwPhaseViscCoef * stack.oneSidedVolFlux[iFaceLoc];

        // local derivatives
        stack.dDivMassFluxes_dElemVars[ic][0] += dt_upwPhaseViscCoef * stack.dOneSidedVolFlux_dPres[iFaceLoc];
        stack.dDivMassFluxes_dElemVars[ic][0] += ( elemDofNumber == stack.upwViscDofNumber )
                                                 * m_dt * stack.dUpwPhaseViscCoef_dPres[ip][ic] * stack.oneSidedVolFlux[iFaceLoc];
        for( integer jc = 0; jc < numComp; ++jc )
        {
          stack.dDivMassFluxes_dElemVars[ic][jc+1] += dt_upwPhaseViscCoef * stack.dOneSidedVolFlux_dCompDens[iFaceLoc][jc];
          stack.dDivMassFluxes_dElemVars[ic][jc+1] += ( elemDofNumber == stack.upwViscDofNumber )
                                                      * m_dt * stack.dUpwPhaseViscCoef_dCompDens[ip][ic][jc] * stack.oneSidedVolFlux[iFaceLoc];
        }

        // neighbor derivatives
        stack.dDivMassFluxes_dElemVars[ic][elemVarsOffset] += ( elemDofNumber != stack.upwViscDofNumber )
                                                              * m_dt * stack.dUpwPhaseViscCoef_dPres[ip][ic] * stack.oneSidedVolFlux[iFaceLoc];

        for( integer jc = 0; jc < numComp; ++jc )
        {
          stack.dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] += ( elemDofNumber != stack.upwViscDofNumber )
                                                                     * m_dt * stack.dUpwPhaseViscCoef_dCompDens[ip][ic][jc] * stack.oneSidedVolFlux[iFaceLoc];
        }

        for( integer jFaceLoc = 0; jFaceLoc < numFace; ++jFaceLoc )
        {
          stack.dDivMassFluxes_dFaceVars[ic][jFaceLoc] += dt_upwPhaseViscCoef * stack.dOneSidedVolFlux_dFacePres[iFaceLoc][jFaceLoc];
        }
      }
    }

    // collect the relevant dof numbers
    // this is not done in the setup function to avoid recomputations of cell connectivity
    for( integer idof = 0; idof < numDof; ++idof )
    {
      stack.elemDofColIndices[elemVarsOffset+idof] = neighborDofNumber + idof;
    }
  }

  /**
   * @brief In a given element, compute the term \lambda_p \lambda_m / \lambda_T u_T at the face
   * @param[in] local the indices triplet of the local cell
   * @param[in] neighbor the indices triplet of the neighbor cell
   * @param[in] iFaceLoc the index of the face
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void upwindBuoyancyCoefficient( localIndex const (&local)[3],
                                  localIndex const (&neighbor)[3],
                                  localIndex const iFaceLoc,
                                  StackVariables & stack ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    real64 const transGravCoef = (local[0] != neighbor[0] || local[1] != neighbor[1] || local[2] != neighbor[2])
                                 * stack.transMatrixGrav[iFaceLoc][iFaceLoc] * (m_elemGravCoef[local[2]] - m_mimFaceGravCoef[m_elemToFaces[local[2]][iFaceLoc]]);

    // 1) Compute the driving force: T ( \rho^{avg}_{\ell} - \rho^{avg}_m ) g \Delta z
    computePhaseGravTerm( local,
                          neighbor,
                          transGravCoef,
                          stack );

    // 2) Compute the total mobility: \lambda_T = \sum_{\ell} \lambda_{\ell}
    real64 totalMob = 0.0;
    real64 dTotalMob_dPres[2]{};
    real64 dTotalMob_dCompDens[2][numComp]{};

    computeUpwindedTotalMobility( local,
                                  neighbor,
                                  totalMob,
                                  dTotalMob_dPres,
                                  dTotalMob_dCompDens,
                                  stack );
    real64 const totalMobInv = 1.0 / totalMob;

    // 3) Compute the quantities \x_{up}_{c,p} \rho_p \frac{\lambda_p \lambda_m}{\lambda_T}
    real64 dMobRatio_dPres[2]{};
    real64 dMobRatio_dCompDens[2][numComp]{};
    real64 dDensMobRatio_dPres[2]{};
    real64 dDensMobRatio_dCompDens[2][numComp]{};

    real64 dPhaseDens_dC[numComp]{};
    real64 dPhaseCompFrac_dC[numComp]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      localIndex k = 0;
      for( integer jp = 0; jp < numPhase; ++jp )
      {
        if( ip == jp )
        {
          continue;
        }

        // 3.a) Upwinding using the gravity term
        localIndex eru, esru, eiu, posu; // upwind
        localIndex erd, esrd, eid, posd; // downwind
        setIndicesForMobilityProductUpwinding( local, neighbor,
                                               stack.phaseGravTerm[ip][k],
                                               eru, esru, eiu, posu,
                                               erd, esrd, eid, posd );

        // 3.b) Compute mobility ratio \frac{\lambda_l \lambda_m}{\lambda_T}
        real64 const mobRatio = m_phaseMob[eru][esru][eiu][ip] * m_phaseMob[erd][esrd][eid][jp]
                                * totalMobInv;
        dMobRatio_dPres[posu] = ( m_dPhaseMob[eru][esru][eiu][ip][Deriv::dP] * m_phaseMob[erd][esrd][eid][jp]
                                  - mobRatio * dTotalMob_dPres[posu] ) * totalMobInv;
        dMobRatio_dPres[posd] = ( m_dPhaseMob[erd][esrd][eid][jp][Deriv::dP] * m_phaseMob[eru][esru][eiu][ip]
                                  - mobRatio * dTotalMob_dPres[posd] ) * totalMobInv;

        for( integer ic = 0; ic < numComp; ++ic )
        {
          dMobRatio_dCompDens[posu][ic] = ( m_dPhaseMob[eru][esru][eiu][ip][Deriv::dC+ic] * m_phaseMob[erd][esrd][eid][jp]
                                            - mobRatio * dTotalMob_dCompDens[posu][ic] ) * totalMobInv;
          dMobRatio_dCompDens[posd][ic] = ( m_dPhaseMob[erd][esrd][eid][jp][Deriv::dC+ic] * m_phaseMob[eru][esru][eiu][ip]
                                            - mobRatio * dTotalMob_dCompDens[posd][ic] ) * totalMobInv;
        }

        // 3.c) Compute mobility ratio multiplied by upwinded phase density \rho_l \frac{\lambda_l \lambda_m}{\lambda_T}
        applyChainRule( numComp,
                        m_dCompFrac_dCompDens[eru][esru][eiu],
                        m_dPhaseDens[eru][esru][eiu][0][ip],
                        dPhaseDens_dC,
                        Deriv::dC );
        real64 const densMobRatio = m_phaseDens[eru][esru][eiu][0][ip] * mobRatio;
        dDensMobRatio_dPres[posu] = m_dPhaseDens[eru][esru][eiu][0][ip][Deriv::dP] * mobRatio
                                    + m_phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dPres[posu];
        dDensMobRatio_dPres[posd] = m_phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dPres[posd];
        for( integer ic = 0; ic < numComp; ++ic )
        {
          dDensMobRatio_dCompDens[posu][ic] = dPhaseDens_dC[ic] * mobRatio
                                              + m_phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dCompDens[posu][ic];
          dDensMobRatio_dCompDens[posd][ic] = m_phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dCompDens[posd][ic];
        }

        // 3.d) Compute the final gravity coefficient \x_{up}_{c,p} \rho_p \frac{\lambda_l \lambda_m}{\lambda_T}
        for( integer ic = 0; ic < numComp; ++ic )
        {
          applyChainRule( numComp,
                          m_dCompFrac_dCompDens[eru][esru][eiu],
                          m_dPhaseCompFrac[eru][esru][eiu][0][ip][ic],
                          dPhaseCompFrac_dC,
                          Deriv::dC );
          stack.upwPhaseGravCoef[ip][k][ic] = m_phaseCompFrac[eru][esru][eiu][0][ip][ic] * densMobRatio;
          stack.dUpwPhaseGravCoef_dPres[ip][k][ic][posu] = m_dPhaseCompFrac[eru][esru][eiu][0][ip][ic][Deriv::dP] * densMobRatio
                                                           + m_phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dPres[posu];
          stack.dUpwPhaseGravCoef_dPres[ip][k][ic][posd] = m_phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dPres[posd];

          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dUpwPhaseGravCoef_dCompDens[ip][k][ic][posu][jc] = dPhaseCompFrac_dC[jc] * densMobRatio
                                                                     + m_phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dCompDens[posu][jc];
            stack.dUpwPhaseGravCoef_dCompDens[ip][k][ic][posd][jc] = m_phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dCompDens[posd][jc];
          }
        }
        ++k;
      }
    }
  }

  /**
   * @brief In a given element, compute the term T ( \rho_p - \rho_m ) g dz at the face
   * @param[in] local the indices triplet of the local cell
   * @param[in] neighbor the indices triplet of the neighbor cell
   * @param[in] iFaceLoc the index of the face
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computePhaseGravTerm( localIndex const (&local)[ 3 ],
                             localIndex const (&neighbor)[ 3 ],
                             real64 const & transGravCoef,
                             StackVariables & stack ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    localIndex const er   = local[0];
    localIndex const esr  = local[1];
    localIndex const ei   = local[2];
    localIndex const ern  = neighbor[0];
    localIndex const esrn = neighbor[1];
    localIndex const ein  = neighbor[2];

    real64 dPhaseMassDens_dCLoc[numComp]{};
    real64 dPhaseMassDens_dCNeighbor[numComp]{};
    real64 dPhaseMassDens_dC[numComp]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[er][esr][ei],
                      m_dPhaseMassDens[er][esr][ei][0][ip],
                      dPhaseMassDens_dCLoc,
                      Deriv::dC );
      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[ern][esrn][ein],
                      m_dPhaseMassDens[ern][esrn][ein][0][ip],
                      dPhaseMassDens_dCNeighbor,
                      Deriv::dC );

      localIndex k = 0;
      for( integer jp = 0; jp < numPhase; ++jp )
      {
        if( ip == jp )
        {
          continue;
        }

        stack.phaseGravTerm[ip][k] = -( m_phaseMassDens[er][esr][ei][0][ip] + m_phaseMassDens[ern][esrn][ein][0][ip] );
        stack.phaseGravTerm[ip][k] += ( m_phaseMassDens[er][esr][ei][0][jp] + m_phaseMassDens[ern][esrn][ein][0][jp] );
        stack.phaseGravTerm[ip][k] *= 0.5 * transGravCoef;

        stack.dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] = ( -m_dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP] + m_dPhaseMassDens[er][esr][ei][0][jp][Deriv::dP] );
        stack.dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] *= 0.5 * transGravCoef;

        stack.dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] = ( -m_dPhaseMassDens[ern][esrn][ein][0][ip][Deriv::dP] + m_dPhaseMassDens[ern][esrn][ein][0][jp][Deriv::dP] );
        stack.dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] *= 0.5 * transGravCoef;

        for( integer ic = 0; ic < numComp; ++ic )
        {
          stack.dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] = -0.5 * transGravCoef * dPhaseMassDens_dCLoc[ic];
          stack.dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] = -0.5 * transGravCoef * dPhaseMassDens_dCNeighbor[ic];
        }
        applyChainRule( numComp,
                        m_dCompFrac_dCompDens[er][esr][ei],
                        m_dPhaseMassDens[er][esr][ei][0][jp],
                        dPhaseMassDens_dC,
                        Deriv::dC );
        for( integer ic = 0; ic < numComp; ++ic )
        {
          stack.dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] += 0.5 * transGravCoef * dPhaseMassDens_dC[ic];
        }
        applyChainRule( numComp,
                        m_dCompFrac_dCompDens[ern][esrn][ein],
                        m_dPhaseMassDens[ern][esrn][ein][0][jp],
                        dPhaseMassDens_dC,
                        Deriv::dC );
        for( integer ic = 0; ic < numComp; ++ic )
        {
          stack.dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] += 0.5 * transGravCoef * dPhaseMassDens_dC[ic];
        }
        ++k;
      }
    }
  }


  /**
   * @brief In a given element, compute the upwinded total mobility \lambda_T for the gravity term at the face
   * @param[in] local the indices triplet of the local cell
   * @param[in] neighbor the indices triplet of the neighbor cell
   * @param[out] totalMob the upwinded total mobility
   * @param[out] dTotalMob_dPres the upwinded total mobility derivative wrt pressure
   * @param[out] dTotalMob_dCompDens the upwinded total mobility derivative wrt component density
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeUpwindedTotalMobility( localIndex const (&local)[3],
                                     localIndex const (&neighbor)[3],
                                     real64 & totalMob,
                                     real64 ( & dTotalMob_dPres )[2],
                                     real64 ( & dTotalMob_dCompDens )[2][numComp],
                                     StackVariables & stack ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    localIndex totalMobIds[numPhase][3]{};
    localIndex totalMobPos[numPhase]{};

    // get the indices with which the mobilities are evaluated in the total mobility
    setIndicesForTotalMobilityUpwinding( local,
                                         neighbor,
                                         totalMobIds,
                                         totalMobPos,
                                         stack );

    // evaluate the total mobility
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      localIndex const er  = totalMobIds[ip][0];
      localIndex const esr = totalMobIds[ip][1];
      localIndex const ei  = totalMobIds[ip][2];
      localIndex const pos = totalMobPos[ip];

      totalMob += m_phaseMob[er][esr][ei][ip];
      dTotalMob_dPres[pos] += m_dPhaseMob[er][esr][ei][pos][Deriv::dP];
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dTotalMob_dCompDens[pos][ic] += m_dPhaseMob[er][esr][ei][ip][Deriv::dC+ic];
      }
    }

    // check that the total mobility is above the min
    if( totalMob < minTotalMobility )
    {
      totalMob = minTotalMobility;
    }
  }

  /**
   * @brief In a given element, grab the indices used for the upwinding of the mobility product \lambda_m \lambda_p
   * @param[in] local the indices triplet of the local cell
   * @param[in] neighbor the indices triplet of the neighbor cell
   * @param[in] gravTerm the gravity term ( \rho_p - \rho_m ) g dz at the face
   * @param[in] eru the region index of the upwind element
   * @param[in] ersu the subRegion index of the upwind element
   * @param[in] eiu the upwind element index
   * @param[in] posu the position of the upwind element (local or neighbor)
   * @param[in] erd the region index of the downwind element
   * @param[in] ersd the subRegion index of the downwind element
   * @param[in] eid the downwind element index
   * @param[in] posd the position of the downwind element (local or neighbor)x
   */
  GEOS_HOST_DEVICE
  void setIndicesForMobilityProductUpwinding( localIndex const (&local)[3],
                                              localIndex const (&neighbor)[3],
                                              real64 const & gravTerm,
                                              localIndex & eru, localIndex & esru, localIndex & eiu, localIndex & posu,
                                              localIndex & erd, localIndex & esrd, localIndex & eid, localIndex & posd ) const
  {
    if( gravTerm > 0 )
    {
      eru  = local[0]; esru = local[1]; eiu  = local[2];
      posu = Pos::LOCAL;
      erd  = neighbor[0]; esrd = neighbor[1]; eid  = neighbor[2];
      posd = Pos::NEIGHBOR;
    }
    else
    {
      eru  = neighbor[0]; esru = neighbor[1]; eiu  = neighbor[2];
      posu = Pos::NEIGHBOR;
      erd  = local[0]; esrd = local[1]; eid  = local[2];
      posd = Pos::LOCAL;
    }
  }

  /**
   * @brief In a given element, grab the indices used for the upwinding of the total mobility \lambda_T in the gravity term
   * @param[in] local the indices triplet of the local cell
   * @param[in] neighbor the indices triplet of the neighbor cell
   * @param[in] totalMobIds the indices for the evaluation of the upwinded total mobility
   * @param[in] totalMobPos the position of the elements in the upwinded total mobility (local or neighbor)
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setIndicesForTotalMobilityUpwinding( localIndex const (&local)[3],
                                            localIndex const (&neighbor)[3],
                                            localIndex ( & totalMobIds )[numPhase][3],
                                            localIndex ( & totalMobPos )[numPhase],
                                            StackVariables & stack ) const
  {
    if( numPhase == 2 )
    {
      if( stack.phaseGravTerm[0][0] > 0 )
      {
        totalMobIds[0][0] = local[0]; totalMobIds[0][1] = local[1]; totalMobIds[0][2] = local[2];
        totalMobPos[0] = Pos::LOCAL;
        totalMobIds[1][0] = neighbor[0]; totalMobIds[1][1] = neighbor[1]; totalMobIds[1][2] = neighbor[2];
        totalMobPos[1] = Pos::NEIGHBOR;
      }
      else
      {
        totalMobIds[0][0] = neighbor[0]; totalMobIds[0][1] = neighbor[1]; totalMobIds[0][2] = neighbor[2];
        totalMobPos[0] = Pos::NEIGHBOR;
        totalMobIds[1][0] = local[0]; totalMobIds[1][1] = local[1]; totalMobIds[1][2] = local[2];
        totalMobPos[1] = Pos::LOCAL;
      }
    }
    else if( numPhase == 3 )
    {
      // TODO Francois: this should be improved
      // currently this implements the algorithm proposed by SH Lee
      for( integer ip = 0; ip < numPhase; ++ip )
      {
        if( ( stack.phaseGravTerm[ip][0] >= 0 && stack.phaseGravTerm[ip][1] >= 0 ) || // includes the no-buoyancy case
            ( ( LvArray::math::abs( stack.phaseGravTerm[ip][0] ) >= LvArray::math::abs( stack.phaseGravTerm[ip][1] ) ) && stack.phaseGravTerm[ip][1] >= 0 ) ||
            ( ( LvArray::math::abs( stack.phaseGravTerm[ip][1] ) >= LvArray::math::abs( stack.phaseGravTerm[ip][0] ) ) && stack.phaseGravTerm[ip][0] >= 0 ) )
        {
          totalMobIds[ip][0] = local[0]; totalMobIds[ip][1] = local[1]; totalMobIds[ip][2] = local[2];
          totalMobPos[ip] = Pos::LOCAL;
        }
        else
        {
          totalMobIds[ip][0] = neighbor[0]; totalMobIds[ip][1] = neighbor[1]; totalMobIds[ip][2] = neighbor[2];
          totalMobPos[ip] = Pos::NEIGHBOR;
        }
      }
    }
  }

  /**
   * @brief In a given element, compute the T \lambda_p \lambda_m / \lambda_T (\rho_p - \rho_m ) g at the face
   * @param[in] local the indices triplet of the local cell
   * @param[in] neighbor the indices triplet of the neighbor cell
   * @param[in] iFaceLoc the face index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void assembleBuoyancyFlux( localIndex const (&local)[3],
                             localIndex const (&neighbor)[3],
                             localIndex const iFaceLoc,
                             StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( local, neighbor );

    localIndex const elemVarsOffset = numDof*(iFaceLoc+1);

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      for( integer jp = 0; jp < numPhase - 1; ++jp )
      {
        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const dt_upwPhaseGravCoef = m_dt * stack.upwPhaseGravCoef[ip][jp][ic];

          // residual
          stack.divMassFluxes[ic] += dt_upwPhaseGravCoef * stack.phaseGravTerm[ip][jp];

          // local derivatives
          stack.dDivMassFluxes_dElemVars[ic][0] += m_dt * stack.dUpwPhaseGravCoef_dPres[ip][jp][ic][Pos::LOCAL] * stack.phaseGravTerm[ip][jp];
          stack.dDivMassFluxes_dElemVars[ic][0] += dt_upwPhaseGravCoef * stack.dPhaseGravTerm_dPres[ip][jp][Pos::LOCAL];

          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dDivMassFluxes_dElemVars[ic][jc+1] += m_dt * stack.dUpwPhaseGravCoef_dCompDens[ip][jp][ic][Pos::LOCAL][jc] * stack.phaseGravTerm[ip][jp];
            stack.dDivMassFluxes_dElemVars[ic][jc+1] += dt_upwPhaseGravCoef * stack.dPhaseGravTerm_dCompDens[ip][jp][Pos::LOCAL][jc];
          }

          // neighbor derivatives
          stack.dDivMassFluxes_dElemVars[ic][elemVarsOffset] += m_dt * stack.dUpwPhaseGravCoef_dPres[ip][jp][ic][Pos::NEIGHBOR] * stack.phaseGravTerm[ip][jp];
          stack.dDivMassFluxes_dElemVars[ic][elemVarsOffset] += dt_upwPhaseGravCoef * stack.dPhaseGravTerm_dPres[ip][jp][Pos::NEIGHBOR];

          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] += m_dt * stack.dUpwPhaseGravCoef_dCompDens[ip][jp][ic][Pos::NEIGHBOR][jc] * stack.phaseGravTerm[ip][jp];
            stack.dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] += dt_upwPhaseGravCoef * stack.dPhaseGravTerm_dCompDens[ip][jp][Pos::NEIGHBOR][jc];
          }
        }
      }
    }
  }

  /**
   * @brief In a given element, assemble the mass conservation equation
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFluxDivergence( localIndex const ei,
                              StackVariables & stack ) const
  {
    // for each element, loop over the one-sided faces
    for( integer iFaceLoc = 0; iFaceLoc < numFace; ++iFaceLoc )
    {

      // 1) Find if there is a neighbor, and if there is, grab the indices of the neighbor element
      localIndex local[3] = { m_er, m_esr, ei };
      localIndex neighbor[3] = { m_er, m_esr, ei };
      hybridFVMKernels::CellConnectivity::
        isNeighborFound( local,
                         iFaceLoc,
                         m_elemRegionList,
                         m_elemSubRegionList,
                         m_elemList,
                         m_regionFilter,
                         m_elemToFaces[ei],
                         neighbor );

      // 2) *************** Assemble viscous terms ******************

      // 2.a) Compute the upwinded x_{c, \ell} \rho_{\ell} \frac{\lambda_{\ell}}{\lambda_T} for each phase at this face
      upwindViscousCoefficient( local, neighbor, iFaceLoc, stack );

      // 2.b) Add the \x_{c,\ell} \rho_{\ell} \frac{\lambda_{\ell}}{\lambda_T} q_T of this face to the divergence of the flux in this cell
      assembleViscousFlux( local, neighbor, iFaceLoc, stack );

      // 3) *************** Assemble buoyancy terms ******************

      // 3.a) Compute the upwinded x_{c, \ell} \rho_{\ell} \frac{\lambda_{\ell}\lambda_m}{\lambda_T}
      //      and (\rho_{\ell} - \rho_m) g \Delta z for each phase at this face
      upwindBuoyancyCoefficient( local, neighbor, iFaceLoc, stack );

      // 3.b) Add the buoyancy term of this face to the divergence of the flux in this cell
      assembleBuoyancyFlux( local, neighbor, iFaceLoc, stack );

    }
  }

  /**
   * @brief Compute the fluxes contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack,
                FUNC && kernelOp = NoOpFunc{} ) const
  {
    GEOS_UNUSED_VAR( ei, stack, kernelOp );

    real64 const perm[ 3 ] = { m_elemPerm[ei][0][0], m_elemPerm[ei][0][1], m_elemPerm[ei][0][2] };

    // recompute the local transmissibility matrix at each iteration
    // we can decide later to precompute transMatrix if needed
    IP::template compute< numFace >( m_nodePosition,
                                     m_transMultiplier,
                                     m_faceToNodes,
                                     m_elemToFaces[ei],
                                     m_elemCenter[ei],
                                     m_elemVolume[ei],
                                     perm,
                                     m_lengthTolerance,
                                     stack.transMatrix );

    // currently the gravity term in the transport scheme is treated as in MRST, that is, always with TPFA
    // this is why below we have to recompute the TPFA transmissibility in addition to the transmissibility matrix above
    // TODO: treat the gravity term with a consistent inner product
    mimeticInnerProduct::
      TPFAInnerProduct::compute< numFace >( m_nodePosition,
                                            m_transMultiplier,
                                            m_faceToNodes,
                                            m_elemToFaces[ei],
                                            m_elemCenter[ei],
                                            m_elemVolume[ei],
                                            perm,
                                            m_lengthTolerance,
                                            stack.transMatrixGrav );

    /*
     * compute auxiliary quantities at the one sided faces of this element:
     * 1) One-sided volumetric fluxes
     * 2) Upwinded mobilities
     */

    // for each one-sided face of the elem,
    // compute the volumetric flux using transMatrix
    computeGradient( ei, stack );

    // at this point, we know the local flow direction in the element
    // so we can upwind the transport coefficients (mobilities) at the one sided faces
    // ** this function needs non-local information **
    if( m_elemGhostRank[ei] < 0 )
    {

      /*
       * perform assembly in this element in two steps:
       * 1) mass conservation equations
       * 2) face constraints
       */

      // use the computed one sided vol fluxes and the upwinded mobilities
      // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
      computeFluxDivergence( ei, stack );
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    using namespace compositionalMultiphaseUtilities;

    // Step 1: assemble cell-centered residual and its derivatives

    // we are ready to assemble the local flux and its derivatives
    // no need for atomic adds - each row is assembled by a single thread

    if( m_elemGhostRank[ei] < 0 )
    {

      // Apply equation/variable change transformation(s)
      real64 work[numDof*(numFace+1)];
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numDof * (numFace+1), stack.dDivMassFluxes_dElemVars, work );
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numFace, stack.dDivMassFluxes_dFaceVars, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( numComp, stack.divMassFluxes );

      // we are ready to assemble the local flux and its derivatives
      // no need for atomic adds - each row is assembled by a single thread

      for( integer ic = 0; ic < numComp; ++ic )
      {
        localIndex const eqnRowLocalIndex =
          LvArray::integerConversion< localIndex >( stack.cellCenteredEqnRowIndex + ic );

        GEOS_ASSERT_GE( eqnRowLocalIndex, 0 );
        GEOS_ASSERT_GT( m_localMatrix.numRows(), eqnRowLocalIndex );

        // residual
        m_localRhs[eqnRowLocalIndex] += stack.divMassFluxes[ic];

        // jacobian -- derivative wrt elem centered vars
        m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                                    &stack.elemDofColIndices[0],
                                                                    &stack.dDivMassFluxes_dElemVars[0][0] + ic * numDof * (numFace+1),
                                                                    numDof * (numFace+1) );

        // jacobian -- derivatives wrt face centered vars
        m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                                    &stack.faceDofColIndices[0],
                                                                    &stack.dDivMassFluxes_dFaceVars[0][0] + ic * numFace,
                                                                    numFace );
      }
    }

    // Step 2: assemble face-centered residuals and their derivatives

    // fluxes
    real64 dFlux_dElemVars[numDof]{};

    // for each element, loop over the local (one-sided) faces
    for( integer iFaceLoc = 0; iFaceLoc < numFace; ++iFaceLoc )
    {
      if( m_faceGhostRank[m_elemToFaces[ei][iFaceLoc]] < 0 )
      {
        GEOS_ASSERT_GE( stack.faceCenteredEqnRowIndex[iFaceLoc], 0 );
        GEOS_ASSERT_GT( m_localMatrix.numRows(), stack.faceCenteredEqnRowIndex[iFaceLoc] );

        // residual
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[stack.faceCenteredEqnRowIndex[iFaceLoc]], stack.oneSidedVolFlux[iFaceLoc] );

        dFlux_dElemVars[0] = stack.dOneSidedVolFlux_dPres[iFaceLoc];
        for( integer ic = 0; ic < numComp; ++ic )
        {
          dFlux_dElemVars[ic+1] = stack.dOneSidedVolFlux_dCompDens[iFaceLoc][ic];
        }

        // jacobian -- derivative wrt local cell centered pressure term
        m_localMatrix.addToRow< parallelDeviceAtomic >( stack.faceCenteredEqnRowIndex[iFaceLoc],
                                                        &stack.elemDofColIndices[0],
                                                        &dFlux_dElemVars[0],
                                                        numDof );

        // jacobian -- derivatives wrt face pressure terms
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( stack.faceCenteredEqnRowIndex[iFaceLoc],
                                                                            &stack.faceDofColIndices[0],
                                                                            stack.dOneSidedVolFlux_dFacePres[iFaceLoc],
                                                                            numFace );
      }
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
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  // struct to specify local and neighbor derivatives
  struct Pos
  {
    static constexpr integer LOCAL = 0;
    static constexpr integer NEIGHBOR = 1;
  };

  /// offset for my MPI rank
  globalIndex const m_rankOffset;

  /// index of the region and sub-region
  localIndex const m_er;
  localIndex const m_esr;

  /// length tolerance
  real64 const m_lengthTolerance;

  /// time step size
  real64 const m_dt;

  /// ghost rank numbers
  arrayView1d< integer const > const m_elemGhostRank;
  ElementViewConst< arrayView1d< globalIndex const > > const m_elemDofNumber;
  arrayView1d< integer const > const m_faceGhostRank;
  arrayView1d< globalIndex const > const m_faceDofNumber;

  /// topological and geometrical data
  arrayView2d< localIndex const > const m_elemToFaces;
  arrayView2d< real64 const > const m_elemCenter;
  arrayView1d< real64 const > const m_elemVolume;
  arrayView1d< real64 const > const m_elemGravCoef;
  ArrayOfArraysView< localIndex const > const m_faceToNodes;
  arrayView1d< real64 const > const m_faceGravCoef;
  arrayView1d< real64 const > const m_mimFaceGravCoef;

  SortedArrayView< localIndex const > const m_regionFilter;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodePosition;
  arrayView2d< localIndex const > const m_elemRegionList;
  arrayView2d< localIndex const > const m_elemSubRegionList;
  arrayView2d< localIndex const > const m_elemList;

  /// permeability
  arrayView3d< real64 const > const m_elemPerm;
  arrayView1d< real64 const > const m_transMultiplier;

  /// pressure and fluid data
  arrayView1d< real64 const > const m_elemPres;
  arrayView1d< real64 const > const m_facePres;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseMob;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseMob;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const m_dCompFrac_dCompDens;
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseDens;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const m_dPhaseDens;
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseMassDens;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const m_dPhaseMassDens;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const m_phaseCompFrac;
  ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const m_dPhaseCompFrac;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

};

class ElementBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComps number of fluid components
   * @param[in] numPhases number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] er the index of the region
   * @param[in] esr the index of the subregion
   * @param[in] lengthTolerance tolerance used in the transmissibility computations
   * @param[in] elemDofKey the string key to retrieve the element dof numbers
   * @param[in] faceDofKey the string key to retrieve the face dof numbers
   * @param[in] solverName the name of the solver
   * @param[in] nodeManager the node manager
   * @param[in] faceManager the face manager
   * @param[in] elemManager the element region manager
   * @param[in] subRegion the element sub-region
   * @param[in] mimeticInnerProductBase the inner product to dispatch
   * @param[in] permeability the permeability model
   * @param[in] regionFilter the region filter
   * @param[in] dt the time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   localIndex const er,
                   localIndex const esr,
                   real64 const lengthTolerance,
                   string const elemDofKey,
                   string const faceDofKey,
                   string const solverName,
                   NodeManager const & nodeManager,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   CellElementSubRegion const & subRegion,
                   mimeticInnerProduct::MimeticInnerProductBase const & mimeticInnerProductBase,
                   constitutive::PermeabilityBase const & permeability,
                   SortedArrayView< localIndex const > const & regionFilter,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    mimeticInnerProductDispatch( mimeticInnerProductBase,
                                 [&] ( auto const mimeticInnerProduct )
    {
      using IP = TYPEOFREF( mimeticInnerProduct );

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + elemDofKey );

      #define LAUNCH( NC, NP ) \
        using kernelType = ElementBasedAssemblyKernel< NF(), NC, NP, IP >; \
        typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName ); \
        typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName ); \
        kernelType \
        kernel( rankOffset, er, esr, lengthTolerance, faceDofKey, nodeManager, faceManager, \
                subRegion, dofNumberAccessor, compFlowAccessors, multiFluidAccessors, permeability, \
                regionFilter, dt, localMatrix, localRhs ); \
        kernelType::template launch< POLICY >( subRegion.size(), kernel )


      // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
      // So, to avoid code duplication, we use the macro defined above
      if( numPhases == 2 )
      {
        if( numComps == 2 )
        {
          hybridFVMKernels::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NF )
          {
            LAUNCH( 2, 2 );
          } );
        }
        else if( numComps == 3 )
        {
          hybridFVMKernels::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NF )
          {
            LAUNCH( 3, 2 );
          } );
        }
        else if( numComps == 4 )
        {
          hybridFVMKernels::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NF )
          {
            LAUNCH( 4, 2 );
          } );
        }
        else if( numComps == 5 )
        {
          hybridFVMKernels::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NF )
          {
            LAUNCH( 5, 2 );
          } );
        }
        else
        {
          GEOS_ERROR( "Unsupported number of components: " << numComps );
        }
      }
      else if( numPhases == 3 )
      {
        if( numComps == 2 )
        {
          hybridFVMKernels::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NF )
          {
            LAUNCH( 2, 3 );
          } );
        }
        else if( numComps == 3 )
        {
          hybridFVMKernels::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NF )
          {
            LAUNCH( 3, 3 );
          } );
        }
        else if( numComps == 4 )
        {
          hybridFVMKernels::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NF )
          {
            LAUNCH( 4, 3 );
          } );
        }
        else if( numComps == 5 )
        {
          hybridFVMKernels::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NF )
          {
            LAUNCH( 5, 3 );
          } );
        }
        else
        {
          GEOS_ERROR( "Unsupported number of components: " << numComps );
        }
      }
      else
      {
        GEOS_ERROR( "Unsupported number of phases: " << numPhases );
      }
    } );
  }

};

/******************************** PhaseMobilityKernel ********************************/

/**
 * @class PhaseMobilityKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase mobilities
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseMobilityKernel : public isothermalCompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  PhaseMobilityKernel( ObjectManagerBase & subRegion,
                       MultiFluidBase const & fluid,
                       RelativePermeabilityBase const & relperm )
    : Base(),
    m_phaseVolFrac( subRegion.getField< fields::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::flow::dPhaseVolumeFraction >() ),
    m_dCompFrac_dCompDens( subRegion.getField< fields::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseVisc( fluid.phaseViscosity() ),
    m_dPhaseVisc( fluid.dPhaseViscosity() ),
    m_phaseRelPerm( relperm.phaseRelPerm() ),
    m_dPhaseRelPerm_dPhaseVolFrac( relperm.dPhaseRelPerm_dPhaseVolFraction() ),
    m_phaseMob( subRegion.getField< fields::flow::phaseMobility >() ),
    m_dPhaseMob( subRegion.getField< fields::flow::dPhaseMobility >() )
  {}

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
    using Deriv = multifluid::DerivativeOffset;

    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const phaseVisc = m_phaseVisc[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc = m_dPhaseVisc[ei][0];
    arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > const phaseRelPerm = m_phaseRelPerm[ei][0];
    arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > const dPhaseRelPerm_dPhaseVolFrac = m_dPhaseRelPerm_dPhaseVolFrac[ei][0];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dPhaseVolFrac[ei];
    arraySlice1d< real64, compflow::USD_PHASE - 1 > const phaseMob = m_phaseMob[ei];
    arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const dPhaseMob = m_dPhaseMob[ei];

    real64 dRelPerm_dC[numComp]{};
    real64 dVisc_dC[numComp]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      // compute the phase mobility only if the phase is present
      bool const phaseExists = (phaseVolFrac[ip] > 0);
      if( !phaseExists )
      {
        phaseMob[ip] = 0.;
        for( integer jc = 0; jc < numComp + 2; ++jc )
        {
          dPhaseMob[ip][jc] = 0.;
        }
        continue;
      }

      real64 const viscosity = phaseVisc[ip];
      real64 const dVisc_dP = dPhaseVisc[ip][Deriv::dP];
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseVisc[ip], dVisc_dC, Deriv::dC );

      real64 const relPerm = phaseRelPerm[ip];
      real64 dRelPerm_dP = 0.0;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dRelPerm_dC[ic] = 0.0;
      }

      for( integer jp = 0; jp < numPhase; ++jp )
      {
        real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[ip][jp];
        dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac[jp][Deriv::dP];

        for( integer jc = 0; jc < numComp; ++jc )
        {
          dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac[jp][Deriv::dC+jc];
        }
      }

      real64 const mobility = relPerm / viscosity;

      phaseMob[ip] = mobility;
      dPhaseMob[ip][Deriv::dP] = dRelPerm_dP / viscosity
                                 - mobility * dVisc_dP / viscosity;

      // compositional derivatives
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseMob[ip][Deriv::dC+jc] = dRelPerm_dC[jc] / viscosity
                                      - mobility * dVisc_dC[jc] / viscosity;
      }

      // call the lambda in the phase loop to allow the reuse of the relperm, viscosity, and mobility
      // possible use: assemble the derivatives wrt temperature
      phaseMobilityKernelOp( ip, phaseMob[ip], dPhaseMob[ip] );
    }
  }


protected:

  // inputs

  /// Views on the phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dPhaseVolFrac;
  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

  /// Views on the phase viscosities
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseVisc;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseVisc;

  /// Views on the phase relative permeabilities
  arrayView3d< real64 const, relperm::USD_RELPERM > m_phaseRelPerm;
  arrayView4d< real64 const, relperm::USD_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;

  // outputs

  /// Views on the phase mobilities
  arrayView2d< real64, compflow::USD_PHASE > m_phaseMob;
  arrayView3d< real64, compflow::USD_PHASE_DC > m_dPhaseMob;

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
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   MultiFluidBase const & fluid,
                   RelativePermeabilityBase const & relperm )
  {
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityKernel< NUM_COMP, 2 > kernel( subRegion, fluid, relperm );
        PhaseMobilityKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityKernel< NUM_COMP, 3 > kernel( subRegion, fluid, relperm );
        PhaseMobilityKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public solverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = solverBaseKernels::ResidualNormKernelBase< 1 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using CompFlowAccessors =
    StencilAccessors< fields::elementVolume >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< constitutive::MultiFluidBase,
                              fields::multifluid::totalDensity_n >;
  using PorosityAccessors =
    StencilMaterialAccessors< constitutive::PorosityBase,
                              fields::porosity::porosity_n >;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      SortedArrayView< localIndex const > const & regionFilter,
                      FaceManager const & faceManager,
                      CompFlowAccessors const & compFlowAccessors,
                      MultiFluidAccessors const & multiFluidAccessors,
                      PorosityAccessors const & porosityAccessors,
                      real64 const & dt,
                      real64 const & minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_dt( dt ),
    m_regionFilter( regionFilter ),
    m_elemRegionList( faceManager.elementRegionList() ),
    m_elemSubRegionList( faceManager.elementSubRegionList() ),
    m_elemList( faceManager.elementList() ),
    m_volume( compFlowAccessors.get( fields::elementVolume {} ) ),
    m_porosity_n( porosityAccessors.get( fields::porosity::porosity_n {} ) ),
    m_totalDens_n( multiFluidAccessors.get( fields::multifluid::totalDensity_n {} ) )
  {}

  GEOS_HOST_DEVICE
  void computeMassNormalizer( localIndex const kf,
                              real64 & massNormalizer ) const
  {
    integer elemCounter = 0;

    for( integer k = 0; k < m_elemRegionList.size( 1 ); ++k )
    {
      localIndex const er  = m_elemRegionList[kf][k];
      localIndex const esr = m_elemSubRegionList[kf][k];
      localIndex const ei  = m_elemList[kf][k];
      bool const onBoundary = (er == -1 || esr == -1 || ei == -1);
      bool const isInTarget = m_regionFilter.contains( er );

      // if not on boundary, increment the normalizer
      if( !onBoundary && isInTarget )
      {
        massNormalizer += m_totalDens_n[er][esr][ei][0] * m_porosity_n[er][esr][ei][0] * m_volume[er][esr][ei];
        elemCounter++;
      }
    }
    massNormalizer /= elemCounter; // average mass in the adjacent cells at the previous converged time step
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const kf,
                            LinfStackVariables & stack ) const override
  {
    // if the face is adjacent to target region, compute the local values
    if( m_dofNumber[kf] >= 0 )
    {
      real64 massNormalizer = 0.0;
      computeMassNormalizer( kf, massNormalizer );

      // scaled residual to be in mass units (needed because element and face residuals are blended in a single norm)
      stack.localValue[0] += LvArray::math::abs( m_localResidual[stack.localRow] * m_dt ) / LvArray::math::max( m_minNormalizer, massNormalizer );
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const kf,
                          L2StackVariables & stack ) const override
  {
    // if the face is adjacent to target region, compute the local values
    if( m_dofNumber[kf] >= 0 )
    {
      real64 massNormalizer = 0;
      computeMassNormalizer( kf, massNormalizer );

      // scaled residual to be in mass units (needed because element and face residuals are blended in a single norm)
      real64 const valMass = m_localResidual[stack.localRow] * m_dt;
      stack.localValue[0] += valMass * valMass;
      stack.localNormalizer[0] += massNormalizer;
    }
  }


protected:

  /// Time step size
  real64 const m_dt;

  /// Filter to identify the target regions of the solver
  SortedArrayView< localIndex const > const m_regionFilter;

  /// Views on the maps face to elements
  arrayView2d< localIndex const > const m_elemRegionList;
  arrayView2d< localIndex const > const m_elemSubRegionList;
  arrayView2d< localIndex const > const m_elemList;

  /// View on the volume
  ElementViewConst< arrayView1d< real64 const > > const m_volume;

  /// View on porosity at the previous converged time step
  ElementViewConst< arrayView2d< real64 const > > const m_porosity_n;

  /// View on total mass/molar density at the previous converged time step
  ElementViewConst< arrayView2d< real64 const, multifluid::USD_FLUID > > const m_totalDens_n;

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
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] regionFilter filter to identify the target regions of the solver
   * @param[in] solverName the name of the solver
   * @param[in] elemManager reference to the element region manager
   * @param[in] faceManager reference to the face manager
   * @param[in] dt time step size
   * @param[in] minNormalizer the min normalizer
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( solverBaseKernels::NormType const normType,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   SortedArrayView< localIndex const > const & regionFilter,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   FaceManager const & faceManager,
                   real64 const & dt,
                   real64 const & minNormalizer,
                   real64 (& residualNorm)[1],
                   real64 (& residualNormalizer)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = faceManager.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = faceManager.ghostRank();

    using kernelType = ResidualNormKernel;
    typename kernelType::CompFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::MultiFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::PorosityAccessors poroAccessors( elemManager, solverName );

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank,
                               regionFilter, faceManager, flowAccessors, fluidAccessors, poroAccessors, dt, minNormalizer );
    if( normType == solverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( faceManager.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( faceManager.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{

  template< typename POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & facePres,
          real64 const scalingFactor )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const iface )
    {
      if( ghostRank[iface] < 0 && dofNumber[iface] >= 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[iface] - rankOffset );
        {
          real64 const newFacePres = facePres[iface] + scalingFactor * localSolution[localRow];
          check.min( newFacePres >= 0.0 );
        }
      }
    } );
    return check.get();
  }

};

/******************************** PrecomputeKernel ********************************/

struct PrecomputeKernel
{

  template< typename IP_TYPE, integer NF >
  static void
  launch( localIndex const subRegionSize,
          localIndex const faceManagerSize,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView2d< real64 const > const & elemCenter,
          arrayView1d< real64 const > const & elemVolume,
          arrayView3d< real64 const > const & elemPerm,
          arrayView1d< real64 const > const & elemGravCoef,
          arrayView2d< localIndex const > const & elemToFaces,
          arrayView1d< real64 const > const & transMultiplier,
          real64 const & lengthTolerance,
          arrayView1d< RAJA::ReduceSum< serialReduce, real64 > > const & mimFaceGravCoefNumerator,
          arrayView1d< RAJA::ReduceSum< serialReduce, real64 > > const & mimFaceGravCoefDenominator,
          arrayView1d< real64 > const & mimFaceGravCoef )
  {
    forAll< serialPolicy >( subRegionSize, [=] ( localIndex const ei )
    {
      stackArray2d< real64, NF *NF > transMatrix( NF, NF );

      real64 const perm[ 3 ] = { elemPerm[ei][0][0], elemPerm[ei][0][1], elemPerm[ei][0][2] };

      IP_TYPE::template compute< NF >( nodePosition,
                                       transMultiplier,
                                       faceToNodes,
                                       elemToFaces[ei],
                                       elemCenter[ei],
                                       elemVolume[ei],
                                       perm,
                                       lengthTolerance,
                                       transMatrix );

      for( integer ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
      {
        mimFaceGravCoefNumerator[elemToFaces[ei][ifaceLoc]] += elemGravCoef[ei] * transMatrix[ifaceLoc][ifaceLoc];
        mimFaceGravCoefDenominator[elemToFaces[ei][ifaceLoc]] += transMatrix[ifaceLoc][ifaceLoc];
      }
    } );

    forAll< serialPolicy >( faceManagerSize, [=] ( localIndex const iface )
    {
      if( !isZero( mimFaceGravCoefDenominator[iface].get() ) )
      {
        mimFaceGravCoef[iface] = mimFaceGravCoefNumerator[iface].get() / mimFaceGravCoefDenominator[iface].get();
      }
    } );
  }
};

} // namespace isothermalCompositionalMultiphaseHybridFVMKernels

} // namespace geosx

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
