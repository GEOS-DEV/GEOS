/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file DruckerPragerExtended.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGEREXTENDED_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGEREXTENDED_HPP_

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class DruckerPragerExtendedUpdates
 *
 * Class to provide material updates that may be called from a kernel function.
 */
class DruckerPragerExtendedUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] initialFriction The ArrayView holding the initial friction data for each element.
   * @param[in] residualFriction The ArrayView holding the residual friction data for each element.
   * @param[in] dilationRatio The ArrayView holding the ratio between dilation and friction data for each element.
   * @param[in] pressureIntercept The ArrayView holding the pressure intercept for each element.
   * @param[in] hardening The ArrayView holding the hardening data for each element.
   * @param[in] state The ArrayView holding the state data for each element.
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] stress The ArrayView holding the stress data for each quadrature point.
   */
  DruckerPragerExtendedUpdates( arrayView1d< real64 const > const & initialFriction,
                                arrayView1d< real64 const > const & residualFriction,
                                arrayView1d< real64 const > const & dilationRatio,
                                arrayView1d< real64 const > const & pressureIntercept,
                                arrayView1d< real64 const > const & hardening,
                                arrayView2d< real64 > const & state,
                                arrayView1d< real64 const > const & bulkModulus,
                                arrayView1d< real64 const > const & shearModulus,
                                arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                arrayView3d< real64, solid::STRESS_USD > const & oldStress ):
    ElasticIsotropicUpdates( bulkModulus, shearModulus, newStress, oldStress ),
    m_initialFriction( initialFriction ),
    m_residualFriction( residualFriction ),
    m_dilationRatio( dilationRatio ),
    m_pressureIntercept( pressureIntercept ),
    m_hardening( hardening ),
    m_state( state )
  {}

  /// Default copy constructor
  DruckerPragerExtendedUpdates( DruckerPragerExtendedUpdates const & ) = default;

  /// Default move constructor
  DruckerPragerExtendedUpdates( DruckerPragerExtendedUpdates && ) = default;

  /// Deleted default constructor
  DruckerPragerExtendedUpdates() = delete;

  /// Deleted copy assignment operator
  DruckerPragerExtendedUpdates & operator=( DruckerPragerExtendedUpdates const & ) = delete;

  /// Deleted move assignment operator
  DruckerPragerExtendedUpdates & operator=( DruckerPragerExtendedUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; 

  // Bring in base implementations to prevent hiding warnings
  using ElasticIsotropicUpdates::smallStrainUpdate;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  real64 ( &stiffness )[6][6] ) const override final;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const final;

  real64 yield( localIndex const k,
                localIndex const GEOSX_UNUSED_PARAM( q ),
                real64 const invP,
                real64 const invQ,
                real64 const friction ) const
  {
    return  invQ + friction * ( invP - m_pressureIntercept[k] );
  }

  void yieldDerivatives( localIndex const k,
                         localIndex const GEOSX_UNUSED_PARAM( q ),
                         real64 const invP,
                         real64 const GEOSX_UNUSED_PARAM( invQ ),
                         real64 const friction,
                         real64 (& dF)[3] ) const
  {

    // The yield function is: invQ + friction * ( invP - m_pressureIntercept[k] )

    real64 dF_dP = friction;
    real64 dF_dQ = 1.0;
    real64 dF_dFriction = invP - m_pressureIntercept[k];

    dF[0] = dF_dP;
    dF[1] = dF_dQ;
    dF[2] = dF_dFriction;
  }

  void potentialDerivatives( localIndex const k,
                             localIndex const GEOSX_UNUSED_PARAM( q ),
                             real64 const GEOSX_UNUSED_PARAM( invP ),
                             real64 const GEOSX_UNUSED_PARAM( invQ ),
                             real64 const friction,
                             real64 (& dG)[8] ) const
  {

    // The plastic potential function is: G = invQ + invP * m_dilationRatio[k] * friction

    real64 dG_dP = m_dilationRatio[k] * friction;
    real64 dG_dQ = 1.0;

    real64 dG_dP_dP = 0.0;
    real64 dG_dP_dQ = 0.0;
    real64 dG_dP_dH = m_dilationRatio[k];

    real64 dG_dQ_dP = 0.0;
    real64 dG_dQ_dQ = 0.0;
    real64 dG_dQ_dH = 0.0;

    dG[0] = dG_dP;
    dG[1] = dG_dQ;

    dG[2] = dG_dP_dP;
    dG[3] = dG_dP_dQ;
    dG[4] = dG_dP_dH;

    dG[5] = dG_dQ_dP;
    dG[6] = dG_dQ_dQ;
    dG[7] = dG_dQ_dH;
  }

  real64 hardening( localIndex const k,
                    localIndex const GEOSX_UNUSED_PARAM( q ),
                    real64 const state ) const
  {
    if( state<1e-9 )
    {
      return m_initialFriction[k];
    }
    else
    {
      return m_initialFriction[k] + ( m_residualFriction[k] - m_initialFriction[k] ) * state / ( m_hardening[k] + state );
    }
  }

  real64 hardeningDerivatives( localIndex const k,
                               localIndex const GEOSX_UNUSED_PARAM( q ),
                               real64 const state ) const
  {
    if( state<1e-9 )
    {
      return 0.0;
    }
    else
    {
      return ( m_residualFriction[k] - m_initialFriction[k] ) * m_hardening[k] / ( m_hardening[k] + state ) / ( m_hardening[k] + state );
    }
  }

  real64 tmpState( localIndex const k,
                   localIndex const q,
                   real64 const GEOSX_UNUSED_PARAM( invP ),
                   real64 const GEOSX_UNUSED_PARAM( invQ ),
                   real64 const plasticMultiplier ) const
  {

    // The temporal state variable is updated inside the Newton loops by: state += plasticMultiplier
    // starting from the previous saved state

    return m_state[k][q] + plasticMultiplier;
  }

  real64 stateDerivatives( localIndex const GEOSX_UNUSED_PARAM( k ),
                           localIndex const GEOSX_UNUSED_PARAM( q ),
                           real64 const GEOSX_UNUSED_PARAM( invP ),
                           real64 const GEOSX_UNUSED_PARAM( invQ ),
                           real64 const GEOSX_UNUSED_PARAM( plasticMultiplier ) ) const
  {
    return 1.0;
  }

  real64 getStateVariable( localIndex const k,
                           localIndex const q ) const
  {
    return m_state[k][q];
  }

 void saveStateVariable( localIndex const k,
                         localIndex const q,
                         real64 const state ) const
  {
    m_state[k][q] = state;
  }

private:

  /// A reference to the ArrayView holding the initial friction coefficient for each element.
  arrayView1d< real64 const > const m_initialFriction;

  /// A reference to the ArrayView holding the residual friction coefficient for each element.
  arrayView1d< real64 const > const m_residualFriction;

  /// A reference to the ArrayView holding the dilation ratio for each element.
  arrayView1d< real64 const > const m_dilationRatio;

  /// A reference to the ArrayView holding the pressure intercept for each element.
  arrayView1d< real64 const > const m_pressureIntercept;

  /// A reference to the ArrayView holding the hardening parameter for each element.
  arrayView1d< real64 const > const m_hardening;

  /// A reference to the ArrayView holding the state variable for each integration point
  arrayView2d< real64 > const m_state;
};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerExtendedUpdates::smallStrainUpdate( localIndex const k,
                                                      localIndex const q,
                                                      real64 const ( &strainIncrement )[6],
                                                      real64 ( & stress )[6],
                                                      real64 ( & stiffness )[6][6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)

  ElasticIsotropicUpdates::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

  // decompose into mean (P) and von mises (Q) stress invariants

  real64 trialP;
  real64 trialQ;
  real64 deviator[6];

  twoInvariant::stressDecomposition( stress,
                                     trialP,
                                     trialQ,
                                     deviator );

  // Plastic functions and their derivatives

  real64 dF[3], dG[8], hardeningParam, stateVariable, dHardeningParam_dLambda;

  // check yield function F <= 0, using old state

  stateVariable = getStateVariable( k, q );

  hardeningParam = hardening( k, q, stateVariable );

  if( yield( k, q, trialP, trialQ, hardeningParam ) < 1e-9 ) // elasticity
  {
    return;
  }

  // else, plasticity (trial stress point lies outside yield surface)
  // the return mapping can in general be written as a newton iteration.

  real64 solution[3], residual[3], delta[3];
  real64 jacobian[3][3] = {{}}, jacobianInv[3][3] = {{}};

  solution[0] = trialP; // initial guess for newP
  solution[1] = trialQ; // initial guess for newQ
  solution[2] = 1e-5;   // initial guess for plastic multiplier

  real64 norm, normZero = 1e30;

  // begin newton loop

  for( localIndex iter=0; iter<20; ++iter )
  {

    // Derivatives of the yield function
    // dF_dP = dF[0], dF_dQ = dF[1], dF_dHardeningParam = dF[2]

    yieldDerivatives( k, q, solution[0], solution[1], hardeningParam, dF );

    // Derivatives of the plastic potential function
    // dG_dP = dG[0], dG_dQ = dG[1]
    // dG_dP_dP = dG[2], dG_dP_dQ = dG[3], dG_dP_dHardeningParam = dG[4] 
    // dG_dQ_dP = dG[5], dG_dQ_dQ = dG[6], dG_dQ_dHardeningParam = dG[7]

    potentialDerivatives( k, q, solution[0], solution[1], hardeningParam, dG );
 
    // Hardening parameter and its derivative to the plastic multiplier
     
    stateVariable = tmpState( k, q, solution[0], solution[1], solution[2] );
   
    hardeningParam = hardening( k, q, stateVariable );

    dHardeningParam_dLambda  = hardeningDerivatives( k, q, stateVariable );
    dHardeningParam_dLambda *= stateDerivatives( k, q, solution[0], solution[1], solution[2] );

    // assemble residual system
    // resid1 = P - trialP + dlambda*bulkMod*dG/dP = 0
    // resid2 = Q - trialQ + dlambda*3*shearMod*dG/dQ = 0
    // resid3 = F = 0

    residual[0] = solution[0] - trialP + solution[2] * m_bulkModulus[k] * dG[0];
    residual[1] = solution[1] - trialQ + solution[2] * 3.0 * m_shearModulus[k] * dG[1];
    residual[2] = yield( k, q, solution[0], solution[1], hardeningParam );

    // check for convergence

    norm = LvArray::tensorOps::l2Norm< 3 >( residual );

    if( iter==0 )
    {
      normZero = norm;
    }

    if( norm < 1e-8*(normZero+1))
    {
      break;
    }

    // solve Newton system

    jacobian[0][0] = 1.0 + solution[2] * m_bulkModulus[k] * dG[2];
    jacobian[0][1] = solution[2] * m_bulkModulus[k] * dG[3];
    jacobian[0][2] = m_bulkModulus[k] *dG[0] + solution[2] * m_bulkModulus[k] * dG[4] * dHardeningParam_dLambda;
    jacobian[1][0] = solution[2] * 3.0 * m_shearModulus[k] * dG[5];
    jacobian[1][1] = 1.0 + solution[2] * 3.0 * m_shearModulus[k] * dG[6];
    jacobian[1][2] = 3.0 * m_shearModulus[k] * dG[1] + solution[2] * 3.0 * m_shearModulus[k] * dG[7] * dHardeningParam_dLambda;
    jacobian[2][0] = dF[0];
    jacobian[2][1] = dF[1];
    jacobian[2][2] = dF[2] * dHardeningParam_dLambda;

    LvArray::tensorOps::invert< 3 >( jacobianInv, jacobian );
    LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( delta, jacobianInv, residual );

    for( localIndex i=0; i<3; ++i )
    {
      solution[i] -= delta[i];
    }
  }

  // re-construct stress = P*eye + sqrt(2/3)*Q*nhat

  twoInvariant::stressRecomposition( solution[0],
                                     solution[1],
                                     deviator,
                                     stress );

  // construct consistent tangent stiffness
  // note: if trialQ = 0, we will get a divide by zero error below,
  // but this is an unphysical (zero-strength) state anyway

  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );

  real64 c1 = 2 * m_shearModulus[k] * solution[1] / trialQ;
  real64 c2 = jacobianInv[0][0] * m_bulkModulus[k] - c1 / 3;
  real64 c3 = sqrt( 2./3 ) * 3 * m_shearModulus[k] * jacobianInv[0][1];
  real64 c4 = sqrt( 2./3 ) * m_bulkModulus[k] * jacobianInv[1][0];
  real64 c5 = 2 * jacobianInv[1][1] * m_shearModulus[k] - c1;

  real64 identity[6];

  for( localIndex i=0; i<3; ++i )
  {
    stiffness[i][i] = c1;
    stiffness[i+3][i+3] = 0.5 * c1;
    identity[i] = 1.0;
    identity[i+3] = 0.0;
  }

  for( localIndex i=0; i<6; ++i )
  {
    for( localIndex j=0; j<6; ++j )
    {
      stiffness[i][j] +=   c2 * identity[i] * identity[j]
                         + c3 * identity[i] * deviator[j]
                         + c4 * deviator[i] * identity[j]
                         + c5 * deviator[i] * deviator[j];
    }
  }

  // save stress, state and return
  saveStateVariable( k, q, stateVariable );
  saveStress( k, q, stress );
  return;
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerExtendedUpdates::smallStrainUpdate( localIndex const k,
                                                      localIndex const q,
                                                      real64 const ( &strainIncrement )[6],
                                                      real64 ( & stress )[6],
                                                      DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, strainIncrement, stress, stiffness.m_c );
}



/**
 * @class DruckerPragerExtended
 *
 * Extended Drucker-Prager material model.
 */
class DruckerPragerExtended : public ElasticIsotropic
{
public:

  /// @typedef Alias for DruckerPragerExtendedUpdates
  using KernelWrapper = DruckerPragerExtendedUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  DruckerPragerExtended( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~DruckerPragerExtended() override;


  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "DruckerPragerExtended";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default initial friction angle
    static constexpr auto defaultInitialFrictionAngleString = "defaultInitialFrictionAngle";

    /// string/key for default initial friction angle
    static constexpr auto defaultResidualFrictionAngleString = "defaultResidualFrictionAngle";

    /// string/key for default dilation angle
    static constexpr auto defaultDilationRatioString = "defaultDilationRatio";

    /// string/key for default hardening rate
    static constexpr auto defaultHardeningString = "defaultHardening";

    /// string/key for default cohesion
    static constexpr auto defaultCohesionString = "defaultCohesion";

    /// string/key for initial friction angle
    static constexpr auto initialFrictionString  = "initialFriction";

    /// string/key for final friction angle
    static constexpr auto residualFrictionString  = "residualFriction";

    /// string/key for dilation angle
    static constexpr auto dilationRatioString  = "dilationRatio";

    /// string/key for pressure intercept
    static constexpr auto pressureInterceptString  = "pressureIntercept";

    /// string/key for cohesion
    static constexpr auto hardeningString  = "hardening";

    /// string/key for state variable
    static constexpr auto stateString  = "stateVariable";
  };

  /**
   * @brief Create a instantiation of the DruckerPragerExtendedUpdate class that refers to the data in this.
   * @return An instantiation of DruckerPragerExtendedUpdate.
   */
  DruckerPragerExtendedUpdates createKernelUpdates() const
  {
    return DruckerPragerExtendedUpdates( m_initialFriction,
                                         m_residualFriction,
                                         m_dilationRatio,
                                         m_pressureIntercept,
                                         m_hardening,
                                         m_state,
                                         m_bulkModulus,
                                         m_shearModulus,
                                         m_newStress,
                                         m_oldStress );
  }

protected:
  virtual void postProcessInput() override;

  /// Material parameter: The default value of the initial yield surface slope
  real64 m_defaultInitialFrictionAngle;

  /// Material parameter: The default value of the final yield surface slope
  real64 m_defaultResidualFrictionAngle;

  /// Material parameter: The default value of the plastic potential slope ratio
  real64 m_defaultDilationRatio;

  /// Material parameter: The default value of the initial cohesion
  real64 m_defaultCohesion;

  /// Material parameter: The default value of the hardening rate
  real64 m_defaultHardening;

  /// Material parameter: The initial yield surface slope param for each element
  array1d< real64 > m_initialFriction;

  /// Material parameter: The final yield surface slope param for each element
  array1d< real64 > m_residualFriction;

  /// Material parameter: The plastic potential slope param (a ratio w.r.t. current yield surface)
  array1d< real64 > m_dilationRatio;

  /// Material parameter: The pressure intercept (location of cone vertex) for each element
  array1d< real64 > m_pressureIntercept;

  /// Material parameter: The hyperbolic hardening parameter for each element
  array1d< real64 > m_hardening;

  /// State variable: The current equivalent plastic shear strain for each quadrature point
  array2d< real64 > m_state;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGEREXTENDED_HPP_ */

