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
 *  @file ElastoPlasticUpdates.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_ELASTOPLASTICUPDATES_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_ELASTOPLASTICUPDATES_HPP_

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
 * @class ElastoPlasticUpdates
 *
 * Class to provide basic elastoplastic material updates
 */
class ElastoPlasticUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  ElastoPlasticUpdates( arrayView1d< real64 const > const & bulkModulus,
                        arrayView1d< real64 const > const & shearModulus,
                        arrayView3d< real64, solid::STRESS_USD > const & newStress,
                        arrayView3d< real64, solid::STRESS_USD > const & oldStress ):// TODO tmp stress[6] can be considered 
                                                                                     // in the Elasto-Plastic Newton loops
                                                                                     // to avoid holding both new and old stress 
                                                                                     // on the system
    ElasticIsotropicUpdates( bulkModulus, shearModulus, newStress, oldStress )
  {}

  /// Default copy constructor
  ElastoPlasticUpdates( ElastoPlasticUpdates const & ) = default;

  /// Default move constructor
  ElastoPlasticUpdates( ElastoPlasticUpdates && ) = default;

  /// Deleted default constructor
  ElastoPlasticUpdates() = delete;

  /// Deleted copy assignment operator
  ElastoPlasticUpdates & operator=( ElastoPlasticUpdates const & ) = delete;

  /// Deleted move assignment operator
  ElastoPlasticUpdates & operator=( ElastoPlasticUpdates && ) =  delete;

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

private:

  // The yield function that defines the elastic-plastic limit.

  virtual real64 yield( localIndex const k,
                        localIndex const q,
                        real64 const invP,
                        real64 const invQ,
                        real64 const hardeningParam ) const
  {
    GEOSX_UNUSED_VAR( k );  
    GEOSX_UNUSED_VAR( q );  
    GEOSX_UNUSED_VAR( invP );  
    GEOSX_UNUSED_VAR( invQ );  
    GEOSX_UNUSED_VAR( hardeningParam );  

    return 0;
  }

  // Derivatives of the yield function to the stress invariants and the hardening parameter.

  virtual void yieldDerivatives( localIndex const k,
                                 localIndex const q,
                                 real64 const invP,
                                 real64 const invQ,
                                 real64 const hardeningParam,
                                 real64 (& dF)[3] ) const
  {
    GEOSX_UNUSED_VAR( k );  
    GEOSX_UNUSED_VAR( q );  
    GEOSX_UNUSED_VAR( invP );  
    GEOSX_UNUSED_VAR( invQ );  
    GEOSX_UNUSED_VAR( hardeningParam );  
    GEOSX_UNUSED_VAR( dF );  
  }

  // Derivatives of the plastic potential to the stress invariants and the hardening parameter.

  virtual void potentialDerivatives( localIndex const k,
                                     localIndex const q,
                                     real64 const invP,
                                     real64 const invQ,
                                     real64 const hardeningParam,
                                     real64 (& dG)[8] ) const
  {
    GEOSX_UNUSED_VAR( k );  
    GEOSX_UNUSED_VAR( q );  
    GEOSX_UNUSED_VAR( invP );  
    GEOSX_UNUSED_VAR( invQ );  
    GEOSX_UNUSED_VAR( hardeningParam );  
    GEOSX_UNUSED_VAR( dG );  
  }

  // The hardening function that defines the relationship between the hardening parameter
  // and the state variable (ex. volumetric and/or deviatoric plastic strain).

  virtual real64 hardening( localIndex const k,
                            localIndex const q,
                            real64 const state ) const
  {
    GEOSX_UNUSED_VAR( k );  
    GEOSX_UNUSED_VAR( q );  
    GEOSX_UNUSED_VAR( state );  

    return 0;
  }

  // Derivative of the hardening parameter to the state variable.
  
  virtual real64 hardeningDerivatives( localIndex const k,
                                       localIndex const q,
                                       real64 const state ) const
  {
    GEOSX_UNUSED_VAR( k );  
    GEOSX_UNUSED_VAR( q );  
    GEOSX_UNUSED_VAR( state );  

    return 0;
  }

  // A function to update temporally the state variable in the elasto-plastic Newton loops.

  virtual real64 tmpState( localIndex const k,
                           localIndex const q,
                           real64 const invP,
                           real64 const invQ,
                           real64 const plasticMultiplier ) const
  {
    GEOSX_UNUSED_VAR( k );  
    GEOSX_UNUSED_VAR( q );  
    GEOSX_UNUSED_VAR( invP );  
    GEOSX_UNUSED_VAR( invQ );  
    GEOSX_UNUSED_VAR( plasticMultiplier );  

    return 0; 
  }

  // Derivative of the state variable to the plastic multiplier.

  virtual real64 stateDerivatives( localIndex const k,
                                   localIndex const q,
                                   real64 const invP,
                                   real64 const invQ,
                                   real64 const plasticMultiplier ) const
  {
    GEOSX_UNUSED_VAR( k );  
    GEOSX_UNUSED_VAR( q );  
    GEOSX_UNUSED_VAR( invP );  
    GEOSX_UNUSED_VAR( invQ );  
    GEOSX_UNUSED_VAR( plasticMultiplier );  

    return 0; 
  }

  // Get the state variable that was saved from the previous elasto-plastic update.

  virtual real64 getStateVariable( localIndex const k,
                                   localIndex const q ) const
  {
    GEOSX_UNUSED_VAR( k );  
    GEOSX_UNUSED_VAR( q );  

    return 0;
  }

  // Save the converged state variable after the elasto-plastic Newton loops.
 
  virtual void saveStateVariable( localIndex const k,
                                  localIndex const q,
                                  real64 const state ) const
  {
    GEOSX_UNUSED_VAR( k );  
    GEOSX_UNUSED_VAR( q );  
    GEOSX_UNUSED_VAR( state );  
  }

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ElastoPlasticUpdates::smallStrainUpdate( localIndex const k,
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
void ElastoPlasticUpdates::smallStrainUpdate( localIndex const k,
                                              localIndex const q,
                                              real64 const ( &strainIncrement )[6],
                                              real64 ( & stress )[6],
                                              DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, strainIncrement, stress, stiffness.m_c );
}

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_ELASTOPLASTICUPDATES_HPP_ */
