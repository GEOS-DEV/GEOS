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
 *  @file DruckerPrager.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGER_HPP
#define GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGER_HPP

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class DruckerPragerUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class DruckerPragerUpdates : public ElasticIsotropicUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] stress The ArrayView holding the stress data for each quadrature point.
   */
  DruckerPragerUpdates( arrayView1d< real64 const > const & tanFrictionAngle,
                        arrayView1d< real64 const > const & tanDilationAngle,
                        arrayView1d< real64 const > const & hardeningRate,
                        arrayView2d< real64 > const & newCohesion,
                        arrayView2d< real64 > const & oldCohesion,
                        arrayView1d< real64 const > const & bulkModulus,
                        arrayView1d< real64 const > const & shearModulus,
                        arrayView3d< real64, solid::STRESS_USD > const & newStress,
                        arrayView3d< real64, solid::STRESS_USD > const & oldStress ):
    ElasticIsotropicUpdates( bulkModulus, shearModulus, newStress,oldStress ),
    m_tanFrictionAngle( tanFrictionAngle ),
    m_tanDilationAngle( tanDilationAngle ),
    m_hardeningRate( hardeningRate ),
    m_newCohesion( newCohesion ),
    m_oldCohesion( oldCohesion )
  {}

  /// Default copy constructor
  DruckerPragerUpdates( DruckerPragerUpdates const & ) = default;

  /// Default move constructor
  DruckerPragerUpdates( DruckerPragerUpdates && ) = default;

  /// Deleted default constructor
  DruckerPragerUpdates() = delete;

  /// Deleted copy assignment operator
  DruckerPragerUpdates & operator=( DruckerPragerUpdates const & ) = delete;

  /// Deleted move assignment operator
  DruckerPragerUpdates & operator=( DruckerPragerUpdates && ) =  delete;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( & strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  real64 ( & stiffness )[6][6] ) override final;
  
  GEOSX_HOST_DEVICE
  virtual void saveConvergedState() override final;
    
private:
  /// A reference to the ArrayView holding the friction angle for each element.
  arrayView1d< real64 const > const m_tanFrictionAngle;
  
  /// A reference to the ArrayView holding the dilation angle for each element.
  arrayView1d< real64 const > const m_tanDilationAngle;
  
  /// A reference to the ArrayView holding the hardening rate for each element.
  arrayView1d< real64 const > const m_hardeningRate;
  
  /// A reference to the ArrayView holding the new cohesion for each integration point
  arrayView2d< real64 > const m_newCohesion;
  
  /// A reference to the ArrayView holding the old cohesion for each integration point
  arrayView2d< real64 > const m_oldCohesion;
  
};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::smallStrainUpdate( localIndex const k,
                                              localIndex const q,
                                              real64 const ( & strainIncrement )[6],
                                              real64 ( & stress )[6],
                                              real64 ( & stiffness )[6][6] )
{
  // elastic predictor (assume strainIncrement is all elastic)
  
  ElasticIsotropicUpdates::smallStrainUpdate( k, q, strainIncrement, stress, stiffness);
  
  // decompose into mean (P) and von mises (Q) stress invariants
  // could switch to strain invariant formulation using getElasticStrain() if needed
    
  real64 trialP;
  real64 trialQ;
  real64 deviator[6];
  
  twoInvariant::stressDecomposition(stress,
                                    trialP,
                                    trialQ,
                                    deviator);
  
  // check yield function F <= 0, using old hardening variable state
    
  real64 yield = trialQ + m_tanFrictionAngle[k] * trialP - m_oldCohesion[k][q];
  
  if(yield < 1e-9) // elasticity
  {
    return;
  }
  
  // else, plasticity (trial stress point lies outside yield surface)

  // the return mapping can in general be written as a newton iteration.
  // here we have a linear problem, so the algorithm will converge in one
  // iteration, but this is a template for more general models with either
  // nonlinear hardening or yield surfaces.
    
  real64 solution[3], residual[3], delta[3];
  real64 jacobian[3][3] = {{}}, jacobianInv[3][3];
  
  solution[0] = trialP; // initial guess for newP
  solution[1] = trialQ; // initial guess for newQ
  solution[2] = 0;      // initial guess for plastic multiplier
  
  real64 norm,normZero = 1e30;
  
  real64 const & shear       = m_shearModulus[k];
  real64 const & bulk        = m_bulkModulus[k];
  real64 const & friction    = m_tanFrictionAngle[k];
  real64 const & dilation    = m_tanDilationAngle[k];
  real64 const & hardening   = m_hardeningRate[k];
  real64 const & oldCohesion = m_oldCohesion[k][q];
  real64       & newCohesion = m_newCohesion[k][q];
  
  // begin newton loop
  
  for(localIndex iter=0; iter<20; ++iter)
  {
    // apply a linear cohesion decay model.
    
    newCohesion = oldCohesion-solution[2]*hardening;
    real64 cohesionDeriv = -hardening;
    
    if(newCohesion < 0)  // check for complete cohesion loss
    {
      newCohesion = 0;
      cohesionDeriv = 0;
    }
    
    // assemble residual system
    
    residual[0] = solution[0] - trialP + solution[2]*bulk*dilation; // P - trialP + lambda*dG/dP = 0
    residual[1] = solution[1] - trialQ + solution[2]*3*shear;       // Q - trialQ + lambda*dG/dQ = 0
    residual[2] = solution[1] + friction*solution[0] - newCohesion; // F = 0
    
    // check for convergence
    
    norm = LvArray::tensorOps::l2Norm<3>(residual);
        
    if(iter==0)
    {
      normZero = norm;
    }
    
    if(norm < 1e-8*(normZero+1))
    {
      break;
    }
    
    // solve Newton system
    
    jacobian[0][0] = 1;
    jacobian[0][2] = bulk*dilation;
    jacobian[1][1] = 1;
    jacobian[1][2] = 3*shear;
    jacobian[2][0] = friction;
    jacobian[2][1] = 1;
    jacobian[2][2] = cohesionDeriv;
    
    LvArray::tensorOps::invert<3>(jacobianInv,jacobian);
    LvArray::tensorOps::AijBj<3,3>(delta,jacobianInv,residual);
   
    for(localIndex i=0; i<3; ++i)
    {
      solution[i] -= delta[i];
    }
  }
  
  // construct stress = P*eye + sqrt(2/3)*Q*nhat
  
  twoInvariant::stressRecomposition(solution[0],
                                    solution[1],
                                    deviator,
                                    stress);
                                        
  // construct consistent tangent operator
  
  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );
  
  real64 c1 = 2*shear*solution[1]/trialQ;   // divide by zero possible, but only in unphysical zero-strength state
  real64 c2 = jacobianInv[0][0]*bulk - c1/3;
  real64 c3 = sqrt(2./3)*3*shear*jacobianInv[0][1];
  real64 c4 = sqrt(2./3)*bulk*jacobianInv[1][0];
  real64 c5 = 2*jacobianInv[1][1]*shear - c1;
  
  real64 identity[6];
  
  for(localIndex i=0; i<3; ++i)
  {
    stiffness[i][i] = c1;
    stiffness[i+3][i+3] = 0.5*c1;
    identity[i] = 1.0;
    identity[i+3] = 0.0;
  }
  
  for(localIndex i=0; i<6; ++i)
  {
    for(localIndex j=0; j<6; ++j)
    {
      stiffness[i][j] +=   c2 * identity[i] * identity[j]
                         + c3 * identity[i] * deviator[j]
                         + c4 * deviator[i] * identity[j]
                         + c5 * deviator[i] * deviator[j];
    }
  }
  
  // save new stress and return
  
  saveStress( k, q, stress );
  return;
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::saveConvergedState()
{
  ElasticIsotropicUpdates::saveConvergedState();
  m_oldCohesion.setValues< serialPolicy >( m_newCohesion );
}


/**
 * @class DruckerPrager
 *
 * Drucker-Prager material model.
 */
class DruckerPrager : public ElasticIsotropic
{
public:

  /// @typedef Alias for DruckerPragerUpdates
  using KernelWrapper = DruckerPragerUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  DruckerPrager( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~DruckerPrager() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;
  
  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "DruckerPrager";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string CatalogName() { return m_catalogNameString; }

  virtual string GetCatalogName() override { return CatalogName(); }

  ///@}

  /**
   * @struct Set of "char const *" and keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default friction angle
    static constexpr auto defaultTanFrictionAngleString = "defaultTanFrictionAngle";
    
    /// string/key for default dilation angle
    static constexpr auto defaultTanDilationAngleString = "defaultTanDilationAngle";
    
    /// string/key for default hardening rate
    static constexpr auto defaultHardeningRateString = "defaultHardeningRate";
    
    /// string/key for default cohesion
    static constexpr auto defaultCohesionString = "defaultCohesion";
        
    /// string/key for friction angle
    static constexpr auto tanFrictionAngleString  = "TanFrictionAngle";
    
    /// string/key for dilation angle
    static constexpr auto tanDilationAngleString  = "TanDilationAngle";
    
    /// string/key for cohesion
    static constexpr auto hardeningRateString  = "HardeningRate";
    
    /// string/key for cohesion
    static constexpr auto newCohesionString  = "NewCohesion";
    
        /// string/key for cohesion
    static constexpr auto oldCohesionString  = "OldCohesion";
  };

  /**
   * @brief Create a instantiation of the DruckerPragerUpdate class that refers to the data in this.
   * @return An instantiation of DruckerPragerUpdate.
   */
  DruckerPragerUpdates createKernelWrapper() const
  {
    return DruckerPragerUpdates( m_tanFrictionAngle,
                                 m_tanDilationAngle,
                                 m_hardeningRate,
                                 m_newCohesion,
                                 m_oldCohesion,
                                 m_bulkModulus,
                                 m_shearModulus,
                                 m_newStress,
                                 m_oldStress );
  }

protected:
  virtual void PostProcessInput() override;
    
  /// Material parameter: The default value of yield surface slope
  real64 m_defaultTanFrictionAngle;
  
    /// Material parameter: The default value of plastic potential slope
  real64 m_defaultTanDilationAngle;
  
  /// Material parameter: The default value of the initial cohesion
  real64 m_defaultCohesion;
  
  /// Material parameter: The default value of the hardening rate
  real64 m_defaultHardeningRate;
    
  /// Material parameter: The yield surface slope for each element
  array1d< real64 > m_tanFrictionAngle;
  
  /// Material parameter: The plastic potential slope for each element
  array1d< real64 > m_tanDilationAngle;
  
  /// Material parameter: The hardening rate each element
  array1d< real64 > m_hardeningRate;
  
  /// History variable: The current cohesion value for each quadrature point
  array2d< real64 > m_newCohesion;
  
  /// History variable: The previous cohesion value for each quadrature point
  array2d< real64 > m_oldCohesion;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGER_HPP_ */


