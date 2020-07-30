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
 *  @file DelftEgg.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DELFTEGG_HPP
#define GEOSX_CONSTITUTIVE_SOLID_DELFTEGG_HPP

#include "SolidBase.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace constitutive
{

/**
 * @class DelftEggUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class DelftEggUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] stress The ArrayView holding the stress data for each quadrature point.
   */
  DelftEggUpdates( arrayView1d< real64 const > const & bulkModulus,
                        arrayView1d< real64 const > const & shearModulus,
                        arrayView2d< real64 const > const & refPressure,
                        arrayView2d< real64 const > const & refStrainVol,
                        arrayView1d< real64 const > const & recompressionIndex,
                        arrayView1d< real64 const > const & virginCompressionIndex,
                        arrayView1d< real64 const > const & cslSlope,
                        arrayView1d< real64 const > const & shapeParameter,
                        arrayView2d< real64 > const & newPreConsolidationPressure,
                        arrayView2d< real64 > const & oldPreConsolidationPressure,
                        arrayView3d< real64, solid::STRESS_USD > const & newStress,
                        arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                        arrayView3d< real64, solid::STRESS_USD > const & stress):
    SolidBaseUpdates( stress ),
    m_bulkModulus( bulkModulus ),
    m_shearModulus( shearModulus ),
    m_refPressure( refPressure ),
    m_refStrainVol( refStrainVol ),
    m_recompressionIndex( recompressionIndex ),
    m_virginCompressionIndex( virginCompressionIndex ),
    m_cslSlope( cslSlope ),
    m_shapeParameter( shapeParameter ),
    m_newPreConsolidationPressure( newPreConsolidationPressure ),
    m_oldPreConsolidationPresure( oldPreConsolidationPressure ),
    m_newStress( newStress ),
    m_oldStress( oldStress )
  {}

  /// Default copy constructor
  DelftEggUpdates( DelftEggUpdates const & ) = default;

  /// Default move constructor
  DelftEggUpdates( DelftEggUpdates && ) = default;

  /// Deleted default constructor
  DelftEggUpdates() = delete;

  /// Deleted copy assignment operator
  DelftEggUpdates & operator=( DelftEggUpdates const & ) = delete;

  /// Deleted move assignment operator
  DelftEggUpdates & operator=( DelftEggUpdates && ) =  delete;

  GEOSX_HOST_DEVICE
  virtual void SmallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  arraySlice1d< real64 const > const & strainIncrement,
                                  arraySlice1d< real64 > const & stress,
                                  arraySlice2d< real64 > const & stiffness ) override final;
  
  GEOSX_HOST_DEVICE
  virtual void SaveConvergedState( localIndex const k,
                                   localIndex const q ) override final;
    
private:
  /// A reference to the ArrayView holding the bulk modulus for each element.
  arrayView1d< real64 const > const m_bulkModulus;

  /// A reference to the ArrayView holding the shear modulus for each element.
  arrayView1d< real64 const > const m_shearModulus;
  
  /// A reference to the ArrayView holding the reference pressure (p0) for each element.
  arrayView2d< real64 const > const m_refPressure;

    /// A reference to the ArrayView holding the shear modulus for each element.
  arrayView2d< real64 const > const m_refStrainVol;

  /// A reference to the ArrayView holding the recompression index for each element.
  arrayView1d< real64 const > const m_recompressionIndex;
  
  /// A reference to the ArrayView holding the virgin compression index for each element.
  arrayView1d< real64 const > const m_virginCompressionIndex;
  
  /// A reference to the ArrayView holding the slope of the critical state line for each element.
  arrayView1d< real64 const > const m_cslSlope;
  
  /// A reference to the ArrayView holding the shape parameter for each integration point
  arrayView2d< real64 const > const m_shapeParameter;

  /// A reference to the ArrayView holding the old preconsolidation pressure for each integration point
  arrayView2d< real64 > const m_oldPreConsolidationPressure;

  /// A reference to the ArrayView holding the new preconsolidation pressure for each integration point
  arrayView2d< real64 > const m_newPreConsolidationPressure;
  
  /// A reference to the ArrayView holding the new stress for each integration point
  arrayView3d< real64, solid::STRESS_USD > const m_newStress;
  
  /// A reference to the ArrayView holding the old stress for each integration point
  arrayView3d< real64, solid::STRESS_USD > const m_oldStress;
};

///////////////////////// new proposal /////////////////////////////////////////

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DelftEggUpdates::SmallStrainUpdate( localIndex const k,
                                              localIndex const q,
                                              arraySlice1d< real64 const > const & strainIncrement,
                                              arraySlice1d< real64 > const & stress,
                                              arraySlice2d< real64 > const & stiffness )
{
  real64 const pc     = m_oldPreConsolidationPressure[k][q]; //pre-consolidation pressure
  real64 const mu     = m_shearModulus[k];
  real64 const p0     = m_refPressure[k][q];
  real64 const eps_v0 = m_refStrainVol[k][q];
  
  real64 const M      = m_cslSlope[k];
  real64 const Cr     = m_recompressionIndex[k];
  real64 const Cc     = m_virginCompressionIndex[k];
  real64 const alpha  = m_shapeParameter[k];

  // two-invariant decomposition of strain increment

  real64 strainIncrementVol = 0; //volumetric part of strain increment 
  for(localIndex i=0; i<3; ++i)
  {
    strainIncrementVol += strainIncrement[i];
  }

  real64 temp = strainIncrementVol/3
  array1d< real64 > deviator(6);  // array allocation
  for(localIndex i=0; i<3; ++i)
  {
    deviator[i] = strainIncrement[i]-temp;
    deviator[i+3] = strainIncrement[i+3];
  }

  real64 strainIncrementDev = 0; //deviatoric part of strain increment
  for(localIndex i=0; i<3; ++i)
  {
    strainIncrementDev += deviator[i]*deviator[i];
    strainIncrementDev += deviator[i+3]*deviator[i+3]/2; //because of Voigt form of strain increment vector
  }
  strainIncrementDev = strainIncrementDev + 1e-15; // (NOT needed for now?) perturbed to avoid divide by zero when strainIncrementDev=0;

  strainIncrementDev *= std::sqrt(2./3.);

  // two-invariant decomposition of old stress in P-Q space (mean & deviatoric stress)

  real64 oldP = 0;
  real64 oldQ = 0;
  for(localIndex i=0; i<6; ++i)
  {
    stress[i] = m_oldStress[k][q][i];
  }

  for(localIndex i=0; i<3; ++i)
  {
    oldP += stress[i];
  }
  oldP /= 3;
  
  array1d< real64 > oldDeviator(6);  // array allocation
  for(localIndex i=0; i<3; ++i)
  {
    oldDeviator[i] = stress[i]-oldP;
    oldDeviator[i+3] = stress[i+3];
  }
  
  for(localIndex i=0; i<3; ++i)
  {
    oldQ += oldDeviator[i]*oldDeviator[i];
    oldQ += 2*oldDeviator[i+3]*oldDeviator[i+3];
  }


  oldQ = std::sqrt(oldQ) + 1e-15; // perturbed to avoid divide by zero when Q=0;

    for(localIndex i=0; i<6; ++i)
  {
    oldDeviator[i] /= oldQ; // normalized deviatoric direction, "nhat" from previous step
  }

  oldQ *= std::sqrt(3./2.);

  // Recover elastic strains from the previous step, based on stress from the previous step
  // [Note: in order to minimize data transfer, we are not storing and passing elastic strains] 

  real64 oldElasticStrainVol = std::log(oldP/p0) * Cr * (-1) + eps_v0 ;
  real64 oldElasticStrainDev = oldQ/3/mu;  

  // Now recover the old strain tensor from the strain invariants. 
  // Note that we need the deviatoric direction (n-hat) from the previous step.

  
  
  // elastic predictor
  // newP= oldP * exp(-1/Cr* strainIncrementVol)
  // newQ = oldQ + 3 * mu * strainIncrementDev

  real64 eps_v_trial = oldElasticStrainVol + strainIncrementVol; 
  real64 eps_s_trial = oldElasticStrainDev + strainIncrementDev;

  real64 trialP = oldP * std::exp(-1/Cr* strainIncrementVol);
  real64 trialQ = oldQ + 3 * mu * strainIncrementDev;

  // Calculate the normalized deviatoric direction, "nhat"

  for(localIndex i=0; i<6; ++i)
  {
    stress[i] = m_oldStress[k][q][i];
    
    for(localIndex j=0; j<6; ++j)
    {
      stress[i] += stiffness[i][j]*strainIncrement[j];
    }
  }

  // two-invariant decomposition in P-Q space (mean & deviatoric stress)
  real64 trialP = 0;
  for(localIndex i=0; i<3; ++i)
  {
    trialP += stress[i];
  }
  trialP /= 3;
  
  array1d< real64 > deviator(6);  // array allocation
  for(localIndex i=0; i<3; ++i)
  {
    deviator[i] = stress[i]-trialP;
    deviator[i+3] = stress[i+3];
  }
  
  real64 trialQ = 0;
  for(localIndex i=0; i<3; ++i)
  {
    trialQ += deviator[i]*deviator[i];
    trialQ += deviator[i+3]*deviator[i+3];
  }
  trialQ = std::sqrt(trialQ) + 1e-15; // perturbed to avoid divide by zero when Q=0;
  
  for(localIndex i=0; i<6; ++i)
  {
    deviator[i] /= trialQ; // normalized deviatoric direction, "nhat"
  }
  trialQ *= std::sqrt(3./2.);
  // check yield function F <= 0
  
  real64 yield = trialQ + friction*trialP - cohesion;
  
    // set stiffness to elastic predictor
  
  for(localIndex i=0; i<6; ++i)
  {
    for(localIndex j=0; j<6; ++j)
    {
      stiffness[i][j] = 0;
    }
  }
  
  stiffness[0][0] = lame + 2*shear;
  stiffness[0][1] = lame;
  stiffness[0][2] = lame;

  stiffness[1][0] = lame;
  stiffness[1][1] = lame + 2*shear;
  stiffness[1][2] = lame;

  stiffness[2][0] = lame;
  stiffness[2][1] = lame;
  stiffness[2][2] = lame + 2*shear;

  stiffness[3][3] = shear;
  stiffness[4][4] = shear;
  stiffness[5][5] = shear;

  if(yield > 1e-9) // plasticity branch
  {
    // the return mapping can in general be written as a newton iteration.
    // here we have a linear problem, so the algorithm will converge in one
    // iteration, but this is a template for more general models with either
    // nonlinear hardening or curved yield surfaces.
    
    // for GPU we can simplify the DP model to avoid Newton, but Cam-Clay and
    // others will need it.
    
    array1d< real64 > solution(3), residual(3), delta(3);
    array2d< real64 > jacobian(3,3), jacobianInv(3,3);
    
    solution[0] = trialP; // initial guess for newP
    solution[1] = trialQ; // initial guess for newQ
    solution[2] = 0;      // initial guess for plastic multiplier
    
    real64 norm,normZero = 1e30;
    jacobian = 0;
    
    for(localIndex iter=0; iter<10; ++iter) // could be fixed at one iter
    {
      // apply a linear cohesion decay model.  this requires an if() check for negative cohesion.
      // a log decay model would avoid this, at the price of nonlinearity
      
      cohesion = oldCohesion-solution[2]*hardening;
      real64 cohesionDeriv = -hardening;
      
      if(cohesion < 0)  // branch
      {
        cohesion = 0;
        cohesionDeriv = 0;
      }
      
      // assemble residual system
      
      residual[0] = solution[0] - trialP + solution[2]*bulk*dilation; // P - trialP + lambda*dG/dP = 0
      residual[1] = solution[1] - trialQ + solution[2]*3*shear;       // Q - trialQ + lambda*dG/dQ = 0
      residual[2] = solution[1] + friction*solution[0] - cohesion;    // F = 0
      
      // check for convergence (can be avoided for linear model)
      
      norm = LvArray::tensorOps::l2Norm<3>(residual);
      
      //residual.L2_Norm();  //std::cout << iter << " " << norm << std::endl;
      
      if(iter==0)
      {
        normZero = norm;
      }
      
      if(norm < 1e-8*(normZero+1)) 
      {
        break;
      }
      
      // solve Newton system
      
      jacobian(0,0) = 1;
      jacobian(0,2) = bulk*dilation;
      jacobian(1,1) = 1;
      jacobian(1,2) = 3*shear;
      jacobian(2,0) = friction;
      jacobian(2,1) = 1;
      jacobian(2,2) = cohesionDeriv;
      
      LvArray::tensorOps::invert<3>(jacobianInv,jacobian);
      LvArray::tensorOps::AijBj<3,3>(delta,jacobianInv,residual);
     
      for(localIndex i=0; i<3; ++i)
      {
        solution[i] -= delta[i];
      }
    }
    
    // construct stress = P*eye + sqrt(2/3)*Q*nhat
    
    for(localIndex i=0; i<6; ++i)
    {
      stress[i] = sqrt(2./3.)*solution[1]*deviator[i];
    }
    for(localIndex i=0; i<3; ++i)
    {
      stress[i] += solution[0];
    }
    
    // construct consistent tangent operator
    
    for(localIndex i=0; i<6; ++i)
    {
      for(localIndex j=0; j<6; ++j)
      {
        stiffness[i][j] = 0;
      }
    }
    
    real64 c1 = 2*shear*solution(1)/trialQ; // TODO: confirm trialQ != 0 for linear DP
    real64 c2 = jacobianInv(0,0)*bulk - c1/3;
    real64 c3 = sqrt(2./3)*3*shear*jacobianInv(0,1);
    real64 c4 = sqrt(2./3)*bulk*jacobianInv(1,0);
    real64 c5 = 2*jacobianInv(1,1)*shear - c1;
    
    array1d< real64 > identity(6);
    
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
    
  } // end plastic branch
  
  // remember history variables before returning
  
  m_newCohesion[k][q] = cohesion;
  
  for(localIndex i=0; i<6; ++i)
  {
    m_newStress[k][q][i] = stress[i];
  }
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DelftEggUpdates::SaveConvergedState( localIndex const k,
                                               localIndex const q )
{
  m_oldCohesion[k][q] = m_newCohesion[k][q];
  
  for(localIndex i=0; i<6; ++i)
  {
    m_oldStress[k][q][i] = m_newStress[k][q][i];
  }
}


/**
 * @class DelftEgg
 *
 * Drucker-Prager material model.
 */
class DelftEgg : public SolidBase
{
public:

  /// @typedef Alias for DelftEggUpdates
  using KernelWrapper = DelftEggUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  DelftEgg( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~DelftEgg() override;

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
  static constexpr auto m_catalogNameString = "DelftEgg";

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
    /// string/key for default bulk modulus
    static constexpr auto defaultBulkModulusString  = "defaultBulkModulus";

    /// string/key for default shear modulus
    static constexpr auto defaultShearModulusString = "defaultShearModulus";

    /// string/key for default friction angle
    static constexpr auto defaultTanFrictionAngleString = "defaultTanFrictionAngle";
    
    /// string/key for default dilation angle
    static constexpr auto defaultTanDilationAngleString = "defaultTanDilationAngle";
    
    /// string/key for default hardening rate
    static constexpr auto defaultHardeningRateString = "defaultHardeningRate";
    
    /// string/key for default cohesion
    static constexpr auto defaultCohesionString = "defaultCohesion";
    
    /// string/key for bulk modulus
    static constexpr auto bulkModulusString  = "BulkModulus";
    
    /// string/key for poisson ratio
    static constexpr auto shearModulusString = "ShearModulus";
    
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
    
        /// string/key for cohesion
    static constexpr auto newStressString  = "NewStress";
    
        /// string/key for cohesion
    static constexpr auto oldStressString  = "OldStress";
  };

  /**
   * @brief Create a instantiation of the DelftEggUpdate class that refers to the data in this.
   * @return An instantiation of DelftEggUpdate.
   */
  DelftEggUpdates createKernelWrapper()
  {
    return DelftEggUpdates( m_bulkModulus,
                                 m_shearModulus,
                                 m_tanFrictionAngle,
                                 m_tanDilationAngle,
                                 m_hardeningRate,
                                 m_newCohesion,
                                 m_oldCohesion,
                                 m_newStress,
                                 m_oldStress,
                                 m_stress ); // TODO: m_stress is redundant with m_newStress
  }

protected:
  virtual void PostProcessInput() override;

private:
  //TODO: maybe define some helper structs for defaults, arrays, arrayViews, etc.
  
  /// Material parameter: The default value of the bulk modulus
  real64 m_defaultBulkModulus;

  /// Material parameter: The default value of the shear modulus
  real64 m_defaultShearModulus;
  
  /// Material parameter: The default value of yield surface slope
  real64 m_defaultTanFrictionAngle;
  
    /// Material parameter: The default value of plastic potential slope
  real64 m_defaultTanDilationAngle;
  
  /// Material parameter: The default value of the initial cohesion
  real64 m_defaultCohesion;
  
  /// Material parameter: The default value of the hardening rate
  real64 m_defaultHardeningRate;
  
  /// Material parameter: The bulk modulus for each element
  array1d< real64 > m_bulkModulus;

  /// Material parameter: The shear modulus for each element
  array1d< real64 > m_shearModulus;
  
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
  
  /// History variable: The current stress for each quadrature point
  array3d< real64, solid::STRESS_PERMUTATION > m_newStress; //TODO: redundant storage with base m_stress
  
  /// History variable: The previous stress for each quadrature point.
  array3d< real64, solid::STRESS_PERMUTATION > m_oldStress;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DELFTEGG_HPP_ */


