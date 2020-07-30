/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior
 * University Copyright (c) 2018-2019 Total, S.A Copyright (c) 2019-     GEOSX
 * Contributors All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS
 * files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file ModifiedCamClay.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_MODIFIEDCAMCLAY_HPP
#define GEOSX_CONSTITUTIVE_SOLID_MODIFIEDCAMCLAY_HPP

#include "LvArray/src/tensorOps.hpp"
#include "SolidBase.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"

namespace geosx {

namespace constitutive {

/**
 * @class ModifiedCamClayUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class ModifiedCamClayUpdates : public SolidBaseUpdates 
{
public:
  /** Check how to modifiy this?
   * @brief Constructor
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] stress The ArrayView holding the stress data for each quadrature point.
   */
  ModifiedCamClayUpdates( arrayView2d<real64 const> const & refPInvariant,
                          arrayView2d<real64 const> const & refStrainVol,
                          arrayView2d<real64 const> const & refShearModulus,
                          arrayView1d<real64 const> const & alpha,
                          arrayView1d<real64 const> const & Cc,
                          arrayView1d<real64 const> const & Cr,
                          arrayView1d<real64 const> const & M,
                          arrayView1d<real64 const> const & associativity,
                          arrayView2d<real64> const & newPreconsPressure,
                          arrayView2d<real64> const & oldPreconsPressure,
                          arrayView3d<real64, solid::STRESS_USD> const & newElasticStrain, 
                          arrayView3d<real64, solid::STRESS_USD> const & oldElasticStrain, 
                          arrayView3d<real64, solid::STRESS_USD> const & stress):
    SolidBaseUpdates( stress ), 
    m_referencePInvariant( refPInvariant ),
    m_refElasticStrainVolumetric( refStrainVol ), 
    m_referenceShearModulus( refShearModulus ),
    m_shearModulusEvolution( alpha ), 
    m_virginCompressionIndex( Cc ),
    m_recompressionIndex( Cr ), 
    m_criticalStateSlope( M ),
    m_associativity( associativity ),
    m_newPreconsolidationPressure( newPreconsPressure ),
    m_oldPreconsolidationPressure( oldPreconsPressure ),
    m_newElasticStrain( newElasticStrain ),
    m_oldElasticStrain( oldElasticStrain ) 
  {}

  /// Default copy constructor
  ModifiedCamClayUpdates( ModifiedCamClayUpdates const & ) = default;

  /// Default move constructor
  ModifiedCamClayUpdates( ModifiedCamClayUpdates && ) = default;

  /// Deleted default constructor
  ModifiedCamClayUpdates() = delete;

  /// Deleted copy assignment operator
  ModifiedCamClayUpdates &operator=( ModifiedCamClayUpdates const & ) = delete;

  /// Deleted move assignment operator
  ModifiedCamClayUpdates &operator=( ModifiedCamClayUpdates && ) = delete;

  GEOSX_HOST_DEVICE
  virtual void SmallStrainUpdate( localIndex const k, 
                                  localIndex const q,
                                  arraySlice1d< real64 const > const & strainIncrement,
                                  arraySlice1d< real64 > const & stress,
                                  arraySlice2d< real64 > const & stiffness) override final;

  GEOSX_HOST_DEVICE
  virtual void SaveConvergedState( localIndex const k,
                                   localIndex const q ) override final;
  GEOSX_HOST_DEVICE
  virtual void DeleteElasticStrain( localIndex const k,
                                    localIndex const q );

private:
  /// A reference to the ArrayView holding the reference p invariant for each integration point.
  arrayView2d< real64 const > const m_referencePInvariant;

  /// A reference to the ArrayView holding the reference for the volumetric
  /// strain invariant of the elastic strain for each integration point.
  arrayView2d< real64 const > const m_refElasticStrainVolumetric;

  /// A reference to the ArrayView holding the reference shear modulus for each integration point.
  arrayView2d< real64 const > const m_referenceShearModulus;

  /// A reference to the ArrayView holding the parameter for the shear modulus evolution for 
  /// each element.
  arrayView1d< real64 const > const m_shearModulusEvolution;

  /// A reference to the ArrayView holding the virgin compression index Cc for each element.
  arrayView1d< real64 const > const m_virginCompressionIndex;

  /// A reference to the ArrayView holding the recompression index Cr for each element.
  arrayView1d< real64 const > const m_recompressionIndex;

  /// A reference to the ArrayView holding the slope M for the critical state line for  
  /// each element.
  arrayView1d< real64 const > const m_criticalStateSlope;

  /// A reference to the ArrayView holding the associativity parameter (to accomodate 
  /// non-associated plasticity) for each element
  arrayView1d< real64 const > const m_associativity;

  /// A reference to the ArrayView holding the new preconsolidation pressure for each 
  /// integration point
  arrayView2d< real64 > const m_newPreconsolidationPressure;

  /// A reference to the ArrayView holding the old preconsolidation pressure for each 
  /// integration point
  arrayView2d< real64 > const m_oldPreconsolidationPressure;

  /// A reference to the ArrayView holding the new strain for each integration point
  arrayView3d< real64, solid::STRESS_USD > const m_newElasticStrain;

  /// A reference to the ArrayView holding the old strain for each integration point
  arrayView3d< real64, solid::STRESS_USD > const m_oldElasticStrain;
};

////////////////////////////////////////////////////////////////////////////////////

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ModifiedCamClayUpdates::SmallStrainUpdate( localIndex const k, 
                                              localIndex const q,
                                              arraySlice1d< real64 const > const & strainIncrement,
                                              arraySlice1d< real64 > const & stress,
                                              arraySlice2d< real64 > const & stiffness ) 
{
  real64 const refPInvariant = m_referencePInvariant[k][q];
  real64 const refStrainVol = m_refElasticStrainVolumetric[k][q];
  real64 const refShearModulus = m_referenceShearModulus[k][q];
  real64 const alpha = m_shearModulusEvolution[k];

  real64 const oldPreconsPressure = m_oldPreconsolidationPressure[k][q];
  real64 const Cc = m_virginCompressionIndex[k];
  real64 const Cr = m_recompressionIndex[k];
  real64 const M = m_criticalStateSlope[k];
  real64 const associativity = m_associativity[k];

  real64 preconsPressure = oldPreconsPressure;

  // elastic predictor
  // newStrainElastic = oldStrain + strainIncrement
  array1d< real64 > strainElastic(6);
  
  for (localIndex i = 0; i < 6; ++i) 
  {
    strainElastic[i] = m_oldElasticStrain[k][q][i] + strainIncrement[i];
  }
  
  // Decompose the strain tensor into the elastic strain invariants strainVol and strainDev
  real64 strainVol = 0;
  for (localIndex i = 0; i < 3; ++i) 
  {
    strainVol += strainElastic[i];
  }

  array1d< real64 > deviator(6);
  for (localIndex i = 0; i < 3; ++i) 
  {
    deviator[i] = strainElastic[i] - 1. / 3. * strainVol;
    deviator[i + 3] = strainElastic[i + 3] / 2.; //divided by 2. due to Voigt notation
  }

  real64 strainDev = 0;
  for (localIndex i = 0; i < 3; ++i) 
  {
    strainDev += deviator[i] * deviator[i];
    strainDev += 2 * deviator[i + 3] * deviator[i + 3];
  }
  strainDev = std::sqrt(strainDev) + 1e-15; // perturbed to avoid divide by zero when strainDev=0;
    
  for (localIndex i = 0; i < 6; ++i) 
  {
    deviator[i] /= strainDev; // normalized deviatoric direction, "nhat"
  }
  strainDev *= std::sqrt(2. / 3.);
    
  // two-invariant decomposition in P-Q space (mean & deviatoric stress)
  real64 expOmega = std::exp((refStrainVol - strainVol) / Cr);
  real64 shear = refShearModulus - alpha * refPInvariant * expOmega;
  real64 Q = 3 * shear * strainDev;
  real64 P = refPInvariant * (1 + 3 * alpha * strainDev * strainDev / (2 * Cr)) * expOmega;
  real64 bulk = -P / Cr;

  // Q: should we add a check for Poisson ration<0?
  // real64 poissonRatio = ( 3 * bulk - 2 * shear) / ( 2 * ( 3 * bulk + shear ) );
  // GEOSX_ASSERT_MSG(poissonRatio >=0, "Negative poisson ratio produced");

  // construct stress = P*eye + sqrt(2/3)*Q*nhat

  for (localIndex i = 0; i < 6; ++i) 
  {
    stress[i] = sqrt(2. / 3.) * Q * deviator[i];
  }
  for (localIndex i = 0; i < 3; ++i) 
  {
    stress[i] += P;
  }

  // set stiffness to elastic predictor
  // stiffness = c1*I + c2*eye_dyad_eye + c3*(eye_dyad_nhat+nhat_dyad_eye)

  // constants
  real64 c1 = 2 * shear;
  real64 c2 = bulk - c1 / 3. ;
  real64 coupling = 3 * refPInvariant * alpha * strainDev * expOmega / Cr;
  real64 c3 = std::sqrt(2. / 3.) * coupling;

  for (localIndex i = 0; i < 6; ++i) 
  {
    for (localIndex j = 0; j < 6; ++j) 
    {
      stiffness[i][j] = 0;
    }
  }
  
  array1d<real64> identity(6);

  for (localIndex i = 0; i < 3; ++i) 
  {
    stiffness[i][i] = c1;
    stiffness[i + 3][i + 3] = 0.5 * c1;
    identity[i] = 1.0;
    identity[i + 3] = 0.0;
  }

  for (localIndex i = 0; i < 6; ++i) 
  {
    for (localIndex j = 0; j < 6; ++j) 
    {
      stiffness[i][j] += c2 * identity[i] * identity[j] + c3 * identity[i] * deviator[j]
                       + c3 * deviator[i] * identity[j];
    }
  }

  // check yield function F <= 0

  real64 yield = Q * Q / (M * M) + P * (P - oldPreconsPressure);

  if (yield > 1e-9) // plasticity branch
  {
    array1d<real64> solution(3), residual(3), delta(3);
    array2d<real64> jacobian(3, 3), jacobianInv(3, 3), hessianHyper(2, 2), hessiansMult(2, 2);
    real64 yieldDerivP, yieldDerivQ, plasticPotDerivP, plasticPotDerivQ, preconsPressureDeriv;

    solution[0] = strainVol; // initial guess for newStrainVol
    solution[1] = strainDev; // initial guess for newStrainDev
    solution[2] = 0;         // initial guess for plastic multiplier

    real64 norm, normZero = 1e30;
    jacobian = 0;

    for (localIndex iter = 0; iter < 10; ++iter) 
    {
      // apply an exponential hardening model for the preconsolidation pressure
      preconsPressure = oldPreconsPressure * std::exp((solution[0] - strainVol) / (Cc - Cr));
      preconsPressureDeriv = preconsPressure / (Cc - Cr);

      // plastic potential = associativity * Q*Q/(M*M) + P*(P-oldPreconsPressure) - check source
      yieldDerivP = 2 * P - preconsPressure;
      yieldDerivQ = 2 * Q / (M * M);
      plasticPotDerivP = yieldDerivP;
      plasticPotDerivQ = associativity * yieldDerivQ;

      // assemble residual system

      // r[0] = Eps_v - trialEps_v + lambda*dG/dP = 0
      // r[1] = Eps_s - trialEps_v + lambda*dG/dQ = 0
      // r[2] = F = 0
      residual[0] = solution[0] - strainVol + solution[2] * plasticPotDerivP;
      residual[1] = solution[1] - strainDev + solution[2] * plasticPotDerivQ;
      residual[2] = Q * Q / (M * M) + P * (P - preconsPressure);

      norm = LvArray::tensorOps::l2Norm<3>(residual);

      //residual.L2_Norm();
//      std::cout << "iter: " << iter << " " << "norm: " << norm << std::endl;

      if (iter == 0) 
      {
        normZero = norm;
      }

      if (norm < 1e-8 * (normZero + 1)) 
      {
        break;
      }

      // solve Newton system

      // Jacobian matrix depends on the components of two matrices:
      // Hessian matrix of the stored energy function for hyperelasticity = hessianHyper 
      // Hessian matrix of the plastic potential function = hessianPlastPot 
      // and the result of the multiplication: hessianPlastPot*hessianHyper = hessiansMult

      hessianHyper(0, 0) = -P / Cr;
      hessianHyper(0, 1) = 3 * refPInvariant * alpha * solution[1] * expOmega / Cr;
      hessianHyper(1, 0) = hessianHyper(0, 1);
      hessianHyper(1, 1) = 3 * shear;

      hessiansMult(0, 0) = 2 * hessianHyper(0, 0);
      hessiansMult(0, 1) = 2 * hessianHyper(0, 1);
      hessiansMult(1, 0) = 2 * hessianHyper(0, 1) * associativity / (M * M);
      hessiansMult(1, 1) = 2 * hessianHyper(1, 1) * associativity / (M * M);

      jacobian(0, 0) = 1 + solution[2] * (hessiansMult(0, 0) - preconsPressureDeriv);
      jacobian(0, 1) = solution[2] * hessiansMult(0, 1);
      jacobian(0, 2) = plasticPotDerivP;
      jacobian(1, 0) = solution[2] * hessiansMult(1, 0);
      jacobian(1, 1) = 1 + solution[2] * hessiansMult(1, 1);
      jacobian(1, 2) = plasticPotDerivQ;
      jacobian(2, 0) = hessianHyper(0, 0) * yieldDerivP + hessianHyper(1, 0) * yieldDerivQ 
                     - preconsPressureDeriv * P;
      jacobian(2, 1) = hessianHyper(0, 1) * yieldDerivP + hessianHyper(1, 1) * yieldDerivQ;
      jacobian(2, 2) = 0.0;

      LvArray::tensorOps::invert<3>(jacobianInv, jacobian);
      LvArray::tensorOps::AijBj<3, 3>(delta, jacobianInv, residual);

      for (localIndex i = 0; i < 3; ++i) 
      {
        solution[i] -= delta[i];
      }

      // update P and Q used in the residual equation:
      expOmega = std::exp((refStrainVol - solution[0]) / Cr);
      shear = refShearModulus - alpha * refPInvariant * expOmega;
      Q = 3 * shear * solution[1];
      P = refPInvariant * (1 + 3 * alpha * solution[1] * solution[1] / (2 * Cr)) * expOmega;
      bulk = -P / Cr;

      // Q: should we assert for negative Poisson?
      // real64 poissonRatio = ( 3 * bulk - 2 * shear) / ( 2 * ( 3 * bulk + shear ) );
      // GEOSX_ASSERT_MSG(poissonRatio >=0, "Negative poisson ratio produced");
    }

    // construct elastic strain
    // plasticPotDerivStress
    array1d<real64> plasticPotDerivStress(6);

    for (localIndex i = 0; i < 3; ++i)
    {
      plasticPotDerivStress[i] = std::sqrt(3. / 2.) * plasticPotDerivQ * deviator[i];
      plasticPotDerivStress[i+3] = std::sqrt(3. / 2.) * plasticPotDerivQ * 2 * deviator[i+3];
    }
    for (localIndex i = 0; i < 3; ++i) 
    {
      plasticPotDerivStress[i] += 1. / 3. * plasticPotDerivP;
    }

    for (localIndex i = 0; i < 3; ++i) 
    {
      strainElastic[i] -= solution[2] * plasticPotDerivStress[i];
      strainElastic[i+3] -= solution[2] * plasticPotDerivStress[i+3];
    }
      
    // construct stress = P*eye + sqrt(2/3)*Q*nhat

    for (localIndex i = 0; i < 6; ++i) 
    {
      stress[i] = sqrt(2. / 3.) * Q * deviator[i];
    }
    for (localIndex i = 0; i < 3; ++i) 
    {
      stress[i] += P;
    }

    // construct consistent tangent operator

    for (localIndex i = 0; i < 6; ++i) 
    {
      for (localIndex j = 0; j < 6; ++j) 
      {
        stiffness[i][j] = 0;
      }
    }

    // use of two auxiliar matrices bTilda and dBar to construct the cto
    // dBar uses the inverse of the jacobian and other constants
    real64 preconsPressureDerivTrial = - preconsPressureDeriv;
    real64 residual0DerivStrainVol  = - 1 - solution[2] * preconsPressureDerivTrial;
    real64 residual2DerivStrainVol  = - P * preconsPressureDerivTrial;
    array2d<real64> bTilda(2, 2), dBar(2, 2);
    
    bTilda(0, 0) = - ( residual0DerivStrainVol * jacobianInv(0,0) 
                     + residual2DerivStrainVol * jacobianInv(0,2) );
    bTilda(0, 1) = jacobianInv(0,1);
    bTilda(1, 0) = - ( residual0DerivStrainVol * jacobianInv(1,0) 
                     + residual2DerivStrainVol * jacobianInv(1,2) );
    bTilda(1, 1) = jacobianInv(1,1);

    //How to perform matrix multiplication with LvArray?

    for (localIndex i = 0; i < 2; ++i)
    {
      for (localIndex j = 0; j < 2; ++j)
      {
        dBar[i][j] = 0;
      }
    }

    for (localIndex i = 0; i < 2; ++i)
    {
      for (localIndex j = 0; j < 2; ++j)
      {
        for (localIndex l = 0; l < 2; ++l) 
        {
          dBar[i][j] += hessianHyper(i, l) * bTilda(l, j);
        }
      }
    }

    c1 = (std::fabs(strainDev) > 1e-10) ? 2 * Q / ( 3 * strainDev ) : 2 * shear;
    c2 = dBar(0, 0) - c1 / 3.;
    c3 = std::sqrt(2. / 3) * dBar(0, 1);
    real64 c4 = std::sqrt(2. / 3) * dBar(1, 0);
    real64 c5 = (2. / 3) * dBar(1, 1) - c1;

    for (localIndex i = 0; i < 3; ++i) 
    {
      stiffness[i][i] = c1;
      stiffness[i + 3][i + 3] = 0.5 * c1;
    }

    for (localIndex i = 0; i < 6; ++i) 
    {
      for (localIndex j = 0; j < 6; ++j) 
      {
        stiffness[i][j] += c2 * identity[i] * identity[j] + c3 * identity[i] * deviator[j]
                         + c4 * deviator[i] * identity[j] + c5 * deviator[i] * deviator[j];
      }
    }

  } // end plastic branch

  // remember history variables before returning
  m_newPreconsolidationPressure[k][q] = preconsPressure;
  for (localIndex i = 0; i < 6; ++i) 
  {
    m_newElasticStrain[k][q][i] = strainElastic[i];
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ModifiedCamClayUpdates::SaveConvergedState( localIndex const k, localIndex const q ) 
{
  m_oldPreconsolidationPressure[k][q] = m_newPreconsolidationPressure[k][q];

  for (localIndex i = 0; i < 6; ++i) 
  {
    m_oldElasticStrain[k][q][i] = m_newElasticStrain[k][q][i];
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void ModifiedCamClayUpdates::DeleteElasticStrain( localIndex const k, localIndex const q )
{
  for (localIndex i = 0; i < 6; ++i)
  {
    m_oldElasticStrain[k][q][i] = 0.0;
  }
}

/**
 * @class ModifiedCamClay
 *
 * ModifiedCamClay material model.
 */
class ModifiedCamClay : public SolidBase 
{
public:
  /// @typedef Alias for ModifiedCamClayUpdates
  using KernelWrapper = ModifiedCamClayUpdates;

  /**
   * constructor
   * @param[in] name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  ModifiedCamClay(string const &name, Group *const parent);

  /**
   * Default Destructor
   */
  virtual ~ModifiedCamClay() override;

  virtual void
  DeliverClone( string const &name, 
                Group *const parent,
                std::unique_ptr<ConstitutiveBase> &clone ) const override;

  virtual void AllocateConstitutiveData( 
                                  dataRepository::Group *const parent,
                                  localIndex const numConstitutivePointsPerParentIndex) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ModifiedCamClay";

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
    /// string/key for default reference value for the stress invariant P
    static constexpr auto defaultRefPInvariantString = "defaultRefPInvariant";

    /// string/key for default reference value for the elastic volumetric strain
    static constexpr auto defaultRefElasticStrainVolumetricString = "defaultRefElasticStrainVolumetric";

    /// string/key for default reference shear modulus
    static constexpr auto defaultRefShearModulusString = "defaultRefShearModulus";

    /// string/key for default shear modulus evolution
    static constexpr auto defaultShearModulusEvolutionString = "defaultShearModulusEvolution";

    /// string/key for default virgin compression index
    static constexpr auto defaultVirginCompressionIndexString = "defaultVirginCompressionIndex";

    /// string/key for default recompression index
    static constexpr auto defaultRecompressionIndexString = "defaultRecompressionIndex";

    /// string/key for default critical state slope
    static constexpr auto defaultCriticalStateSlopeString = "defaultCriticalStateSlope";

    /// string/key for default associativity
    static constexpr auto defaultAssociativityString = "defaultAssociativity";

    /// string/key for default preconsolidation pressure
    static constexpr auto defaultPreconsolidationPressureString = "defaultPreconsolidationPressure";

    /// string/key for reference value for the stress invariant P
    static constexpr auto refPInvariantString = "RefPInvariant";

    /// string/key for reference value for the elastic volumetric strain
    static constexpr auto refElasticStrainVolumetricString = "RefElasticStrainVolumetric";

    /// string/key for reference shear modulus
    static constexpr auto refShearModulusString = "RefShearModulus";

    /// string/key for shear modulus evolution
    static constexpr auto shearModulusEvolutionString = "ShearModulusEvolution";

    /// string/key for virgin compression index
    static constexpr auto virginCompressionIndexString = "VirginCompressionIndex";

    /// string/key for recompression index 
    static constexpr auto recompressionIndexString = "RecompressionIndex";

    /// string/key for critical state slope
    static constexpr auto criticalStateSlopeString = "CriticalStateSlope";

    /// string/key for associativity
    static constexpr auto associativityString = "Associativity";

    /// string/key for new value of preconsolidation pressure
    static constexpr auto newPreconsolidationPressureString = "NewPreconsolidationPressure";

    /// string/key for old value of preconsolidation pressure
    static constexpr auto oldPreconsolidationPressureString = "OldPreconsolidationPressure";

    /// string/key for the current elastic strain
    static constexpr auto newElasticStrainString = "NewElasticStrain";

    /// string/key for the previous elastic strain
    static constexpr auto oldElasticStrainString = "OldElasticStrain";
  };

  /**
   * @brief Create a instantiation of the ModifiedCamClayUpdate class that refers to the data in this.
   * @return An instantiation of ModifiedCamClayUpdate.
   */
  ModifiedCamClayUpdates createKernelWrapper() 
  {
    return ModifiedCamClayUpdates( m_referencePInvariant, 
                                   m_refElasticStrainVolumetric, 
                                   m_referenceShearModulus, 
                                   m_shearModulusEvolution, 
                                   m_virginCompressionIndex, 
                                   m_recompressionIndex, 
                                   m_criticalStateSlope,
                                   m_associativity, 
                                   m_newPreconsolidationPressure, 
                                   m_oldPreconsolidationPressure,
                                   m_newElasticStrain, 
                                   m_oldElasticStrain,
                                   m_stress);
                                  // m_strainElastic );
  }

protected:
  virtual void PostProcessInput() override;

private:
  // TODO: maybe define some helper structs for defaults, arrays, arrayViews, etc.

  /// Material parameter: The default value of the bulk modulus
  real64 m_defaultRefPInvariant;

  /// Material parameter: The default value of the shear modulus
  real64 m_defaultRefElasticStrainVolumetric;

  /// Material parameter: The default value of the yield surface slope
  real64 m_defaultRefShearModulus;

  /// Material parameter: The default value of the parameter setting the shear modulus evolution
  real64 m_defaultShearModulusEvolution;

  /// Material parameter: The default value of the virgin compression index 
  real64 m_defaultVirginCompressionIndex;

  /// Material parameter: The default value of the recompression index 
  real64 m_defaultRecompressionIndex;

  /// Material parameter: The default value of the critical state line slope 
  real64 m_defaultCriticalStateSlope;

  /// Material parameter: The default value of the critical state line slope 
  real64 m_defaultAssociativity;

  /// History variable: The default value of the preconsolidation pressure
  real64 m_defaultPreconsolidationPressure;

  /// Material parameter: The reference value for the stress invariant P for quadrature point
  array2d<real64> m_referencePInvariant;

  /// Material parameter: The reference value for the volumetric part of the elastic strain for each quadrature point
  array2d<real64> m_refElasticStrainVolumetric;

  /// Material parameter: The reference value for the shear modulus for each quadrature point
  array2d<real64> m_referenceShearModulus;

  /// Material parameter: The shear modulus evolution for each element
  array1d<real64> m_shearModulusEvolution;

  /// Material parameter: The virgin compression index (Cc) for each element
  array1d<real64> m_virginCompressionIndex;

  /// Material parameter: The recompression index (Cr) for each element
  array1d<real64> m_recompressionIndex;

  /// Material parameter: The critical state line slope for each element
  array1d<real64> m_criticalStateSlope;
  
  /// Material parameter: Parameter to accommodate non-associated plasticity for each element
  array1d<real64> m_associativity;

  /// History variable: The current preconsolidation pressure for each quadrature point
  array2d<real64> m_newPreconsolidationPressure;

  /// History variable: The previous preconsolidation pressure for each quadrature point
  array2d<real64> m_oldPreconsolidationPressure;

  /// History variable: The current elastic strain for each quadrature point
  array3d<real64, solid::STRESS_PERMUTATION> m_newElasticStrain;

  /// History variable: The previous elastic strain for each quadrature point.
  array3d<real64, solid::STRESS_PERMUTATION> m_oldElasticStrain;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_MODIFIEDCAMCLAY_HPP_ */
