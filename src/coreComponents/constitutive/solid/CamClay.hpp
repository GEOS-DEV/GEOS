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
 *  @file CamClay.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_CAMCLAY_HPP
#define GEOSX_CONSTITUTIVE_SOLID_CAMCLAY_HPP

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
 * @class CamClayUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class CamClayUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] stress The ArrayView holding the stress data for each quadrature point.
   */
  CamClayUpdates( arrayView1d< real64 const > const & refPressure,
                  arrayView1d< real64 const > const & refStrainVol,
                  arrayView1d< real64 const > const & recompressionIndex,
                  arrayView1d< real64 const > const & virginCompressionIndex,
                  arrayView1d< real64 const > const & cslSlope,
                  arrayView1d< real64 const > const & shapeParameter,
                  arrayView2d< real64 > const & newPreConsolidationPressure,
                  arrayView2d< real64 > const & oldPreConsolidationPressure,
                  arrayView1d< real64 const > const & bulkModulus,
                  arrayView1d< real64 const > const & shearModulus,
                  arrayView3d< real64, solid::STRESS_USD > const & newStress,
                  arrayView3d< real64, solid::STRESS_USD > const & oldStress ):
    ElasticIsotropicUpdates( bulkModulus, shearModulus, newStress, oldStress ),
    m_refPressure( refPressure ),
    m_refStrainVol( refStrainVol ),
    m_recompressionIndex( recompressionIndex ),
    m_virginCompressionIndex( virginCompressionIndex ),
    m_cslSlope( cslSlope ),
    m_shapeParameter( shapeParameter ),
    m_newPreConsolidationPressure( newPreConsolidationPressure ),
    m_oldPreConsolidationPressure( oldPreConsolidationPressure )
  {}

  /// Default copy constructor
  CamClayUpdates( CamClayUpdates const & ) = default;

  /// Default move constructor
  CamClayUpdates( CamClayUpdates && ) = default;

  /// Deleted default constructor
  CamClayUpdates() = delete;

  /// Deleted copy assignment operator
  CamClayUpdates & operator=( CamClayUpdates const & ) = delete;

  /// Deleted move assignment operator
  CamClayUpdates & operator=( CamClayUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using ElasticIsotropicUpdates::smallStrainUpdate;
  using ElasticIsotropicUpdates::saveStress;

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

  /// A reference to the ArrayView holding the reference pressure for each quadrature point.
  arrayView1d< real64 const > const m_refPressure;

  /// A reference to the ArrayView holding the reference volumetric strain for each quadrature point.
  arrayView1d< real64 const > const m_refStrainVol;

  /// A reference to the ArrayView holding the recompression index  for each element.
  arrayView1d< real64 const > const m_recompressionIndex;

  /// A reference to the ArrayView holding the virgin compression index for each element.
  arrayView1d< real64 const > const m_virginCompressionIndex;

  /// A reference to the ArrayView holding the slope of the critical state line for each element.
  arrayView1d< real64 const > const m_cslSlope;

  /// A reference to the ArrayView holding the shape parameter for each element.
  arrayView1d< real64 const > const m_shapeParameter;

  /// A reference to the ArrayView holding the new preconsolidation pressure for each integration point
  arrayView2d< real64 > const m_newPreConsolidationPressure;

  /// A reference to the ArrayView holding the old preconsolidation presure for each integration point
  arrayView2d< real64 > const m_oldPreConsolidationPressure;

};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CamClayUpdates::smallStrainUpdate( localIndex const k,
                                        localIndex const q,
                                        real64 const ( &strainIncrement )[6],
                                        real64 ( & stress )[6],
                                        real64 ( & stiffness )[6][6] ) const
{

  // Rename variables for easier implementation

  real64 const oldPc  = m_oldPreConsolidationPressure[k][q];   //pre-consolidation pressure
  real64 const mu     = m_shearModulus[k];
  real64 const p0     = m_refPressure[k];

  real64 const eps_v0 = m_refStrainVol[k];
  real64 const M      = m_cslSlope[k];
  real64 const Cr     = m_recompressionIndex[k];
  real64 const Cc     = m_virginCompressionIndex[k];
  real64 const alpha  = m_shapeParameter[k];

  real64 pc    = oldPc;
  real64 bulkModulus  = -p0/Cr;


  // two-invariant decomposition of old stress in P-Q space (mean & deviatoric stress)

  real64 oldP;
  real64 oldQ;
  real64 oldDeviator[6];
  real64 deviator[6];
  real64 oldStrainElastic[6];
  real64 strainElasticTrial[6];
  real64 eps_s_trial;
  real64 eps_v_trial;

  twoInvariant::stressDecomposition( stress,
                                     oldP,
                                     oldQ,
                                     oldDeviator );



  // Recover elastic strains from the previous step, based on stress from the previous step
  // [Note: in order to minimize data transfer, we are not storing and passing elastic strains]

  real64 oldElasticStrainVol = std::log( oldP/p0 ) * Cr * (-1.0) + eps_v0;
  real64 oldElasticStrainDev = oldQ/3./mu;

  // Now recover the old strain tensor from the strain invariants.
  // Note that we need the deviatoric direction (n-hat) from the previous step.

  twoInvariant::strainRecomposition( oldElasticStrainVol,
                                     oldElasticStrainDev,
                                     oldDeviator,
                                     oldStrainElastic );

  // elastic predictor (assume strainIncrement is all elastic)
  // TODO: define ElasticIsotropicPressureUpdates to call here
  //ElasticIsotropicUpdates::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

  for( localIndex i=0; i<6; ++i )
  {
    strainElasticTrial[i] = oldStrainElastic[i] + strainIncrement[i];
  }
  // two-invariant decomposition of trial elastic strain

  twoInvariant::strainDecomposition( strainElasticTrial,
                                     eps_v_trial,
                                     eps_s_trial,
                                     deviator );

  // Calculate trial mean and deviatoric stress

  real64 trialP = p0 * std::exp( -1./Cr* (eps_v_trial-eps_v0));
  real64 trialQ = 3. * mu * eps_s_trial;

  twoInvariant::stressRecomposition( trialP,
                                     trialQ,
                                     deviator,
                                     stress );

  // set stiffness to elastic predictor

  bulkModulus = -trialP/Cr;
  real64 lame = bulkModulus - 2./3. * mu;

  for( localIndex i=0; i<6; ++i )
  {
    for( localIndex j=0; j<6; ++j )
    {
      stiffness[i][j] = 0;
    }
  }

  stiffness[0][0] = lame + 2.*mu;
  stiffness[0][1] = lame;
  stiffness[0][2] = lame;

  stiffness[1][0] = lame;
  stiffness[1][1] = lame + 2.*mu;
  stiffness[1][2] = lame;

  stiffness[2][0] = lame;
  stiffness[2][1] = lame;
  stiffness[2][2] = lame + 2.*mu;

  stiffness[3][3] = mu;
  stiffness[4][4] = mu;
  stiffness[5][5] = mu;

  // check yield function F <= 0
  real64 yield = trialQ*trialQ/(M*M)- alpha*alpha*trialP *(2*alpha/(alpha+1)*pc-trialP)+alpha*alpha*(alpha-1)/(alpha+1)* pc*pc;

  if( yield < 1e-9 ) // elasticity
  {
      std::cout << "elastic" <<  "\n " << std::endl;

    return;
  }

  // else, plasticity (trial stress point lies outside yield surface)
  std::cout << "plastic " <<  "\n " << std::endl;


  real64 solution[3], residual[3], delta[3];
  real64 jacobian[3][3] = {{}}, jacobianInv[3][3] = {{}};

  solution[0] = eps_v_trial; // initial guess for elastic volumetric strain
  solution[1] = eps_s_trial; // initial guess for elastic deviatoric strain
  solution[2] = 1e-5;   // initial guess for plastic multiplier

  real64 norm, normZero = 1e30;

  // begin Newton loop

  for( localIndex iter=0; iter<20; ++iter )
  {
    trialP = p0 * std::exp( -1./Cr* (solution[0] - eps_v0));
    trialQ = 3. * mu * solution[1];
    bulkModulus = -trialP/Cr;
    pc = oldPc * std::exp( -1./(Cc-Cr)*(eps_v_trial-solution[0]));

    yield = trialQ*trialQ/(M*M)- alpha*alpha*trialP *(2.*alpha/(alpha+1.)*pc-trialP)+alpha*alpha*(alpha-1.)/(alpha+1.)* pc*pc;

    // derivatives of yield surface
    real64 alphaTerm = 2. * alpha*alpha*alpha / (alpha+1.);
    real64 df_dp = -alphaTerm * pc + 2. * alpha * alpha* trialP;
    real64 df_dq = 2. * trialQ /(M*M);
    real64 df_dpc = 2. * alpha*alpha*(alpha-1.) /(alpha+1.) * pc - alphaTerm * trialP;
    real64 dpc_dve = -1./(Cc-Cr) * pc;   //TODO: Check negative or positive

    real64 df_dp_dve = 2. * alpha * alpha * bulkModulus - alphaTerm * dpc_dve;
    real64 df_dq_dse = 2. /(M*M) * 3. * mu;
    //real64 df_dpc_dve = -alphaTerm * bulkModulus + 2*alpha*alpha*(alpha-1) /(alpha+1) * dpc_dve;


    // assemble residual system
    residual[0] = solution[0] - eps_v_trial + solution[2]*df_dp;   // strainElasticDev - strainElasticTrialDev + dlambda*dG/dPQ = 0
    residual[1] = solution[1] - eps_s_trial + solution[2]*df_dq;         // strainElasticVol - strainElasticTrialVol + dlambda*dG/dQ = 0
    residual[2] = yield;      // F = 0


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

    jacobian[0][0] = 1. + solution[2] * df_dp_dve;
    jacobian[0][2] = df_dp;
    jacobian[1][1] = 1. + solution[2]*df_dq_dse;
    jacobian[1][2] = df_dq;
    jacobian[2][0] = bulkModulus * df_dp - dpc_dve * df_dpc;
    jacobian[2][1] = 3.0 * mu * df_dq;
    jacobian[2][2] = 0.0;

    LvArray::tensorOps::invert< 3 >( jacobianInv, jacobian );
    LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( delta, jacobianInv, residual );

    for( localIndex i=0; i<3; ++i )
    {
      solution[i] -= delta[i];
    }
  }

  // re-construct stress = P*eye + sqrt(2/3)*Q*nhat

  twoInvariant::stressRecomposition( trialP,
                                     trialQ,
                                     deviator,
                                     stress );

  // construct consistent tangent stiffness

  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );
  real64 BB[2][2] = {{}};

  real64 dpc_dve = -1./(Cc-Cr) * pc;
  real64 a1= 1. + solution[2]*dpc_dve;
  real64 a2 = trialP * dpc_dve;

  bulkModulus = -trialP/Cr;

  BB[0][0] = bulkModulus*(a1*jacobianInv[0][0]+a2*jacobianInv[0][2]);
  BB[0][1] =bulkModulus*jacobianInv[0][1];
  BB[1][0] =3. * mu*(a1*jacobianInv[1][0]+a2*jacobianInv[1][2]);
  BB[1][1] = 3. * mu*jacobianInv[1][1];

  real64 c1;

  if( eps_s_trial<1e-15 ) // confirm eps_s_trial != 0
  {
    c1 = 2. * mu;
  }
  else
  {
    c1 = 2. * trialQ/(3. * eps_s_trial);
  }

  real64 c2 = BB[0][0] - c1/3.;
  real64 c3 = std::sqrt( 2./3. ) * BB[0][1];
  real64 c4 = std::sqrt( 2./3. ) * BB[1][0];
  real64 c5 = 2./3. * BB[1][1] - c1;

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
  // remember history variables before returning

  m_newPreConsolidationPressure[k][q] = pc;

  // save new stress and return
  saveStress( k, q, stress );
  return;
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CamClayUpdates::smallStrainUpdate( localIndex const k,
                                        localIndex const q,
                                        real64 const ( &strainIncrement )[6],
                                        real64 ( & stress )[6],
                                        DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, strainIncrement, stress, stiffness.m_c );
}



/**
 * @class CamClay
 *
 * Modified Cam-Clay and Delft-Egg material model.
 */
class CamClay : public ElasticIsotropic
{
public:

  /// @typedef Alias for CamClayUpdates
  using KernelWrapper = CamClayUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  CamClay( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~CamClay() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "CamClay";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default friction angle
    static constexpr char const * defaultRefPressureString() { return "defaultRefPressure"; }

    /// string/key for default dilation angle
    static constexpr char const * defaultRefStrainVolString() { return "defaultRefStrainVol"; }

    /// string/key for default hardening rate
    static constexpr char const * defaultRecompressionIndexString() { return "defaultRecompressionIndex"; }

    /// string/key for default cohesion
    static constexpr char const * defaultVirginCompressionIndexString() { return "defaultVirginCompressionIndex"; }

    /// string/key for default cohesion
    static constexpr char const * defaultCslSlopeString() { return "defaultCslSlope"; }

    /// string/key for default cohesion
    static constexpr char const * defaultShapeParameterString() { return "defaultShapeParameter"; }

    /// string/key for default cohesion
    static constexpr char const * defaultPreConsolidationPressureString() { return "defaultPreConsolidationPressure"; }

    /// string/key for friction angle
    static constexpr char const * refPressureString() { return "refPressure"; }

    /// string/key for dilation angle
    static constexpr char const * refStrainVolString() { return "refStrainVol"; }

    /// string/key for cohesion
    static constexpr char const * recompressionIndexString() { return "recompressionIndex"; }

    /// string/key for cohesion
    static constexpr char const * virginCompressionIndexString() { return "virginCompressionIndex"; }

    /// string/key for cohesion
    static constexpr char const * cslSlopeString() { return "cslSlope"; }

    /// string/key for cohesion
    static constexpr char const * shapeParameterString() { return "shapeParameter"; }

    /// string/key for cohesion
    static constexpr char const * newPreConsolidationPressureString() { return "preConsolidationPressure"; }

    /// string/key for cohesion
    static constexpr char const * oldPreConsolidationPressureString() { return "oldPreConsolidationPressure"; }
  };

  /**
   * @brief Create a instantiation of the CamClayUpdate class that refers to the data in this.
   * @return An instantiation of CamClayUpdate.
   */
  CamClayUpdates createKernelUpdates() const
  {
    return CamClayUpdates( m_refPressure,
                           m_refStrainVol,
                           m_recompressionIndex,
                           m_virginCompressionIndex,
                           m_cslSlope,
                           m_shapeParameter,
                           m_newPreConsolidationPressure,
                           m_oldPreConsolidationPressure,
                           m_bulkModulus,
                           m_shearModulus,
                           m_newStress,
                           m_oldStress );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_refPressure,
                          m_refStrainVol,
                          m_recompressionIndex,
                          m_virginCompressionIndex,
                          m_cslSlope,
                          m_shapeParameter,
                          m_newPreConsolidationPressure,
                          m_oldPreConsolidationPressure,
                          m_bulkModulus,
                          m_shearModulus,
                          m_newStress,
                          m_oldStress );
  }


protected:
  virtual void postProcessInput() override;

  /// Material parameter: The default value of reference pressure
  real64 m_defaultRefPressure;

  /// Material parameter: The default value of reference volumetric strain
  real64 m_defaultRefStrainVol;

  /// Material parameter: The default value of the recompression index
  real64 m_defaultRecompressionIndex;

  /// Material parameter: The default value of the virgin compression index
  real64 m_defaultVirginCompressionIndex;

  /// Material parameter: The default value of the slope of the critical state line
  real64 m_defaultCslSlope;

  /// Material parameter: The default value of the shape parameter of the yield surface
  real64 m_defaultShapeParameter;

  /// Material parameter: The default value of the preconsolidation pressure
  real64 m_defaultPreConsolidationPressure;

  /// Material parameter: The reference pressure for each quadrature point
  array1d< real64 > m_refPressure;

  /// Material parameter: The reference volumetric strain for each quadrature point
  array1d< real64 > m_refStrainVol;

  /// Material parameter: The recompression index for each element
  array1d< real64 > m_recompressionIndex;

  /// Material parameter: The virgin compression index for each element
  array1d< real64 > m_virginCompressionIndex;

  /// Material parameter: The slope of the critical state line for each element
  array1d< real64 > m_cslSlope;

  /// Material parameter: Thehape parameter of the yield surface for each element
  array1d< real64 > m_shapeParameter;

  /// State variable: The current preconsolidation pressure for each quadrature point
  array2d< real64 > m_newPreConsolidationPressure;

  /// State variable: The previous preconsolidation pressure for each quadrature point
  array2d< real64 > m_oldPreConsolidationPressure;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_CAMCLAY_HPP_ */
