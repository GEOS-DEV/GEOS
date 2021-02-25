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
    m_refPressure(refPressure),
    m_refStrainVol(refStrainVol),
    m_recompressionIndex(recompressionIndex),
    m_virginCompressionIndex(virginCompressionIndex),
    m_cslSlope(cslSlope),
    m_shapeParameter(shapeParameter),
    m_newPreConsolidationPressure(newPreConsolidationPressure),
    m_oldPreConsolidationPressure(oldPreConsolidationPressure)
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
    
  /// A reference to the ArrayView holding the reference pressure for each element.
  arrayView1d< real64 const > const m_refPressure;
    
  /// A reference to the ArrayView holding the reference volumetric strain for each element.
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
  // elastic predictor (assume strainIncrement is all elastic)

  ElasticIsotropicUpdates::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

  // decompose into mean (P) and von Mises (Q) stress invariants

  real64 trialP;
  real64 trialQ;
  real64 deviator[6];

  twoInvariant::stressDecomposition( stress,
                                     trialP,
                                     trialQ,
                                     deviator );

  // check yield function F <= 0, using old hardening variable state

  real64 yield = trialQ + m_friction[k] * trialP - m_oldCohesion[k][q];

  if( yield < 1e-9 ) // elasticity
  {
    return;
  }

  // else, plasticity (trial stress point lies outside yield surface)

  // the return mapping can in general be written as a Newton iteration.
  // here we have a linear problem, so the algorithm will converge in one
  // iteration, but this is a template for more general models with either
  // nonlinear hardening or yield surfaces.

  real64 solution[3], residual[3], delta[3];
  real64 jacobian[3][3] = {{}}, jacobianInv[3][3] = {{}};

  solution[0] = trialP; // initial guess for newP
  solution[1] = trialQ; // initial guess for newQ
  solution[2] = 1e-5;   // initial guess for plastic multiplier

  real64 norm, normZero = 1e30;

  // begin Newton loop

  for( localIndex iter=0; iter<20; ++iter )
  {
    // apply a linear cohesion decay model,
    // then check for complete cohesion loss

    m_newCohesion[k][q] = m_oldCohesion[k][q] + solution[2] * m_hardening[k];
    real64 cohesionDeriv = m_hardening[k];

    if( m_newCohesion[k][q] < 0 )
    {
      m_newCohesion[k][q] = 0;
      cohesionDeriv = 0;
    }

    // assemble residual system
    // resid1 = P - trialP + dlambda*bulkMod*dG/dP = 0
    // resid2 = Q - trialQ + dlambda*3*shearMod*dG/dQ = 0
    // resid3 = F = 0

    residual[0] = solution[0] - trialP + solution[2] * m_bulkModulus[k] * m_dilation[k];
    residual[1] = solution[1] - trialQ + solution[2] * 3 * m_shearModulus[k];
    residual[2] = solution[1] + m_friction[k] * solution[0] - m_newCohesion[k][q];

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

    jacobian[0][0] = 1;
    jacobian[0][2] = m_bulkModulus[k] * m_dilation[k];
    jacobian[1][1] = 1;
    jacobian[1][2] = 3 * m_shearModulus[k];
    jacobian[2][0] = m_friction[k];
    jacobian[2][1] = 1;
    jacobian[2][2] = -cohesionDeriv;

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
    return CamClayUpdates(m_refPressure
                          m_refStrainVol
                          m_recompressionIndex
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
                          m_refPressure
                          m_refStrainVol
                          m_recompressionIndex
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
    
  /// Material parameter: The reference pressure for each element
  array1d< real64 > m_refPressure;
    
  /// Material parameter: The reference volumetric strain for each element
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
