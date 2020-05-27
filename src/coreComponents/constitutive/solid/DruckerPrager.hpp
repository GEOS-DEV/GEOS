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

#include "SolidBase.hpp"

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
class DruckerPragerUpdates : public SolidBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] stress The ArrayView holding the stress data for each quadrature point.
   */
  DruckerPragerUpdates( arrayView1d< real64 const > const & bulkModulus,
                        arrayView1d< real64 const > const & poissonRatio,
                        arrayView1d< real64 const > const & tanFrictionAngle,
                        arrayView1d< real64 const > const & hardeningRate,
                        arrayView2d< real64 > const & newCohesion,
                        arrayView2d< real64 > const & oldCohesion,
                        arrayView3d< real64, solid::STRESS_USD > const & newStress,
                        arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                        arrayView3d< real64, solid::STRESS_USD > const & stress ):
    SolidBaseUpdates( stress ),
    m_bulkModulus( bulkModulus ),
    m_poissonRatio( poissonRatio ),
    m_tanFrictionAngle( tanFrictionAngle ),
    m_hardeningRate( hardeningRate ),
    m_newCohesion( newCohesion ),
    m_oldCohesion( oldCohesion ),
    m_newStress( newStress ),
    m_oldStress( oldStress )
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

  /**
   * accessor to return the stiffness at a given element
   * @param[in] k the element number
   * @param[in] c the stiffness array
   */
  GEOSX_HOST_DEVICE inline
  virtual void GetStiffness( localIndex const k, real64 (& c)[6][6] ) const override final
  {
    real64 const G = m_poissonRatio[k];
    real64 const Lame = m_bulkModulus[k] - 2.0/3.0 * G;

    memset( c, 0, sizeof( c ) );

    c[0][0] = Lame + 2 * G;
    c[0][1] = Lame;
    c[0][2] = Lame;

    c[1][0] = Lame;
    c[1][1] = Lame + 2 * G;
    c[1][2] = Lame;

    c[2][0] = Lame;
    c[2][1] = Lame;
    c[2][2] = Lame + 2 * G;

    c[3][3] = G;

    c[4][4] = G;

    c[5][5] = G;
  }
  
  GEOSX_HOST_DEVICE
  virtual void SmallStrainNoState( localIndex const k,
                                   real64 const * const GEOSX_RESTRICT voigtStrain,
                                   real64 * const GEOSX_RESTRICT stress ) const override final;

  GEOSX_HOST_DEVICE
  virtual void SmallStrain( localIndex const k,
                            localIndex const q,
                            real64 const * const GEOSX_RESTRICT voigtStrainIncrement ) const override final;

  GEOSX_HOST_DEVICE
  virtual void HypoElastic( localIndex const k,
                            localIndex const q,
                            real64 const * const GEOSX_RESTRICT Ddt,
                            R2Tensor const & Rot ) const override final;

  GEOSX_HOST_DEVICE
  virtual void HyperElastic( localIndex const k,
                             real64 const (&FmI)[3][3],
                             real64 * const GEOSX_RESTRICT stress ) const override final;

  GEOSX_HOST_DEVICE
  virtual void HyperElastic( localIndex const k,
                             localIndex const q,
                             real64 const (&FmI)[3][3] ) const override final;
  
  
private:
  /// A reference to the ArrayView holding the bulk modulus for each element.
  arrayView1d< real64 const > const m_bulkModulus;

  /// A reference to the ArrayView holding the shear modulus for each element.
  arrayView1d< real64 const > const m_poissonRatio;
  
  /// A reference to the ArrayView holding the shear modulus for each element.
  arrayView1d< real64 const > const m_tanFrictionAngle;
  
  /// A reference to the ArrayView holding the hardening rate for each element.
  arrayView1d< real64 const > const m_hardeningRate;
  
  /// A reference to the ArrayView holding the new cohesion for each integration point
  arrayView2d< real64 > const m_newCohesion;
  
  /// A reference to the ArrayView holding the old cohesion for each integration point
  arrayView2d< real64 > const m_oldCohesion;
  
  /// A reference to the ArrayView holding the new stress for each integration point
  arrayView3d< real64, solid::STRESS_USD > const m_newStress;
  
  /// A reference to the ArrayView holding the old stress for each integration point
  arrayView3d< real64, solid::STRESS_USD > const m_oldStress;
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::SmallStrainNoState( localIndex const k,
                                                        real64 const * GEOSX_RESTRICT const voigtStrain,
                                                        real64 * GEOSX_RESTRICT const stress ) const
{
  real64 const lambda = m_bulkModulus[k] - 2.0/3.0 * m_poissonRatio[k];
  real64 const diag = lambda * ( voigtStrain[0] + voigtStrain[1] + voigtStrain[2] );
  real64 const TwoG = 2.0 * m_poissonRatio[k];

  stress[0] = stress[0] + diag + TwoG * voigtStrain[0];
  stress[1] = stress[1] + diag + TwoG * voigtStrain[1];
  stress[2] = stress[2] + diag + TwoG * voigtStrain[2];
  stress[3] = stress[3] + m_poissonRatio[k] * voigtStrain[3];
  stress[4] = stress[4] + m_poissonRatio[k] * voigtStrain[4];
  stress[5] = stress[5] + m_poissonRatio[k] * voigtStrain[5];

}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::SmallStrain( localIndex const k,
                                                 localIndex const q,
                                                 real64 const * const GEOSX_RESTRICT voigtStrainInc ) const
{
  real64 const lambda = m_bulkModulus[k] - 2.0/3.0 * m_poissonRatio[k];
  real64 const volStrain = ( voigtStrainInc[0] + voigtStrainInc[1] + voigtStrainInc[2] );
  real64 const TwoG = 2.0 * m_poissonRatio[k];

  m_stress( k, q, 0 ) =  m_stress( k, q, 0 ) + TwoG * voigtStrainInc[0] + lambda * volStrain;
  m_stress( k, q, 1 ) =  m_stress( k, q, 1 ) + TwoG * voigtStrainInc[1] + lambda * volStrain;
  m_stress( k, q, 2 ) =  m_stress( k, q, 2 ) + TwoG * voigtStrainInc[2] + lambda * volStrain;
  m_stress( k, q, 3 ) =  m_stress( k, q, 3 ) + m_poissonRatio[k] * voigtStrainInc[3];
  m_stress( k, q, 4 ) =  m_stress( k, q, 4 ) + m_poissonRatio[k] * voigtStrainInc[4];
  m_stress( k, q, 5 ) =  m_stress( k, q, 5 ) + m_poissonRatio[k] * voigtStrainInc[5];
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::HypoElastic( localIndex const k,
                                                 localIndex const q,
                                                 real64 const * const GEOSX_RESTRICT Ddt,
                                                 R2Tensor const & Rot ) const
{
  real64 const lambda = m_bulkModulus[k] - 2.0/3.0 * m_poissonRatio[k];
  real64 const volStrain = ( Ddt[0] + Ddt[2] + Ddt[5] );
  real64 const TwoG = 2.0 * m_poissonRatio[k];


  m_stress( k, q, 0 ) =  m_stress( k, q, 0 ) + TwoG * Ddt[0] + lambda * volStrain;
  m_stress( k, q, 1 ) =  m_stress( k, q, 1 ) + TwoG * Ddt[2] + lambda * volStrain;
  m_stress( k, q, 2 ) =  m_stress( k, q, 2 ) + TwoG * Ddt[5] + lambda * volStrain;
  m_stress( k, q, 3 ) =  m_stress( k, q, 3 ) + TwoG * Ddt[4];
  m_stress( k, q, 4 ) =  m_stress( k, q, 4 ) + TwoG * Ddt[3];
  m_stress( k, q, 5 ) =  m_stress( k, q, 5 ) + TwoG * Ddt[1];

  R2SymTensor stress;
  stress = m_stress[k][q];

  R2SymTensor temp;
  real64 const * const pTemp = temp.Data();
  temp.QijAjkQlk( stress, Rot );

  m_stress( k, q, 0 ) = pTemp[0];
  m_stress( k, q, 1 ) = pTemp[2];
  m_stress( k, q, 2 ) = pTemp[5];
  m_stress( k, q, 3 ) = pTemp[4];
  m_stress( k, q, 4 ) = pTemp[3];
  m_stress( k, q, 5 ) = pTemp[1];
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::HyperElastic( localIndex const k,
                                                  real64 const (&FmI)[3][3],
                                                  real64 * const GEOSX_RESTRICT stress ) const
{
  real64 const C1 = 0.5 * m_poissonRatio[k];
  real64 const D1 = 0.5 * m_bulkModulus[k];
  real64 const detFm1 = FmI[0][0] + FmI[1][1] + FmI[2][2]
                        - FmI[1][2]*FmI[2][1] + FmI[1][1]*FmI[2][2]
                        + FmI[0][2]*(-FmI[2][0] - FmI[1][1]*FmI[2][0] + FmI[1][0]*FmI[2][1])
                        + FmI[0][1]*(-FmI[1][0] + FmI[1][2]*FmI[2][0] - FmI[1][0]*FmI[2][2])
                        + FmI[0][0]*( FmI[1][1] - FmI[1][2]*FmI[2][1] + FmI[2][2] + FmI[1][1]*FmI[2][2]);


  real64 const p = -2 * D1 * ( detFm1 + 1.0 ) * detFm1;
  real64 devB[6] = { 1/3 * (2 * FmI[0][0] * (2 + FmI[0][0]) - FmI[1][1] * (2 + FmI[1][1]) - FmI[2][2] * (2 + FmI[2][2]) +
                            2 * FmI[0][1]*FmI[0][1] + 2 * FmI[0][2] * FmI[0][2] - FmI[1][0] * FmI[1][0] - FmI[1][2] * FmI[1][2] -
                            FmI[2][0] * FmI[2][0] - FmI[2][1] * FmI[2][1]),
                     1/3 * (-FmI[0][0] * (2 + FmI[0][0]) + 2 * FmI[1][1] * ( 2 + FmI[1][1]) - FmI[2][2] * (2 + FmI[2][2]) -
                            FmI[0][1]*FmI[0][1] - FmI[0][2]*FmI[0][2] + 2 * FmI[1][0]*FmI[1][0] + 2 * FmI[1][2]*FmI[1][2] - FmI[2][0]*FmI[2][0] - FmI[2][1]*
                            FmI[2][1]),
                     1/3 *(-FmI[0][0] * (2 + FmI[0][0]) - FmI[1][1] * (2 + FmI[1][1]) + 2 * FmI[2][2] * (2 + FmI[2][2]) -
                           FmI[0][1]*FmI[0][1] - FmI[0][2]*FmI[0][2] - FmI[1][0]*FmI[1][0] - FmI[1][2]*FmI[1][2] + 2 * FmI[2][0]*FmI[2][0] + 2 * FmI[2][1]*
                           FmI[2][1]),
                     FmI[1][2] + FmI[1][0]*FmI[2][0] + FmI[2][1] + FmI[1][1]*FmI[2][1] + FmI[1][2]*FmI[2][2],
                     FmI[0][2] + FmI[2][0] + FmI[0][0]*FmI[2][0] + FmI[0][1]*FmI[2][1] + FmI[0][2]*FmI[2][2],
                     FmI[0][1] + FmI[1][0] + FmI[0][0]*FmI[1][0] + FmI[0][1]*FmI[1][1] + FmI[0][2]*FmI[1][2]
  };

  real64 const C = 2 * C1 / pow( detFm1 + 1, 2.0/3.0 );
  stress[0] = -p + C * devB[0];
  stress[1] = -p + C * devB[1];
  stress[2] = -p + C * devB[2];
  stress[3] = C * devB[3];
  stress[4] = C * devB[4];
  stress[5] = C * devB[5];
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::HyperElastic( localIndex const k,
                                                  localIndex const q,
                                                  real64 const (&FmI)[3][3] ) const
{
  real64 stress[6];
  HyperElastic( k, FmI, stress );

  for( localIndex i=0; i<6; ++i )
  {
    m_stress( k, q, i ) = stress[i];
  }
}



/**
 * @class DruckerPrager
 *
 * Drucker-Prager material model.
 */
class DruckerPrager : public SolidBase
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
    /// string/key for default bulk modulus
    static constexpr auto defaultBulkModulusString  = "defaultBulkModulus";

    /// string/key for default shear modulus
    static constexpr auto defaultPoissonRatioString = "defaultPoissonRatio";

    /// string/key for default friction angle
    static constexpr auto defaultTanFrictionAngleString = "defaultTanFrictionAngle";
    
    /// string/key for default hardening rate
    static constexpr auto defaultHardeningRateString = "defaultHardeningRate";
    
    /// string/key for default cohesion
    static constexpr auto defaultCohesionString = "defaultCohesion";
    
    /// string/key for bulk modulus
    static constexpr auto bulkModulusString  = "BulkModulus";
    
    /// string/key for poisson ratio
    static constexpr auto poissonRatioString = "PoissonRatio";
    
    /// string/key for friction angle
    static constexpr auto tanFrictionAngleString  = "TanFrictionAngle";
    
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
   * @brief Create a instantiation of the DruckerPragerUpdate class
   *        that refers to the data in this.
   * @return An instantiation of DruckerPragerUpdate.
   */
  DruckerPragerUpdates createKernelWrapper()
  {
    return DruckerPragerUpdates( m_bulkModulus,
                                 m_poissonRatio,
                                 m_tanFrictionAngle,
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
  
  /// Material parameter: The default value of the bulk modulus
  real64 m_defaultBulkModulus;

  /// Material parameter: The default value of the Poisson ratio
  real64 m_defaultPoissonRatio;
  
  /// Material parameter: The default value of yield surface slope
  real64 m_defaultTanFrictionAngle;
  
  /// Material parameter: The default value of the initial cohesion
  real64 m_defaultCohesion;
  
  /// Material parameter: The default value of the hardening rate
  real64 m_defaultHardeningRate;
  
  /// Material parameter: The bulk modulus for each element
  array1d< real64 > m_bulkModulus;

  /// Material parameter: The shear modulus for each element
  array1d< real64 > m_poissonRatio;
  
  /// Material parameter: The yield surface slope for each element
  array1d< real64 > m_tanFrictionAngle;
  
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
  
  /* TODO: maybe define helper structs for defaults, arrays, arrayViews, etc.
  struct DefaultParameters
  {
    real64 bulkModulus;
    real64 poissonRatio;
    real64 tanFrictionAngle;
    real64 cohesion;
    real64 hardeningRate;
  }
  m_default;
  */
  
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGER_HPP_ */
