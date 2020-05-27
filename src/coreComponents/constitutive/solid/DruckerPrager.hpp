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

  GEOSX_HOST_DEVICE
  virtual void SmallStrain( localIndex const k,
                            localIndex const q,
                            arraySlice1d< real64 const > const & strainIncrement,
                            arraySlice1d< real64 > const & stress,
                            arraySlice2d< real64 > const & stiffness ) override final;
  
  // remaining interface functions not implemented for various reasons
  // but these are all pure virtual so must be included here
  
  // reason: stiffness depends on quadrature point, not just element
  GEOSX_HOST_DEVICE inline
  virtual void GetStiffness( localIndex const k, real64 (& c)[6][6] ) const override final
  {
    GEOSX_UNUSED_VAR(k);
    GEOSX_UNUSED_VAR(c);
    GEOSX_ERROR("Not implemented");
  }
  
  // reason: DP-model requires state
  GEOSX_HOST_DEVICE
  virtual void SmallStrainNoState( localIndex const k,
                                   real64 const * const GEOSX_RESTRICT voigtStrain,
                                   real64 * const GEOSX_RESTRICT stress ) const override final;

  // this is possible to implement, but would just be a shortcut version without stiffness calculation
  GEOSX_HOST_DEVICE
  virtual void SmallStrain( localIndex const k,
                            localIndex const q,
                            real64 const * const GEOSX_RESTRICT voigtStrainIncrement ) const override final;

  // reason: need to fix small strain interface first
  GEOSX_HOST_DEVICE
  virtual void HypoElastic( localIndex const k,
                            localIndex const q,
                            real64 const * const GEOSX_RESTRICT Ddt,
                            R2Tensor const & Rot ) const override final;

  // reason: need to fix small strain interface first
  GEOSX_HOST_DEVICE
  virtual void HyperElastic( localIndex const k,
                             real64 const (&FmI)[3][3],
                             real64 * const GEOSX_RESTRICT stress ) const override final;

  // reason: need to fix small strain interface first
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
void DruckerPragerUpdates::SmallStrain( localIndex const k,
                                        localIndex const q,
                                        arraySlice1d< real64 const > const & strainIncrement,
                                        arraySlice1d< real64 > const & stress,
                                        arraySlice2d< real64 > const & stiffness )
{
  GEOSX_UNUSED_VAR(k);
  GEOSX_UNUSED_VAR(q);
  
  for(localIndex i=0; i<6; ++i)
  {
    stress[i] = strainIncrement[i];
  }
  
  GEOSX_UNUSED_VAR(stiffness);
}
  
// base class defines many pure virtual functions
// TODO: revisit if this is required

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::SmallStrainNoState( localIndex const k,
                                               real64 const * GEOSX_RESTRICT const voigtStrain,
                                               real64 * GEOSX_RESTRICT const stress ) const
{
  GEOSX_UNUSED_VAR(k);
  GEOSX_UNUSED_VAR(voigtStrain);
  GEOSX_UNUSED_VAR(stress);
  GEOSX_ERROR("Not implemented");
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::SmallStrain( localIndex const k,
                                        localIndex const q,
                                        real64 const * const GEOSX_RESTRICT voigtStrainInc ) const
{
  GEOSX_UNUSED_VAR(k);
  GEOSX_UNUSED_VAR(q);
  GEOSX_UNUSED_VAR(voigtStrainInc);
  GEOSX_ERROR("Not implemented");
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::HypoElastic( localIndex const k,
                                        localIndex const q,
                                        real64 const * const GEOSX_RESTRICT Ddt,
                                        R2Tensor const & Rot ) const
{
  GEOSX_UNUSED_VAR(k);
  GEOSX_UNUSED_VAR(q);
  GEOSX_UNUSED_VAR(Ddt);
  GEOSX_UNUSED_VAR(Rot);
  GEOSX_ERROR("Not implemented");
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::HyperElastic( localIndex const k,
                                         real64 const (&FmI)[3][3],
                                         real64 * const GEOSX_RESTRICT stress ) const
{
  GEOSX_UNUSED_VAR(k);
  GEOSX_UNUSED_VAR(FmI);
  GEOSX_UNUSED_VAR(stress);
  GEOSX_ERROR("Not implemented");
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DruckerPragerUpdates::HyperElastic( localIndex const k,
                                         localIndex const q,
                                         real64 const (&FmI)[3][3] ) const
{
  GEOSX_UNUSED_VAR(k);
  GEOSX_UNUSED_VAR(q);
  GEOSX_UNUSED_VAR(FmI);
  GEOSX_ERROR("Not implemented");
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
   * @brief Create a instantiation of the DruckerPragerUpdate class that refers to the data in this.
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
  //TODO: maybe define some helper structs for defaults, arrays, arrayViews, etc.
  
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
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DRUCKERPRAGER_HPP_ */
