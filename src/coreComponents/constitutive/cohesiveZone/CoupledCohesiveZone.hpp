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
 * @file CoupledCohesiveZone.hpp
 */

#ifndef GEOS_CONSTITUTIVE_COUPLEDCOHESIVEZONE_HPP_
#define GEOS_CONSTITUTIVE_COUPLEDCOHESIVEZONE_HPP_

#include "constitutive/cohesiveZone/CohesiveZoneBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{


class CoupledCohesiveZoneUpdates : public CohesiveZoneBaseUpdates
{
public:
    /**
    * @brief constructor
    * @param[in] newTraction The new traction data from the cohesive zone model class.
    * @param[in] oldTraction The old traction data from the cohesive zone model class.
    * @param[in] damage The damage data from the cohesive zone model class
    */
    CoupledCohesiveZoneUpdates( real64 const & characteristicNormalDisplacement,
                                real64 const & characteristicTangentialDisplacement,
                                real64 const & defaultMaxNormalStress,
                                real64 const & defaultMaxShearStress,
                                real64 const & maxNormalDisplacement,
                                real64 const & maxTangentialDisplacement,
                                arrayView1d< real64 > const & maxNormalStress,
                                arrayView1d< real64 > const & maxShearStress,
                                arrayView1d< real64 > const & damage,
                                arrayView2d< real64 > const & newNormalStress,
                                arrayView2d< real64 > const & newShearStress,
                                arrayView2d< real64 > const & oldNormalStress,
                                arrayView2d< real64 > const & oldShearStress ):
        CohesiveZoneBaseUpdates( newNormalStress,
                                 newShearStress,
                                 oldNormalStress,
                                 oldShearStress ),
        m_characteristicNormalDisplacement( characteristicNormalDisplacement  ),
        m_characteristicTangentialDisplacement( characteristicTangentialDisplacement ),
        m_defaultMaxNormalStress( defaultMaxNormalStress ),
        m_defaultMaxShearStress( defaultMaxShearStress ),
        m_maxNormalDisplacement( maxNormalDisplacement ),
        m_maxTangentialDisplacement( maxTangentialDisplacement ),
        m_maxNormalStress( maxNormalStress ),
        m_maxShearStress( maxShearStress ),
        m_damage( damage )
    {}

    /// Deleted default constructor
    CoupledCohesiveZoneUpdates() = delete;

    /**
    * @brief Copy Constructor
    * @param source Object to copy
    */
    CoupledCohesiveZoneUpdates( CoupledCohesiveZoneUpdates const & source ) = default;

    /**
    * @brief Move Constructor
    * @param source Object to move resources from
    */
    CoupledCohesiveZoneUpdates( CoupledCohesiveZoneUpdates && source ) = default;

    /// Deleted copy assignment operator
    CoupledCohesiveZoneUpdates & operator=( CoupledCohesiveZoneUpdates const & ) = delete;

    /// Deleted move assignment operator
    CoupledCohesiveZoneUpdates & operator=( CoupledCohesiveZoneUpdates && ) =  delete;

    GEOS_HOST_DEVICE
    void jumpDisplacementUpdate( localIndex const k,
                                 real64 const & normalDisplacement,
                                 real64 const & tangentialDisplacement,
                                 real64 & normalStress,
                                 real64 & shearStress ) const
    {
        if( ( normalDisplacement >= m_maxNormalDisplacement ) ||
            ( tangentialDisplacement >= m_maxTangentialDisplacement ) )
        {
            m_damage[k] = LvArray::math::max( 1.0, m_damage[k] );
        }

        real64 e = 2.7182818284590452353602874713526624977572470936999; // Euler's number
        real64 r = 0.0;
        real64 q = 1.0;
   
        real64 normalWorkOfSeparation = e * m_maxNormalStress[k] * m_characteristicNormalDisplacement;
        real64 shearWorkOfSeparation = LvArray::math::sqrt( e / 2 ) * m_maxShearStress[k] * m_characteristicTangentialDisplacement;

        real64 normalizedNormalDisplacement = normalDisplacement / m_characteristicNormalDisplacement;
        real64 normalizedTangentialDisplacement = tangentialDisplacement/ m_characteristicTangentialDisplacement;

        normalStress = -( normalWorkOfSeparation / m_characteristicNormalDisplacement ) * exp( -normalizedNormalDisplacement ) *
                        ( normalizedNormalDisplacement *
                        LvArray::math::exp( -LvArray::math::pow( normalizedTangentialDisplacement,
                                    2 ) ) + ( 1.0 - q ) / ( r - 1.0 ) * ( 1.0 - LvArray::math::exp( -pow( normalizedTangentialDisplacement, 2 ) ) ) * ( r - normalizedNormalDisplacement ) );
        shearStress = -( shearWorkOfSeparation / m_characteristicTangentialDisplacement ) * ( 2.0 * m_characteristicNormalDisplacement / m_characteristicTangentialDisplacement ) * normalizedTangentialDisplacement *
                        ( q + ( r - q ) / ( r - 1.0 ) * normalizedNormalDisplacement ) *
                        LvArray::math::exp( -normalizedNormalDisplacement ) * LvArray::math::exp( -LvArray::math::pow( normalizedTangentialDisplacement, 2 ) );
    }

protected:

  // Constants
  real64 const m_characteristicNormalDisplacement;
  real64 const m_characteristicTangentialDisplacement;
  real64 const m_defaultMaxNormalStress;
  real64 const m_defaultMaxShearStress;
  real64 const m_maxNormalDisplacement;
  real64 const m_maxTangentialDisplacement;

  arrayView1d< real64 > const m_maxNormalStress;
  arrayView1d< real64 > const m_maxShearStress;

  /// A reference the current damage at each cohesive zone node.
  arrayView1d< real64 > const m_damage;
};


/**
 * @class CoupledCohesiveZone
 * This class serves as the base class for cohesive zone models.
 */
class CoupledCohesiveZone : public CohesiveZoneBase
{
public:
  /// @typedef Alias for CoupledCohesiveZoneUpdates
  using KernelWrapper = CoupledCohesiveZoneUpdates;

  /**
   * @brief Constructor
   * @param name Name of the CoupledCohesiveZone object in the repository.
   * @param parent The parent group of the CoupledCohesiveZone object.
   */
  CoupledCohesiveZone( string const & name,
                         Group * const parent );

  /**
   * Destructor
   */
  virtual ~CoupledCohesiveZone() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "CoupledCohesiveZone";

  /**
   * @brief Static catalog string
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// Keys for data in this class
  struct viewKeyStruct : public CohesiveZoneBase::viewKeyStruct
  {
    static constexpr char const * characteristicNormalDisplacementString() { return "characteristicNormalDisplacement"; }
    static constexpr char const * characteristicTangentialDisplacementString() { return "characteristicTangentialDisplacement"; }
    static constexpr char const * defaultMaxNormalStressString() { return "defaultMaxNormalStress"; }
    static constexpr char const * defaultMaxShearStressString() { return "defaultMaxShearStress"; }
    static constexpr char const * maxNormalDisplacementString() { return "maxNormalDisplacement"; }
    static constexpr char const * maxTangentialDisplacementString() { return "maxTangentialDisplacement"; }
    static constexpr char const * maxNormalStressString() { return "maxNormalStress"; }
    static constexpr char const * maxShearStressString() { return "maxShearStress"; }
    static constexpr char const * damageString() { return "damage"; }
  };

  /**
   * @name Accessors
   */
  ///@{

  // TODO add getters for constants
//   real64 getMaxNormalStress() const
//   {
//     return m_normalForceConstant;
//   }

//   real64 getShearForceConstant() const
//   {
//     return m_shearForceConstant;
//   }

  /**
   * @brief Non-const/mutable accessor for damage
   * @return Accessor
   */
  arrayView1d< real64 > const getDamage()
  {
    return m_damage;
  }

  /**
   * @brief Const/non-mutable accessor for damage
   * @return Accessor
   */
  arrayView1d< real64 const > const getDamage() const
  {
    return m_damage;
  }

  ///@}

  /**
   * @brief Create a instantiation of the CoupledCohesiveZoneUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of CoupledCohesiveZoneUpdate.
   */
  CoupledCohesiveZoneUpdates createKernelUpdates( bool const includeState = true ) const
  {
    GEOS_UNUSED_VAR( includeState );
    return CoupledCohesiveZoneUpdates( m_characteristicNormalDisplacement,
                                       m_characteristicTangentialDisplacement,
                                       m_defaultMaxNormalStress,
                                       m_defaultMaxShearStress,
                                       m_maxNormalDisplacement,
                                       m_maxTangentialDisplacement,
                                       m_maxNormalStress,
                                       m_maxShearStress,
                                       m_damage,
                                       m_newNormalStress,
                                       m_newShearStress,
                                       m_oldNormalStress,
                                       m_oldShearStress );
  } 

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for
   *   the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams ) const
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_characteristicNormalDisplacement,
                          m_characteristicTangentialDisplacement,
                          m_defaultMaxNormalStress,
                          m_defaultMaxShearStress,
                          m_maxNormalDisplacement,
                          m_maxTangentialDisplacement,
                          m_maxNormalStress,
                          m_maxShearStress,
                          m_damage,
                          m_newNormalStress,
                          m_newShearStress,
                          m_oldNormalStress,
                          m_oldShearStress );
  }

protected:

  /// Post-process XML input
  virtual void postInputInitialization() override;

  // Constants
  real64 m_characteristicNormalDisplacement;
  real64 m_characteristicTangentialDisplacement;
  real64 m_defaultMaxNormalStress;
  real64 m_defaultMaxShearStress;
  real64 m_maxNormalDisplacement;
  real64 m_maxTangentialDisplacement;

  // Stores the damage for each node in cohesive zone
  array1d< real64 > m_maxNormalStress;
  array1d< real64 > m_maxShearStress;
  array1d< real64 > m_damage;

};

} // namespace constitutive
} // namespace geos

#endif /* GEOS_CONSTITUTIVE_COUPLEDCOHESIVEZONE_HPP_ */
