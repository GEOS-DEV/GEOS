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
 * @file UncoupledCohesiveZone.hpp
 */

#ifndef GEOS_CONSTITUTIVE_UNCOUPLEDCOHESIVEZONE_HPP_
#define GEOS_CONSTITUTIVE_UNCOUPLEDCOHESIVEZONE_HPP_

#include "constitutive/cohesiveZone/CohesiveZoneBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

class UncoupledCohesiveZoneUpdates : public CohesiveZoneBaseUpdates
{
public:
    /**
    * @brief constructor
    * @param[in] newTraction The new traction data from the cohesive zone model class.
    * @param[in] oldTraction The old traction data from the cohesive zone model class.
    * @param[in] damage The damage data from the cohesive zone model class
    */
    UncoupledCohesiveZoneUpdates( real64 const & normalForceConstant,
                                  real64 const & shearForceConstant,
                                  real64 const & maxNormalDisplacement,
                                  real64 const & maxTangentialDisplacement,
                                  arrayView1d< real64 > const & damage,
                                  arrayView2d< real64 > const & newNormalStress,
                                  arrayView2d< real64 > const & newShearStress,
                                  arrayView2d< real64 > const & oldNormalStress,
                                  arrayView2d< real64 > const & oldShearStress ):
        CohesiveZoneBaseUpdates( newNormalStress,
                                 newShearStress,
                                 oldNormalStress,
                                 oldShearStress ),
        m_normalForceConstant( normalForceConstant ),
        m_shearForceConstant( shearForceConstant ),
        m_maxNormalDisplacement( maxNormalDisplacement ),
        m_maxTangentialDisplacement( maxTangentialDisplacement ),
        m_damage( damage )
    {}

    /// Deleted default constructor
    UncoupledCohesiveZoneUpdates() = delete;

    /**
    * @brief Copy Constructor
    * @param source Object to copy
    */
    UncoupledCohesiveZoneUpdates( UncoupledCohesiveZoneUpdates const & source ) = default;

    /**
    * @brief Move Constructor
    * @param source Object to move resources from
    */
    UncoupledCohesiveZoneUpdates( UncoupledCohesiveZoneUpdates && source ) = default;

    /// Deleted copy assignment operator
    UncoupledCohesiveZoneUpdates & operator=( UncoupledCohesiveZoneUpdates const & ) = delete;

    /// Deleted move assignment operator
    UncoupledCohesiveZoneUpdates & operator=( UncoupledCohesiveZoneUpdates && ) =  delete;

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
            m_damage[k] = 1.0;
        }

        normalStress = -m_normalForceConstant * normalDisplacement * (1.0 - m_damage[k]);
        shearStress = -m_shearForceConstant * tangentialDisplacement * (1.0 - m_damage[k]);
    }

private:

  // Constants
  real64 const m_normalForceConstant;
  real64 const m_shearForceConstant;
  real64 const m_maxNormalDisplacement;
  real64 const m_maxTangentialDisplacement;

  /// A reference the current damage at each cohesive zone node.
  arrayView1d< real64 > const m_damage;
};


/**
 * @class UncoupledCohesiveZone
 * This class serves as the base class for cohesive zone models.
 */
class UncoupledCohesiveZone : public CohesiveZoneBase
{
public:
  /// @typedef Alias for UncoupledCohesiveZoneUpdates
  using KernelWrapper = UncoupledCohesiveZoneUpdates;

  /**
   * @brief Constructor
   * @param name Name of the UncoupledCohesiveZone object in the repository.
   * @param parent The parent group of the UncoupledCohesiveZone object.
   */
  UncoupledCohesiveZone( string const & name,
                         Group * const parent );

  /**
   * Destructor
   */
  virtual ~UncoupledCohesiveZone() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "UncoupledCohesiveZone";

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
    static constexpr char const * normalForceConstantString() { return "normalForceConstant"; }
    static constexpr char const * shearForceConstantString() { return "shearForceConstant"; }
    static constexpr char const * maxNormalDisplacementString() { return "maxNormalDisplacement"; }
    static constexpr char const * maxTangentialDisplacementString() { return "maxTangentialDisplacement"; }
    static constexpr char const * damageString() { return "damage"; }
  };

  /**
   * @name Accessors
   */
  ///@{

  real64 getNormalForceConstant() const
  {
    return m_normalForceConstant;
  }

  real64 getShearForceConstant() const
  {
    return m_shearForceConstant;
  }

  real64 getMaxNormalDisplacement() const
  {
    return m_maxNormalDisplacement;
  }

  real64 getMaxTangentialDisplacement() const
  {
    return m_maxTangentialDisplacement;
  }

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
   * @brief Create a instantiation of the UncoupledCohesiveZoneUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of UncoupledCohesiveZoneUpdate.
   */
  UncoupledCohesiveZoneUpdates createKernelUpdates( bool const includeState = true ) const
  {
    GEOS_UNUSED_VAR( includeState );
    return UncoupledCohesiveZoneUpdates( m_normalForceConstant,
                                         m_shearForceConstant,
                                         m_maxNormalDisplacement,
                                         m_maxTangentialDisplacement,
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
                          m_normalForceConstant,
                          m_shearForceConstant,
                          m_maxNormalDisplacement,
                          m_maxTangentialDisplacement,
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
  real64 m_normalForceConstant;
  real64 m_shearForceConstant;
  real64 m_maxNormalDisplacement;
  real64 m_maxTangentialDisplacement;

  // Stores the damage for each node in cohesive zone
  array1d< real64 > m_damage;

};

} // namespace constitutive
} // namespace geos

#endif /* GEOS_CONSTITUTIVE_UNCOUPLEDCOHESIVEZONE_HPP_ */
