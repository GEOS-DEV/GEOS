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
 * @file PolymerCohesiveZone.hpp
 */

#ifndef GEOS_CONSTITUTIVE_POLYMERCOHESIVEZONE_HPP_
#define GEOS_CONSTITUTIVE_POLYMERCOHESIVEZONE_HPP_

#include "constitutive/cohesiveZone/CohesiveZoneBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{


class PolymerCohesiveZoneUpdates : public CohesiveZoneBaseUpdates
{
public:
    /**
    * @brief constructor
    * @param[in] newTraction The new traction data from the cohesive zone model class.
    * @param[in] oldTraction The old traction data from the cohesive zone model class.
    * @param[in] damage The damage data from the cohesive zone model class
    */
    PolymerCohesiveZoneUpdates( real64 const & thickness,
                                real64 const & bulkModulus,
                                real64 const & shearModulus,
                                real64 const & yieldStrength0,
                                real64 const & r0,
                                real64 const & r1,
                                real64 const & r2,
                                real64 const & Gr,
                                real64 const & maxStretch,
                                arrayView1d< real64 > const & damage,
                                arrayView2d< real64 > const & newNormalStress,
                                arrayView2d< real64 > const & newShearStress,
                                arrayView2d< real64 > const & oldNormalStress,
                                arrayView2d< real64 > const & oldShearStress ):
        CohesiveZoneBaseUpdates( newNormalStress,
                                 newShearStress,
                                 oldNormalStress,
                                 oldShearStress ),
       m_thickness( thickness ),
       m_bulkModulus( bulkModulus ),
       m_shearModulus( shearModulus ),
       m_yieldStrength0( yieldStrength0 ),
       m_r0( r0 ),
       m_r1( r1 ),
       m_r2( r2 ),
       m_Gr( Gr ),
       m_maxStretch( maxStretch ),
       m_damage( damage )
    {}

    /// Deleted default constructor
    PolymerCohesiveZoneUpdates() = delete;

    /**
    * @brief Copy Constructor
    * @param source Object to copy
    */
    PolymerCohesiveZoneUpdates( PolymerCohesiveZoneUpdates const & source ) = default;

    /**
    * @brief Move Constructor
    * @param source Object to move resources from
    */
    PolymerCohesiveZoneUpdates( PolymerCohesiveZoneUpdates && source ) = default;

    /// Deleted copy assignment operator
    PolymerCohesiveZoneUpdates & operator=( PolymerCohesiveZoneUpdates const & ) = delete;

    /// Deleted move assignment operator
    PolymerCohesiveZoneUpdates & operator=( PolymerCohesiveZoneUpdates && ) =  delete;

    GEOS_HOST_DEVICE
    void jumpDisplacementUpdate( localIndex const k,
                                 real64 const & normalDisplacement,
                                 real64 const & tangentialDisplacement,
                                 real64 & normalStress,
                                 real64 & shearStress ) const
    {
        real64 t = m_thickness;
        real64 tSqr = t * t;
        real64 normalDisplacementSqr = normalDisplacement * normalDisplacement;
        real64 shearDisplacementSqr = tangentialDisplacement * tangentialDisplacement;

        real64 K = m_bulkModulus;
        real64 g = m_shearModulus;
        real64 gSqr = g * g;

        real64 lambda = LvArray::math::sqrt( LvArray::math::pow( normalDisplacement + t, 2 ) + shearDisplacementSqr ) / t; // Instantaneous stretch

        if( lambda > m_maxStretch )
        {
            m_damage[k] = LvArray::math::max( 1.0, m_damage[k] );
        }

        real64 sigma_H = m_Gr * LvArray::math::pow( lambda * lambda - 1/lambda, 2 );

        real64 tau = LvArray::math::sqrt( gSqr*(4*t*normalDisplacement*(normalDisplacementSqr + shearDisplacementSqr) +
                     LvArray::math::pow( normalDisplacementSqr+shearDisplacementSqr, 2 )+tSqr * (4*normalDisplacementSqr+3*shearDisplacementSqr))/(LvArray::math::pow( t*(t+normalDisplacement), 2 )));

        real64 gamma_p = LvArray::math::max( 0.0, (tau-(m_yieldStrength0+sigma_H))/(K + (4.0/3.0)*g) );

        real64 R_gamma = m_r0*LvArray::math::exp( -LvArray::math::pow( gamma_p/m_r1, m_r2 ));

        real64 yieldStrength = m_yieldStrength0 + R_gamma + sigma_H;

        real64 scale  = 1.0;
        if( tau > yieldStrength )
        {
            scale = yieldStrength / tau;
        }

        normalStress = -scale * (normalDisplacement*(K*(t+normalDisplacement)+g*(2*t+normalDisplacement)))/(t*(t+normalDisplacement));
        shearStress = -g * tangentialDisplacement / t;
    }

private:

  // Constants
  real64 const m_thickness;
  real64 const m_bulkModulus;
  real64 const m_shearModulus;
  real64 const m_yieldStrength0;
  real64 const m_r0;
  real64 const m_r1;
  real64 const m_r2;
  real64 const m_Gr;
  real64 const m_maxStretch;

  /// A reference the current damage at each cohesive zone node.
  arrayView1d< real64 > const m_damage;
};


/**
 * @class PolymerCohesiveZone
 * This class serves as the base class for cohesive zone models.
 */
class PolymerCohesiveZone : public CohesiveZoneBase
{
public:
  /// @typedef Alias for PolymerCohesiveZoneUpdates
  using KernelWrapper = PolymerCohesiveZoneUpdates;

  /**
   * @brief Constructor
   * @param name Name of the PolymerCohesiveZone object in the repository.
   * @param parent The parent group of the PolymerCohesiveZone object.
   */
  PolymerCohesiveZone( string const & name,
                         Group * const parent );

  /**
   * Destructor
   */
  virtual ~PolymerCohesiveZone() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "PolymerCohesiveZone";

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
    static constexpr char const * thicknessString() { return "thickness"; }
    static constexpr char const * bulkModulusString() { return "bulkModulus"; }
    static constexpr char const * shearModulusString() { return "shearModulus"; }
    static constexpr char const * yieldStrength0String() { return "yieldStrength0"; }
    static constexpr char const * r0String() { return "r0"; }
    static constexpr char const * r1String() { return "r1"; }
    static constexpr char const * r2String() { return "r2"; }
    static constexpr char const * GrString() { return "Gr"; }
    static constexpr char const * maxStretchString() { return "maxStretch"; }
    static constexpr char const * damageString() { return "damage"; }
  };

  /**
   * @name Accessors
   */
  ///@{

  // TODO add getters for constants

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
   * @brief Create a instantiation of the PolymerCohesiveZoneUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of PolymerCohesiveZoneUpdate.
   */
  PolymerCohesiveZoneUpdates createKernelUpdates( bool const includeState = true ) const
  {
    GEOS_UNUSED_VAR( includeState );
    return PolymerCohesiveZoneUpdates( m_thickness,
                                       m_bulkModulus,
                                       m_shearModulus,
                                       m_yieldStrength0,
                                       m_r0,
                                       m_r1,
                                       m_r2,
                                       m_Gr,
                                       m_maxStretch,
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
                          m_thickness,
                          m_bulkModulus,
                          m_shearModulus,
                          m_yieldStrength0,
                          m_r0,
                          m_r1,
                          m_r2,
                          m_Gr,
                          m_maxStretch,
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
  real64 m_thickness;
  real64 m_bulkModulus;
  real64 m_shearModulus;
  real64 m_yieldStrength0;
  real64 m_r0;
  real64 m_r1;
  real64 m_r2;
  real64 m_Gr;
  real64 m_maxStretch;

  // Stores the damage for each node in cohesive zone
  array1d< real64 > m_damage;

};

} // namespace constitutive
} // namespace geos

#endif /* GEOS_CONSTITUTIVE_POLYMERCOHESIVEZONE_HPP_ */
