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
                                real64 const & bulkModulusA,
                                real64 const & bulkModulusB,
                                real64 const & bulkModulusT0,
                                real64 const & shearModulus,
                                real64 const & shearModulusA,
                                real64 const & shearModulusB,
                                real64 const & shearModulusT0,
                                real64 const & yieldStrength0,
                                real64 const & yieldStrengthA,
                                real64 const & yieldStrengthB,
                                real64 const & yieldStrengthT0,
                                real64 const & r0,
                                real64 const & r0A,
                                real64 const & r0B,
                                real64 const & r0T0,
                                real64 const & r1,
                                real64 const & r2,
                                real64 const & Gr,
                                real64 const & GrA,
                                real64 const & GrB,
                                real64 const & GrT0,
                                real64 const & maxStretch,
                                real64 const & maxStretchA,
                                real64 const & maxStretchB,
                                real64 const & maxStretchT0,
                                arrayView1d< real64 > const & damage,
                                arrayView1d< real64 > const & temperature,
                                arrayView1d< real64 > const & previousLambda,
                                arrayView1d< real64 > const & previousPlasticStrain,
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
       m_bulkModulusA( bulkModulusA ),
       m_bulkModulusB( bulkModulusB ),
       m_bulkModulusT0( bulkModulusT0 ),
       m_shearModulus( shearModulus ),
       m_shearModulusA( shearModulusA ),
       m_shearModulusB( shearModulusB ),
       m_shearModulusT0( shearModulusT0 ),
       m_yieldStrength0( yieldStrength0 ),
       m_yieldStrengthA( yieldStrengthA ),
       m_yieldStrengthB( yieldStrengthB ),
       m_yieldStrengthT0( yieldStrengthT0 ),
       m_r0( r0 ),
       m_r0A( r0A ),
       m_r0B( r0B ),
       m_r0T0( r0T0 ),
       m_r1( r1 ),
       m_r2( r2 ),
       m_Gr( Gr ),
       m_GrA( GrA ),
       m_GrB( GrB ),
       m_GrT0( GrT0 ),
       m_maxStretch( maxStretch ),
       m_maxStretchA( maxStretchA ),
       m_maxStretchB( maxStretchB ),
       m_maxStretchT0( maxStretchT0 ),
       m_damage( damage ),
       m_temperature( temperature ),
       m_previousLambda( previousLambda ),
       m_previousPlasticStrain( previousPlasticStrain )
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

      real64 K = m_bulkModulus * thermalSoftening( m_temperature[k], m_bulkModulusT0, m_bulkModulusA, m_bulkModulusB );
      real64 g = m_shearModulus * thermalSoftening( m_temperature[k], m_shearModulusT0, m_shearModulusA, m_shearModulusB );
      real64 gSqr = g * g;

      real64 lambda = LvArray::math::sqrt( LvArray::math::pow( normalDisplacement + t, 2 ) + shearDisplacementSqr ) / t; // Instantaneous stretch

      real64 maxStretch = m_maxStretch * thermalSoftening( m_temperature[k], m_maxStretchT0, m_maxStretchA, m_maxStretchB );
      if( lambda > maxStretch )
      {
          m_damage[k] = LvArray::math::max( 1.0, m_damage[k] );
      }

      real64 Gr = m_Gr * thermalSoftening( m_temperature[k], m_GrT0, m_GrA, m_GrB );
      real64 sigma_H = Gr * LvArray::math::pow( lambda * lambda - 1/lambda, 2 );

      real64 tau = LvArray::math::sqrt( gSqr*(4*t*normalDisplacement*(normalDisplacementSqr + shearDisplacementSqr) +
                    LvArray::math::pow( normalDisplacementSqr+shearDisplacementSqr, 2 )+tSqr * (4*normalDisplacementSqr+3*shearDisplacementSqr))/(LvArray::math::pow( t*(t+normalDisplacement), 2 )));

      real64 gamma_p = LvArray::math::max( 0.0, (tau-(m_yieldStrength0+sigma_H))/(K + (4.0/3.0)*g) );

      real64 r0 = m_r0 * thermalSoftening( m_temperature[k], m_r0T0, m_r0A, m_r0B );
      real64 R_gamma = r0*LvArray::math::exp( -LvArray::math::pow( gamma_p/m_r1, m_r2 ));

      real64 yieldStrength = m_yieldStrength0 * thermalSoftening( m_temperature[k], m_yieldStrengthT0, m_yieldStrengthA, m_yieldStrengthB ) + R_gamma + sigma_H;

      real64 scale  = 1.0;
      if( tau > yieldStrength )
      {
          scale = yieldStrength / tau;
      }

      normalStress = -scale * (normalDisplacement*(K*(t+normalDisplacement)+g*(2*t+normalDisplacement)))/(t*(t+normalDisplacement));
      shearStress = -scale * g * tangentialDisplacement / t;
    }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 thermalSoftening( const real64 & T,
                           const real64 & T0,
                           const real64 & A,
                           const real64 & B
                           ) const
  {
    if( LvArray::math::abs( A ) > 1.e-16 )
    {
      return 1.0 + A / (1.0 + LvArray::math::exp( B * (T-T0) ) );
    }
    else
    {
      return 1.0;
    }
  }


private:

  // Constants
  real64 m_thickness;

  real64 m_bulkModulus;
  real64 m_bulkModulusA;
  real64 m_bulkModulusB;
  real64 m_bulkModulusT0;

  real64 m_shearModulus;
  real64 m_shearModulusA;
  real64 m_shearModulusB;
  real64 m_shearModulusT0;

  real64 m_yieldStrength0;
  real64 m_yieldStrengthA;
  real64 m_yieldStrengthB;
  real64 m_yieldStrengthT0;
  
  real64 m_r0;
  real64 m_r0A;
  real64 m_r0B;
  real64 m_r0T0;

  real64 m_r1;
  real64 m_r2;

  real64 m_Gr;
  real64 m_GrA;
  real64 m_GrB;
  real64 m_GrT0;
  
  real64 m_maxStretch;
  real64 m_maxStretchA;
  real64 m_maxStretchB;
  real64 m_maxStretchT0;

  /// A reference to the current damage at each cohesive zone node.
  arrayView1d< real64 > const m_damage;

  /// A reference to the current temperature at each cohesive zone node
  arrayView1d< real64 > const m_temperature;

  arrayView1d< real64 > const m_previousLambda;
  arrayView1d< real64 > const m_previousPlasticStrain;
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
    static constexpr char const * bulkModulusAString() { return "bulkModulusA"; }
    static constexpr char const * bulkModulusBString() { return "bulkModulusB"; }
    static constexpr char const * bulkModulusT0String() { return "bulkModulusT0"; }

    static constexpr char const * shearModulusString() { return "shearModulus"; }
    static constexpr char const * shearModulusAString() { return "shearModulusA"; }
    static constexpr char const * shearModulusBString() { return "shearModulusB"; }
    static constexpr char const * shearModulusT0String() { return "shearModulusT0"; }

    static constexpr char const * yieldStrength0String() { return "yieldStrength0"; }
    static constexpr char const * yieldStrengthAString() { return "yieldStrengthA"; }
    static constexpr char const * yieldStrengthBString() { return "yieldStrengthB"; }
    static constexpr char const * yieldStrengthT0String() { return "yieldStrengthT0"; }

    static constexpr char const * r0String() { return "r0"; }
    static constexpr char const * r0AString() { return "r0A"; }
    static constexpr char const * r0BString() { return "r0B"; }
    static constexpr char const * r0T0String() { return "r0T0"; }

    static constexpr char const * r1String() { return "r1"; }
    static constexpr char const * r2String() { return "r2"; }
    
    static constexpr char const * GrString() { return "Gr"; }
    static constexpr char const * GrAString() { return "GrA"; }
    static constexpr char const * GrBString() { return "GrB"; }
    static constexpr char const * GrT0String() { return "GrT0"; }
    

    static constexpr char const * maxStretchString() { return "maximumStretch"; }
    static constexpr char const * maxStretchAString() { return "maximumStretchA"; }
    static constexpr char const * maxStretchBString() { return "maximumStretchB"; }
    static constexpr char const * maxStretchT0String() { return "maximumStretchT0"; }

    static constexpr char const * damageString() { return "damage"; }
    static constexpr char const * temperatureString() { return "temperature"; }
    static constexpr char const * previousLambdaString() { return "previousLambda"; }
    static constexpr char const * previousPlasticStrainString() { return "previousPlasticStrain"; }
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
                                       m_bulkModulusA,
                                       m_bulkModulusB,
                                       m_bulkModulusT0,
                                       m_shearModulus,
                                       m_shearModulusA,
                                       m_shearModulusB,
                                       m_shearModulusT0,
                                       m_yieldStrength0,
                                       m_yieldStrengthA,
                                       m_yieldStrengthB,
                                       m_yieldStrengthT0,
                                       m_r0,
                                       m_r0A,
                                       m_r0B,
                                       m_r0T0,
                                       m_r1,
                                       m_r2,
                                       m_Gr,
                                       m_GrA,
                                       m_GrB,
                                       m_GrT0,
                                       m_maxStretch,
                                       m_maxStretchA,
                                       m_maxStretchB,
                                       m_maxStretchT0,
                                       m_damage,
                                       m_temperature,
                                       m_previousLambda,
                                       m_previousPlasticStrain,
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
                          m_bulkModulusA,
                          m_bulkModulusB,
                          m_bulkModulusT0,
                          m_shearModulus,
                          m_shearModulusA,
                          m_shearModulusB,
                          m_shearModulusT0,
                          m_yieldStrength0,
                          m_yieldStrengthA,
                          m_yieldStrengthB,
                          m_yieldStrengthT0,
                          m_r0,
                          m_r0A,
                          m_r0B,
                          m_r0T0,
                          m_r1,
                          m_r2,
                          m_Gr,
                          m_GrA,
                          m_GrB,
                          m_GrT0,
                          m_maxStretch,
                          m_maxStretchA,
                          m_maxStretchB,
                          m_maxStretchT0,
                          m_damage,
                          m_temperature,
                          m_previousLambda,
                          m_previousPlasticStrain,
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
  real64 m_bulkModulusA;
  real64 m_bulkModulusB;
  real64 m_bulkModulusT0;
  
  real64 m_shearModulus;
  real64 m_shearModulusA;
  real64 m_shearModulusB;
  real64 m_shearModulusT0;

  real64 m_yieldStrength0;
  real64 m_yieldStrengthA;
  real64 m_yieldStrengthB;
  real64 m_yieldStrengthT0;
  
  real64 m_r0;
  real64 m_r0A;
  real64 m_r0B;
  real64 m_r0T0;

  real64 m_r1;
  real64 m_r2;
  
  real64 m_Gr;
  real64 m_GrA;
  real64 m_GrB;
  real64 m_GrT0;
  
  real64 m_maxStretch;
  real64 m_maxStretchA;
  real64 m_maxStretchB;
  real64 m_maxStretchT0;

  // Stores the damage for each node in cohesive zone
  array1d< real64 > m_damage;
  array1d< real64 > m_temperature;
  array1d< real64 > m_previousLambda;
  array1d< real64 > m_previousPlasticStrain;
};

} // namespace constitutive
} // namespace geos

#endif /* GEOS_CONSTITUTIVE_POLYMERCOHESIVEZONE_HPP_ */
