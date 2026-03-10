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
 * @file BicrystalCohesiveZone.hpp
 */

#ifndef GEOS_CONSTITUTIVE_BICRYSTALCOHESIVEZONE_HPP_
#define GEOS_CONSTITUTIVE_BICRYSTALCOHESIVEZONE_HPP_

#include "constitutive/cohesiveZone/CoupledCohesiveZone.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{


class BicrystalCohesiveZoneUpdates : public CoupledCohesiveZoneUpdates
{
public:
    /**
    * @brief constructor
    * @param[in] newTraction The new traction data from the cohesive zone model class.
    * @param[in] oldTraction The old traction data from the cohesive zone model class.
    * @param[in] damage The damage data from the cohesive zone model class
    */
    BicrystalCohesiveZoneUpdates( arrayView3d < real64 > const & misorientation, 
                                  real64 const & characteristicNormalDisplacement,
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
        CoupledCohesiveZoneUpdates( characteristicNormalDisplacement,
                                    characteristicTangentialDisplacement,
                                    defaultMaxNormalStress,
                                    defaultMaxShearStress,
                                    maxNormalDisplacement,
                                    maxTangentialDisplacement,
                                    maxNormalStress,
                                    maxShearStress,
                                    damage,
                                    newNormalStress,
                                    newShearStress,
                                    oldNormalStress,
                                    oldShearStress ),
        m_misorientation( misorientation )
    {}

    /// Deleted default constructor
    BicrystalCohesiveZoneUpdates() = delete;

    /**
    * @brief Copy Constructor
    * @param source Object to copy
    */
    BicrystalCohesiveZoneUpdates( BicrystalCohesiveZoneUpdates const & source ) = default;

    /**
    * @brief Move Constructor
    * @param source Object to move resources from
    */
    BicrystalCohesiveZoneUpdates( BicrystalCohesiveZoneUpdates && source ) = default;

    /// Deleted copy assignment operator
    BicrystalCohesiveZoneUpdates & operator=( BicrystalCohesiveZoneUpdates const & ) = delete;

    /// Deleted move assignment operator
    BicrystalCohesiveZoneUpdates & operator=( BicrystalCohesiveZoneUpdates && ) =  delete;

    GEOS_HOST_DEVICE
    void jumpDisplacementUpdate( localIndex const k,
                                 real64 const & normalDisplacement,
                                 real64 const & tangentialDisplacement,
                                 real64 & normalStress,
                                 real64 & shearStress ) const
    {
      // Compute Euler angles from misorientation
      // X-basis denotes c-axis of crystal
      real64 theta = LvArray::math::acos( m_misorientation[k][0][0] ); // Returns 0 - pi, but should be symmetric about 90 deg
      // real64 phi = LvArray::math::acos( -m_misorientation[k][0][1] * LvArray::math::sin(theta) + m_misorientation[k][1][1] * LvArray::math::cos(theta) );
      
      real64 pi = 3.14159265358979;
      if( theta > pi/2 )
      {
        theta = pi - theta;
      }
      // Clamp theta to we don't get numerical issues
      theta = LvArray::math::max( LvArray::math::min(theta, pi/2 ), 0.0 );


      // Assumes simple Read-Shockley model of interface
      real64 knockdown = 0.0;
      if (!isZero(theta))
      {
        real64 theta_m = 15.0 * pi / 180.0; // Limiting misorientation angle 
        if( LvArray::math::abs(theta) < theta_m )
        {
          real64 theta_p = theta / theta_m; 
          knockdown = theta_p * ( 1.0 - LvArray::math::log(theta_p) );
        }
        else
        {
          knockdown = 1.0;
        }
      }
      
      real64 knockdown_limit = 0.75; // Caps knockdown for highly misoriented interfaces so they still have some strength
      real64 scale = 1.0 - knockdown_limit * knockdown;
      m_maxNormalStress[k] = m_defaultMaxNormalStress * scale;
      m_maxShearStress[k] = m_defaultMaxShearStress * scale;

      CoupledCohesiveZoneUpdates::jumpDisplacementUpdate( k, 
                                                          normalDisplacement,
                                                          tangentialDisplacement, 
                                                          normalStress,
                                                          shearStress );
    } 

private:

  // Constants
  // None currently...

  // A reference the misorientation at each cohesive zone node.
  arrayView3d< real64 > const m_misorientation;
};


/**
 * @class BicrystalCohesiveZone
 * This class serves as the base class for cohesive zone models.
 */
class BicrystalCohesiveZone : public CoupledCohesiveZone
{
public:
  /// @typedef Alias for BicrystalCohesiveZoneUpdates
  using KernelWrapper = BicrystalCohesiveZoneUpdates;

  /**
   * @brief Constructor
   * @param name Name of the BicrystalCohesiveZone object in the repository.
   * @param parent The parent group of the BicrystalCohesiveZone object.
   */
  BicrystalCohesiveZone( string const & name,
                         Group * const parent );

  /**
   * Destructor
   */
  virtual ~BicrystalCohesiveZone() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "BicrystalCohesiveZone";

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
  struct viewKeyStruct : public CoupledCohesiveZone::viewKeyStruct
  {
    static constexpr char const * misorientationString() { return "misorientation"; }
  };

  /**
   * @name Accessors
   */
  ///@{

  /**
   * @brief Non-const/mutable accessor for misorientation
   * @return Accessor
   */
  arrayView3d< real64 > const getMisorientation()
  {
    return m_misorientation;
  }

  /**
   * @brief Const/non-mutable accessor for misorientation
   * @return Accessor
   */
  arrayView3d< real64 const > const getMisorientation() const
  {
    return m_misorientation;
  }

  ///@}

  /**
   * @brief Create a instantiation of the BicrystalCohesiveZoneUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of BicrystalCohesiveZoneUpdate.
   */
  BicrystalCohesiveZoneUpdates createKernelUpdates( bool const includeState = true ) const
  {
    GEOS_UNUSED_VAR( includeState );
    return BicrystalCohesiveZoneUpdates( m_misorientation,
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
                          m_misorientation,
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
  // real64 m_characteristicNormalDisplacement;

  // Stores the misorientation as a rotational transformation for each point
  array3d< real64 > m_misorientation;
};

} // namespace constitutive
} // namespace geos

#endif /* GEOS_CONSTITUTIVE_BICRYSTALCOHESIVEZONE_HPP_ */
