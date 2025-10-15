/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_OEDOMETRICSTRESSPATH_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_OEDOMETRICSTRESSPATH_HPP_

#include "constitutive/prescribedStressPathGeomechanics/PrescribedStressPathGeomechanicsBase.hpp"

namespace geos
{
namespace constitutive
{

class OedometricStressPathUpdates
{
public:

  OedometricStressPathUpdates( real64 const aperture0 )
  {}

  /// Default copy constructor
  OedometricStressPathUpdates( OedometricStressPathUpdates const & ) = default;

  /// Default move constructor
  OedometricStressPathUpdates( OedometricStressPathUpdates && ) = default;

  /// Deleted default constructor
  OedometricStressPathUpdates() = default;

  /// Deleted copy assignment operator
  OedometricStressPathUpdates & operator=( OedometricStressPathUpdates const & ) = delete;

  /// Deleted move assignment operator
  OedometricStressPathUpdates & operator=( OedometricStressPathUpdates && ) =  delete;
  
  /**
   * @brief Evaluate the effective aperture, and its derivative wrt aperture
   * @param[in] aperture the model aperture/gap
   * @param[out] dHydraulicAperture_dAperture the derivative of the effective aperture wrt aperture
   * @return The hydraulic aperture that is always > 0
   */
  GEOS_HOST_DEVICE
  real64 computeHydraulicAperture( real64 const aperture,
                                   real64 const normalTraction,
                                   integer const fractureState,
                                   real64 & dHydraulicAperture_aperture,
                                   real64 & dHydraulicAperture_dNormalStress ) const;
private:
  real64 m_aperture0;
}
/**
 * @class Holds parameters and status for execution of nonlinear solution schemes.
 */
class OedometricStressPath : public PrescribedStressPathGeomechanicsBase
{
public:

  /**
   * @brief Constructor
   * @param[in] name The name of the new instantiation of this Group.
   * @param[in] parent A pointer to the parent of this Group.
   */
  OedometricStressPath( string const & name,
                        Group * const parent );

  /**
   * @brief The name of this object in the catalog.
   * @return A string containing the name of this object in the catalog.
   * The CatalogName is the string that will result in the creation of a new
   * OedometricStressPath object when calling
   * Group::getCatalog()::Allocate().
   */
  static string catalogName() { return "OedometricStressPath"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void postInputInitialization() override;

  void print() const;

  struct viewKeysStruct
  {
    static constexpr char const * biotString()                 { return "biot"; }
    static constexpr char const * poissonString()              { return "poisson"; }
    static constexpr char const * referencePressureString()    { return "referencePressure"; }
    static constexpr char const * referenceTotalStressString() { return "referenceTotalStress"; }
  } viewKeys;

  virtual real64 computeFractureStress( real64 const pressure, R1Tensor const & normal ) const override;
  
private:
  real64 m_biot;
  real64 m_poisson;
  real64 m_referencePressure; // p_0
  
  R1Tensor m_referenceTotalStress; // sigmaT_0 computed analytically
  R1Tensor m_referenceEffectiveStress; // sigma_0

};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_PERMEABILITY_OEDOMETRICSTRESSPATH_HPP_ */
