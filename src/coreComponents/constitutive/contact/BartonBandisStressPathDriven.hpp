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

/**
 * @file BartonBandis.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_BARTONBANDISSTRESSPATHDRIVEN_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_BARTONBANDISSTRESSPATHDRIVEN_HPP_

#include "constitutive/contact/HydraulicApertureBase.hpp"
#include "functions/TableFunction.hpp"
#include "physicsSolvers/solidMechanics/contact/FractureState.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class BartonBandisUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class BartonBandisStressPathDrivenUpdates
{
public:

  BartonBandisStressPathDrivenUpdates( real64 const aperture0,
                                       real64 const biot,
                                       real64 const poisson,
                                       real64 const normalStiffness,
                                       real64 const referencePressure,
                                       R1Tensor const &referenceTotalStress)
    : m_aperture0( aperture0 ),
      m_biot( biot ),
      m_poisson( poisson ),
      m_normalStiffness( normalStiffness ), // Kni
      m_referencePressure( referencePressure ),
      m_referenceTotalStress( referenceTotalStress )
  {
    m_referenceEffectiveStress[0] = m_referenceTotalStress[0] - m_biot*m_referencePressure; 
    m_referenceEffectiveStress[1] = m_referenceTotalStress[1] - m_biot*m_referencePressure; 
    m_referenceEffectiveStress[2] = m_referenceTotalStress[2] - m_biot*m_referencePressure; 
  }

  /// Default copy constructor
  BartonBandisStressPathDrivenUpdates( BartonBandisStressPathDrivenUpdates const & ) = default;

  /// Default move constructor
  BartonBandisStressPathDrivenUpdates( BartonBandisStressPathDrivenUpdates && ) = default;

  /// Deleted default constructor
  BartonBandisStressPathDrivenUpdates() = default;

  /// Deleted copy assignment operator
  BartonBandisStressPathDrivenUpdates & operator=( BartonBandisStressPathDrivenUpdates const & ) = delete;

  /// Deleted move assignment operator
  BartonBandisStressPathDrivenUpdates & operator=( BartonBandisStressPathDrivenUpdates && ) =  delete;

  /**
   * @brief Evaluate the effective aperture, and its derivative wrt aperture
   * @param[in] aperture the model aperture/gap
   * @param[out] dHydraulicAperture_dAperture the derivative of the effective aperture wrt aperture
   * @return The hydraulic aperture that is always > 0
   */
  GEOS_HOST_DEVICE
  real64 computeHydraulicAperture( real64 const pressure, 
                                   array1d< real64 > const & normal) const;


private:
  real64 m_aperture0;
  
  real64 m_biot;
  real64 m_poisson;
  real64 m_normalStiffness; // Kni
  real64 m_referencePressure; // p_0
  
  R1Tensor m_referenceTotalStress; // sigmaT_0 computed analytically
  R1Tensor m_referenceEffectiveStress; // sigma_0

  real64 computeFractureStress( real64 const pressure, array1d< real64 > const & normal ) const;
};


/**
 * @class BartonBandis
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
class BartonBandisStressPathDriven : public HydraulicApertureBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  BartonBandisStressPathDriven( string const & name,
                                Group * const parent );

  static string catalogName() { return "BartonBandisStressPathDriven"; }

  virtual string getCatalogName() const override { return catalogName(); }


  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BartonBandisStressPathDrivenUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  struct viewKeyStruct : public HydraulicApertureBase::viewKeyStruct
  {
    static constexpr char const * biotString()                  { return "biot"; }
    static constexpr char const * poissonString()               { return "poisson"; }
    static constexpr char const * normalStiffnessString()       { return "normalStiffness"; }
    static constexpr char const * referencePressureString()     { return "referencePressure"; }
    static constexpr char const * referenceTotalStressString()  { return "referenceTotalStress"; }
  };

private:
  real64 m_biot;
  real64 m_poisson;
  real64 m_normalStiffness; // Kni  
  real64 m_referencePressure; // p_0
  
  R1Tensor m_referenceTotalStress; // sigmaT_0 computed analytically
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 BartonBandisStressPathDrivenUpdates::computeHydraulicAperture( real64 const pressure, 
                                                                      array1d< real64 > const & normal ) const
{
  real64 const biot_pressure = m_biot * m_referencePressure; // biot is alpha in the equations

  // Computation of maximum fracture closure (Barton-Bandis parameter)
  // Fracture traction via Terzaghi's Principle
  real64 const sigma_c0[3] = { m_referenceTotalStress[0] * normal[0] - biot_pressure * normal[0],
                              m_referenceTotalStress[1] * normal[1] - biot_pressure * normal[1],
                              m_referenceTotalStress[2] * normal[2] - biot_pressure * normal[2] };
  real64 const sigma_n0 = sigma_c0[0]*normal[0] + 
                          sigma_c0[1]*normal[1] + 
                          sigma_c0[2]*normal[2];
  real64 const g0 = (-m_normalStiffness*m_aperture0 + 
                      std::sqrt((m_normalStiffness*m_aperture0)*
                      (m_normalStiffness*m_aperture0) + 
                      4.0*m_normalStiffness*sigma_n0*m_aperture0)) / (2.0*m_normalStiffness);
  real64 const maximumFractureClosure = g0 + m_aperture0; // Vm

  // Normal effective stress on the fracture
  real64 const sigmaN_N = computeFractureStress( pressure, normal );
  real64 const fractureClosure = sigmaN_N*maximumFractureClosure/(m_normalStiffness*maximumFractureClosure + sigmaN_N); // gn_BB

  // Compute the new aperture which is equal to the aperture at the free-stress state 
  // minus the closure from the free-stress state to the current state
  real64 const newHydraulicAperture = maximumFractureClosure - fractureClosure;

  return newHydraulicAperture;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 BartonBandisStressPathDrivenUpdates::computeFractureStress( real64 const pressure,
                                                                   array1d< real64 > const & normal ) const
{  
  // TODO: remove this lambda expression
  auto matmul = [](real64 const (&u)[3], array1d< real64 > const &v, array1d< real64 > &r) -> void
  {
    r[0] = u[0]*v[0];
    r[1] = u[1]*v[1];
    r[2] = u[2]*v[2];
  }; 
  
  real64 const deltaSigmaZ = m_biot * (pressure - m_referencePressure);
  real64 const poisson_deltaSigma = deltaSigmaZ * m_poisson/(1.0 - m_poisson);
  // sigma: matrix diagonal
  real64 effectiveStress[3] = { m_referenceEffectiveStress[0] - poisson_deltaSigma,
                                m_referenceEffectiveStress[1] - poisson_deltaSigma,
                                m_referenceEffectiveStress[2] - deltaSigmaZ };
  array1d< real64 > effectiveStressOnFracture(3); // sigma_c
  matmul(effectiveStress, normal, effectiveStressOnFracture);
  real64 normalComponentOfStressOnFracture = dot(effectiveStressOnFracture, normal); // sigmaN_N
  return normalComponentOfStressOnFracture;
}
} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_BARTONBANDISSTRESSPATHDRIVEN_HPP_ */
