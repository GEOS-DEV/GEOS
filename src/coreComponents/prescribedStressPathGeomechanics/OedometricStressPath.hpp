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

#ifndef GEOS_PHYSICSSOLVERS_OEDOMETRICSTRESSPATH_HPP_
#define GEOS_PHYSICSSOLVERS_OEDOMETRICSTRESSPATH_HPP_

#include "dataRepository/Group.hpp"

namespace geos
{

/**
 * @class Holds parameters and status for execution of nonlinear solution schemes.
 */
class OedometricStressPath : public dataRepository::Group
{
public:

  /**
   * @brief Default constructor.
   */
  OedometricStressPath() = delete;

  /**
   * @brief Constructor
   * @param[in] name The name of the new instantiation of this Group.
   * @param[in] parent A pointer to the parent of this Group.
   */
  OedometricStressPath( string const & name,
                        Group * const parent );

  /**
   * @brief Default destructor
   */
  virtual ~OedometricStressPath() = default;

  /**
   * @brief Default Move Constructor
   * @param The source object of the move.
   */
  OedometricStressPath( OedometricStressPath && ) = default;

  /**
   * @brief Copy Constructor
   * @param The source object.
   */
  OedometricStressPath & operator=( const OedometricStressPath & params )
  {
    m_biot = params.m_biot;
    m_poisson = params.m_poisson;
    m_referencePressure = params.m_referencePressure;
    m_referenceTotalStress = params.m_referenceTotalStress;
    m_referenceEffectiveStress = params.m_referenceEffectiveStress;

    setLogLevel( params.getLogLevel());

    return *this;
  }

  /**
   * @brief The name of this object in the catalog.
   * @return A string containing the name of this object in the catalog.
   * The CatalogName is the string that will result in the creation of a new
   * NonlinearSolverParameters2 object when calling
   * Group::getCatalog()::Allocate().
   */
  static string catalogName() { return "OedometricStressPath"; }

  virtual void postInputInitialization() override;

  void print() const;

  struct viewKeysStruct
  {
    static constexpr char const * biotString()                 { return "biot"; }
    static constexpr char const * poissonString()              { return "poisson"; }
    static constexpr char const * referencePressureString()    { return "referencePressure"; }
    static constexpr char const * referenceTotalStressString() { return "referenceTotalStress"; }
  } viewKeys;

  real64 computeFractureStress( real64 const pressure, R1Tensor const & normal ) const;
  
private:
  real64 m_biot;
  real64 m_poisson;
  real64 m_referencePressure; // p_0
  
  R1Tensor m_referenceTotalStress; // sigmaT_0 computed analytically
  R1Tensor m_referenceEffectiveStress; // sigma_0

};


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_OEDOMETRICSTRESSPATH_HPP_ */
