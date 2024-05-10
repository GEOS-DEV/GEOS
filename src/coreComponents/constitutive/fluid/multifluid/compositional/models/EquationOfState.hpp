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
 * @file EquationOfState.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_EQUATIONOFSTATE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_EQUATIONOFSTATE_HPP_

#include "common/DataTypes.hpp"
#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

enum class EquationOfStateType : integer
{
  PengRobinson,
  SoaveRedlichKwong
};

class EquationOfState final
{
public:
  EquationOfState( string_array const & names );
  ~EquationOfState() = default;
  EquationOfState( const EquationOfState & ) = delete;
  const EquationOfState & operator=( const EquationOfState & ) = delete;

  /**
   * @brief A kernel wrapper for the equation of state parameters
   */
  struct KernelWrapper
  {
    KernelWrapper( arrayView1d< integer const > const & types );

    /**
     * @brief Move the KernelWrapper to the given execution space, optionally touching it.
     * @param space the space to move the KernelWrapper to
     * @param touch whether the KernelWrapper should be touched in the new space or not
     * @note This function exists to enable holding KernelWrapper objects in an ArrayView
     *       and have their contents properly moved between memory spaces.
     */
    void move( LvArray::MemorySpace const space, bool const touch )
    {
      m_types.move( space, touch );
    }

    arrayView1d< integer const > m_types;
  };

  /**
   * @brief Function to create and return a KernelWrapper
   * @return the KernelWrapper object
   */
  KernelWrapper createKernelWrapper() const;

  /**
   * @brief Validate list of equation of state names
   * @details This will raise an error if one of the names is invalid
   * @param[in] names The names of the equations of state
   * @return @c true if valid else @c false
   */
  static bool validateNames( string_array const & names );

  /**
   * @brief Check of an integer EOS type is the same as an enum value
   * @param[in] type The integer type
   * @param[in] eos The @c enum value
   * @return @c true if equal otherwise @c false
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static bool equals( integer const & type, EquationOfStateType const & eos )
  {
    return static_cast< integer >(eos) == type;
  }

  /**
   * @brief Convert equation of state from names to integral types
   * @param[in] names The names of the equations of state
   * @param[out] types The names converted to integer types. This will be resized as necessary
   */
  static void convertNames( string_array const & names, array1d< integer > & types );

private:
  string_array const & m_names;

  array1d< integer > m_types;
};

ENUM_STRINGS( EquationOfStateType,
              "pr",
              "srk" );

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_EQUATIONOFSTATE_HPP_
