/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PVTFunctionHelpers.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONHELPERS_HPP
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONHELPERS_HPP

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

struct PVTFunctionHelpers
{

  /**
   * @brief Look for the expectedNames in the names provided by the user
   * @param[in] inputNames phase or component names provided by the user
   * @param[in] expectedNames expected names that can be accepted by GEOSX
   * @param[inout] id the index of the phase or component that was found
   * @return true is an expectedName was found in the inputNames, false otherwise
   */
  static bool
  findName( array1d< string > const & inputNames,
            array1d< string > const & expectedNames,
            localIndex & id )
  {
    id = -1;
    for( localIndex i = 0; i < inputNames.size(); ++i )
    {
      string const input = inputNames[i];
      for( localIndex j = 0; j < expectedNames.size(); ++j )
      {
        if( input == expectedNames[j] )
        {
          id = i;
          return true;
        }
      }
    }
    return false;
  }

};

} // namespace PVTProps

} // namespace constitutive

} // namespace geosx


#endif
