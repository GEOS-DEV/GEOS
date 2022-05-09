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

#ifndef GEOSX_API_STABLETOXML_HPP
#define GEOSX_API_STABLETOXML_HPP

#include "dataRepository/xmlWrapper.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{
namespace api
{

/**
 * @brief Converts the stable input file to its xml internal representation
 * @param stableInputFileName The yaml input file name.
 * @param doc The output xml document being fed.
 */
void Convert( string const & stableInputFileName, xmlWrapper::xmlDocument & doc );

}
}

#endif //GEOSX_API_STABLETOXML_HPP
