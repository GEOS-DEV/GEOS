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

#ifndef GEOSX_TOTOGAZ_STABLETOXML_HPP
#define GEOSX_TOTOGAZ_STABLETOXML_HPP

#include "dataRepository/xmlWrapper.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{
namespace api
{

void Convert( string const & stableInputFileName, xmlWrapper::xmlDocument & doc );

}
}

#endif //GEOSX_TOTOGAZ_STABLETOXML_HPP
