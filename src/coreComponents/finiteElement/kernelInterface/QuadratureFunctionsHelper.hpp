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
 * @file QuadratureFunctionsHelper.hpp
 */

#ifndef GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_
#define GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
// namespace finiteElement
// {

real64 determinant( real64 const (& J)[3][3] )
{
    real64 const detJ = J[0][0] * (J[1][1] * J[2][2] - J[2][1] * J[1][2])
                      - J[1][0] * (J[0][1] * J[2][2] - J[2][1] * J[0][2])
                      + J[2][0] * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
    return detJ;
}

void adjugate( real64 const (& J)[3][3], real64 (& AdjJ)[3][3] )
{
    AdjJ[0][0] = (J[1][1] * J[2][2]) - (J[1][2] * J[2][1]);
    AdjJ[0][1] = (J[2][1] * J[0][2]) - (J[0][1] * J[2][2]);
    AdjJ[0][2] = (J[0][1] * J[1][2]) - (J[1][1] * J[0][2]);
    AdjJ[1][0] = (J[2][0] * J[1][2]) - (J[1][0] * J[2][2]);
    AdjJ[1][1] = (J[0][0] * J[2][2]) - (J[0][2] * J[2][0]);
    AdjJ[1][2] = (J[1][0] * J[0][2]) - (J[0][0] * J[1][2]);
    AdjJ[2][0] = (J[1][0] * J[2][1]) - (J[2][0] * J[1][1]);
    AdjJ[2][1] = (J[2][0] * J[0][1]) - (J[0][0] * J[2][1]);
    AdjJ[2][2] = (J[0][0] * J[1][1]) - (J[0][1] * J[1][0]);
}

// } // namespace finiteElement
} // namespace geosx



#endif /* GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_ */
