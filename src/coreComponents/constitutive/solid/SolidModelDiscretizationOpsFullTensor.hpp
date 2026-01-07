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
 * @file SolidModelDiscretizationOpsFullTensor.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIAZTIONOPSFULLTENSOR_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIAZTIONOPSFULLTENSOR_HPP_

#include "SolidModelDiscretizationOps.hpp"

namespace geos
{
namespace constitutive
{

struct SolidModelDiscretizationOpsFullTensor : public SolidModelDiscretizationOps
{
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void BTDB( BASIS_GRADIENT const & gradN,
             real64 const & detJxW,
             real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  GEOS_HOST_DEVICE
  inline
  void scaleParams( real64 const scale )
  {
    LvArray::tensorOps::scale< 9, 9 >( m_c, scale );
  }

  real64 m_c[9][9];
};

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsFullTensor::BTDB( BASIS_GRADIENT const & gradN,
                                                  real64 const & detJxW,
                                                  real64 (& elementStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  real64 B[9][NUM_SUPPORT_POINTS * 3];
  
  for (int a = 0; a < NUM_SUPPORT_POINTS; ++a) {
    for (int i = 0; i < 3; ++i) {
      int col = 3 * a + i;
      for (int j = 0; j < 3; ++j) {
        int row = 3 * i + j;
        B[row][col] = gradN[a][j];
      }
    }
  }

  real64 DB[9][NUM_SUPPORT_POINTS * 3];
  LvArray::tensorOps::Rij_eq_AikBkj<9, NUM_SUPPORT_POINTS * 3, 9>(DB, m_c, B);
  LvArray::tensorOps::Rij_eq_AkiBkj<NUM_SUPPORT_POINTS * 3, NUM_SUPPORT_POINTS * 3, 9>(elementStiffness, B, DB);
  LvArray::tensorOps::scale<NUM_SUPPORT_POINTS * 3, NUM_SUPPORT_POINTS * 3>(elementStiffness, detJxW);
}

} // namespace constitutive
} // namespace geos

#endif /* GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIAZTIONOPSFULLTENSOR_HPP_ */