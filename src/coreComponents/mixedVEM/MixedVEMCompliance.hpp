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
 * @file MixedVEMCompliance.hpp
 */

#ifndef GEOS_MIXEDVEM_MIXEDVEMCOMPLIANCE_HPP_
#define GEOS_MIXEDVEM_MIXEDVEMCOMPLIANCE_HPP_

#include "mixedVEM/MixedVEMTypes.hpp"

namespace geos
{

namespace mixedVEM
{

/**
 * @brief Compliance D = C^{-1} of a general anisotropic elasticity tensor.
 * @param[in] stiffness the 6x6 matrix of C on the orthonormal symmetric basis
 * @param[out] compliance the 6x6 matrix of D on the same basis
 *
 * The bilinear form a_E uses D, whereas constitutive models provide C, so the mixed
 * formulation needs this inversion once per material. Host only.
 */
void invertStiffness( real64 const (&stiffness)[NUM_SYM_COMP][NUM_SYM_COMP],
                      real64 ( &compliance )[NUM_SYM_COMP][NUM_SYM_COMP] );

/**
 * @brief Convert a Voigt stiffness with engineering shear strains to the orthonormal basis.
 * @param[in] voigtStiffness the 6x6 Voigt matrix, ordered (xx, yy, zz, yz, xz, xy)
 * @param[out] stiffness the 6x6 matrix on the orthonormal symmetric basis
 *
 * Engineering shear doubles the strain entries, so the shear rows and columns pick up
 * factors of sqrt(2) and 2 when moved to the basis for which pi_a : pi_b = delta_ab.
 */
void convertVoigtStiffness( real64 const (&voigtStiffness)[NUM_SYM_COMP][NUM_SYM_COMP],
                            real64 ( &stiffness )[NUM_SYM_COMP][NUM_SYM_COMP] );

} // namespace mixedVEM

} // namespace geos

#endif // GEOS_MIXEDVEM_MIXEDVEMCOMPLIANCE_HPP_
