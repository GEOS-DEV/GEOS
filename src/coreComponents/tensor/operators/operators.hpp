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
 * @file operators.hpp
 */

#ifndef GEOSX_TENSOR_OPERATORS
#define GEOSX_TENSOR_OPERATORS

/**
 * These are the main mathematical operations applicable on tensors.
 * */

/// Get a tensor slice (lazy operation)
#include "get.hpp"
/// Determinant operators for matrices
#include "determinant.hpp"
/// Point-wise multiplications at quadrature points
#include "point-wise_multiplications/point-wise_multiplications.hpp"
/// Dot product
#include "dot_product.hpp"
/// Norm functions
#include "norms.hpp"
/// Matrix multiplications
#include "matrix_multiplication.hpp"
/// An identity operator (useful when no preconditioner for instance)
#include "identity.hpp"
/// A generic conjugate gradient algorithm
#include "conjugate_gradient.hpp"


/**
 * These are the main mathematical operations using a Basis, Degrees of
 * freedom, and QData.
 * */

/// Basis contractions for the tensors
#include "contractions/contractions.hpp"
/// Interpolation operators at quadrature point, ex: B * u
#include "interpolations/interpolations.hpp"
/// Gradients operators at quadratire point, ex: grad(B) * u
#include "gradients/gradients.hpp"
/// Divergence operators at quadrature point, ex: div(B) * u
#include "divergence/divergence.hpp"
/// Curl operators at quadrature point, ex: curl(B) * u
#include "curl/curl.hpp"
/// Quadrature operators
#include "quadrature_operator.hpp"
/// Element local operators
#include "element_operator.hpp"

#endif // GEOSX_TENSOR_OPERATORS
