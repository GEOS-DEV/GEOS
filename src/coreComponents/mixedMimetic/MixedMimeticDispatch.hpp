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
 * @file MixedMimeticDispatch.hpp
 */

#ifndef GEOS_MIXEDMIMETIC_MIXEDMIMETICDISPATCH_HPP_
#define GEOS_MIXEDMIMETIC_MIXEDMIMETICDISPATCH_HPP_

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/SimpleInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"
#include "common/logger/Logger.hpp"

namespace geos
{

/**
 * @brief Dispatch the mixed-form inner products, i.e. the subset of mimetic inner products
 *        implementing the mass matrix construction (computeM) used by the mixed mimetic solvers.
 * @tparam LAMBDA type of the user-provided generic lambda
 * @param[in] input the mimetic inner product base reference to dispatch
 * @param[in] lambda the user-provided generic lambda
 */
template< typename LAMBDA >
void
mixedMimeticInnerProductDispatch( mimeticInnerProduct::MimeticInnerProductBase const & input,
                                  LAMBDA && lambda )
{
  if( auto const * const ptr1 = dynamic_cast< mimeticInnerProduct::TPFAInnerProduct const * >(&input) )
  {
    lambda( *ptr1 );
  }
  else if( auto const * const ptr2 = dynamic_cast< mimeticInnerProduct::QuasiTPFAInnerProduct const * >(&input) )
  {
    lambda( *ptr2 );
  }
  else if( auto const * const ptr3 = dynamic_cast< mimeticInnerProduct::SimpleInnerProduct const * >(&input) )
  {
    lambda( *ptr3 );
  }
  else if( auto const * const ptr4 = dynamic_cast< mimeticInnerProduct::BdVLMInnerProduct const * >(&input) )
  {
    lambda( *ptr4 );
  }
  else
  {
    GEOS_ERROR( "mixedMimeticInnerProductDispatch() is not implemented for input of "<<LvArray::system::demangleType( input ) );
  }
}

} // namespace geos

#endif //GEOS_MIXEDMIMETIC_MIXEDMIMETICDISPATCH_HPP_
