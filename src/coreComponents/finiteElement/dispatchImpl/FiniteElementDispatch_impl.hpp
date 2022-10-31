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
 * @file FiniteElementDispatch.hpp
 */

#include "finiteElement/FiniteElementDispatch.hpp"

namespace geosx
{
namespace finiteElement
{

//template< typename FE_TYPE >
//void
//dispatch3DImpl( FE_TYPE const & input, std::function< void( FE_TYPE const & ) > & lambda )
//{
//  lambda( input );
//}

template< typename FE_TYPE >
void
dispatch3DImpl( FE_TYPE & input, std::function< void( FE_TYPE & ) > & lambda )
{
  lambda( input );
}

template< typename FE_TYPE >
void
dispatch2DImpl( FE_TYPE const & input, std::function< void( FE_TYPE const & ) > & lambda )
{
  lambda( input );
}

}
}

