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
 * @file AsyncRequest.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_ASYNCREQUEST_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_ASYNCREQUEST_HPP_

#include "common/MpiWrapper.hpp"

namespace geosx
{

template< typename T >
class AsyncRequest 
{
public:

template< typename F >
explicit AsyncRequest( F f )
  : m_function( std::move( f ) )
{
  m_result = std::make_unique< T >(T{});
  m_request = std::make_unique< MPI_Request >();
  m_function( *m_request, *m_result );
}

T const & complete() const
{
  MpiWrapper::wait( m_request.get(), MPI_STATUS_IGNORE );
  return *m_result;
}

private:

  std::unique_ptr< MPI_Request > m_request;
  std::unique_ptr< T > m_result;
  std::function< void(MPI_Request &, T & ) > m_function;

};

} // namespace geosx

#endif // GEOSX_LINEARALGEBRA_UTILITIES_ASYNCREQUEST_HPP_
