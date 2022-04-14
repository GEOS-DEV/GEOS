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

/**
 * @brief Helper class for managing non-blocking operation.
 * @tparam T the type of the result of the non-blocking operation
 */
template< typename T >
class AsyncRequest
{
public:

  /**
   * @brief Constructor.
   * @tparam F the type of the lambda function
   * @param f lambda function that defines the MPI non-blocking operation to be performed
   */
  template< typename F >
  explicit AsyncRequest( F f )
    : m_function( std::move( f ) )
  {
    m_result = std::make_unique< T >( T{} );
    m_request = std::make_unique< MPI_Request >();
    m_function( *m_request, *m_result );
  }

  /**
   * @brief Completes the managed MPI request.
   * @return the result of the non blocking operation
   */
  T const & complete() const
  {
    MpiWrapper::wait( m_request.get(), MPI_STATUS_IGNORE );
    return *m_result;
  }

private:

  std::unique_ptr< MPI_Request > m_request;              ///< Unique pointer to the MPI request
  std::unique_ptr< T > m_result;                         ///< Unique pointer to the result of the non-blocking operation
  std::function< void( MPI_Request &, T & ) > m_function; ///< Lambda function wrapper to a non-blocking MPI function

};

} // namespace geosx

#endif // GEOSX_LINEARALGEBRA_UTILITIES_ASYNCREQUEST_HPP_
