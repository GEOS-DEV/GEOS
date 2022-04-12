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

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_S
#define GEOSX_LINEARALGEBRA_UTILITIES_ARNOLDI_HPP_

namespace geosx
{

template< typename T >
class AsyncRequest 
{
public:

template< typename F >
AsyncRequest( F f )
{
  f(  )
}

T const & complete() const
{
  MpiWrapper::wait( &m_request, MPI_STATUS_IGNORE );
  return m_result;
}

private:

  MPI_Request m_request{};
  T m_result{};

};

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_UTILITIES_ARNOLDI_HPP_
