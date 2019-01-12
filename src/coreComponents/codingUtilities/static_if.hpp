/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * static_if.hpp
 *
 *  Created on: Feb 6, 2018
 *      Author: settgast1
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_STATIC_IF_HPP_
#define SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_STATIC_IF_HPP_


namespace geosx
{

template<bool CONDITION>
struct static_if_wrapper
{
  template<typename LAMBDA_BODY>
  constexpr inline static void if_function(LAMBDA_BODY&&) { }
};

template<>
struct static_if_wrapper<true>
{
  template<typename LAMBDA_BODY>
  constexpr inline static void if_function(LAMBDA_BODY&& lambda)
  {
    lambda();
  }
};

}

#define static_if( condition ) \
  geosx::static_if_wrapper<condition>::if_function( [&] () -> void

#define end_static_if );

#endif /* SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_STATIC_IF_HPP_ */
