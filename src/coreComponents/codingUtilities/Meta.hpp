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

/**
 * @file Meta.hpp
 */

#ifndef GEOSX_META_HPP
#define GEOSX_META_HPP

namespace geosx
{

template<typename T, bool COND>
struct add_const_if
{
  using type = typename std::conditional<COND, typename std::add_const<T>::type, T>::type;
};

template<typename T, bool COND>
using add_const_if_t = typename add_const_if<T, COND>::type;

} // namespace geosx

#endif //GEOSX_META_HPP
