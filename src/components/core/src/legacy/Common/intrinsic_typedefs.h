/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file intrinsic_typedefs.h
 * @author settgast1
 * @date Feb 10, 2011
 */

#ifndef INTRINSIC_TYPEDEFS_H_
#define INTRINSIC_TYPEDEFS_H_

#include <string>
#include <vector>
#include <climits>
#include "common/DataTypes.hpp"

/// set the default iteration type
typedef unsigned int INDEX;

/// this is the typedef to define what a "realT" data is
typedef double realT;

/// this is the typedef to define what type is used for local indexing and local
// ID's
//typedef std::vector<int>::size_type localIndex;

#define LOCALINDEX_MIN 0
#define LOCALINDEX_MAX ULLONG_MAX


/// this is the typedef to define what type is used for global indexing and
// global ID's
//typedef unsigned long long globalIndex;
#define GLOBALINDEX_MIN 0
#define GLOBALINDEX_MAX ULLONG_MAX

/// a vector of ints
typedef std::vector<int> ivector;

/// a vector of reals
typedef std::vector<realT> dvector;


//typedef std::vector<localIndex> lvector;


//typedef std::vector<globalIndex> gvector;



#endif /* INTRINSIC_TYPEDEFS_H_ */
