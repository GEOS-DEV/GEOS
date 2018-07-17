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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */




#ifndef TENSOR_OPS_H_
#define TENSOR_OPS_H_

#include <cmath>

typedef double realT;


inline int e_ijk( const int i, const int j, const int k )
{
  int r_value = -10000;

  // check to see if any of the indices are the same
  if( ( i==j ) || ( i==k ) || ( j==k ) )
    r_value = 0;
  else if( i==1 )
  {
    if( j==2 )
      r_value = 1;
    else
      r_value = -1;
  }
  else if( j==1 )
  {
    if( i==2 )
      r_value = -1;
    else
      r_value = 1;
  }
  else if( k==1 )
  {
    if( i==2 )
      r_value = 1;
    else
      r_value = -1;
  }
  return r_value;
}

inline int d_ij( const int i, const int j )
{
  if( i==j )
    return 1;
  else
    return 0;
}

inline int is_equal( const realT& a, const realT& b, const int num_digits )
{
  int rval=0;
  const realT mean = (a + b)*0.5;
  const realT diff = fabs(a - b);
  if( diff < mean*pow(10.0,-num_digits) )
    rval = 1;

  return rval;
}



#endif /* _CONSTANTS_H_ */
