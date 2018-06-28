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

/*
 * MaterialUtilities.h
 *
 *  Created on: Nov 27, 2013
 *      Author: walsh24
 */

#ifndef MATERIAL_UTILITIES_H_
#define MATERIAL_UTILITIES_H_

#include "Common/typedefs.h"
#include <sys/resource.h>
#include <map>

#include "../Common/GPException.h"

/*
 * FillLinearElasticModuli
 *
 *
 * @param K bulk Modulus
 * @param G Shear Modulus
 * @param E Youngs Modulus
 * @param nu Poissions ratio
 * @param lame Lames constant
 * @param M P-wave modulus
 * @return error code 1 or 2 if failed
 */
inline
int FillLinearElasticModuli(realT& K, realT& G, realT& E, realT& nu,realT& lame, realT& M, const bool forceRecalculate = false){

  if(K <= 0)
  {
    if(G>0 && E >0)
    {
      K = E*G/(3.0*(3.0*G-E));
    }
    else if (G>0  &&  nu >0)
    {
      K = 2.0*G*(1.0+nu)/(3.0*(1.0-2.0*nu));
    }
    else if (G>0  &&  lame >0)
    {
      K = lame + 2.0*G/3.0;
    }
    else if (G>0  &&  M >0)
    {
      K = M- 4.0*G/3.0;
    }
    else if(E>0  &&  nu >0)
    {
      K = E/(3.0*(1.0-2.0*nu));
    }
    else if (lame>0  &&  nu >0)
    {
      K = lame*(1.0+nu)/(3.0*nu);
    }
    else
    {
      return 1;
    }
  }


  if(G<=0)
  {
    if(E>0)
    {
      G = 3.0*K*E/(9.0*K-E);
    }
    else if(lame > 0)
    {
      G = 3.0*(K-lame)/2.0;
    }
    else if(nu > 0)
    {
      G = 3.0*K*(1.0-2*nu)/(2.0*(1+nu));
    }
    else if(M > 0)
    {
      G = 0.75*(M-K);
    }
    else
    {
      return 2;
    }
  }


  //K & G are set - fill in the others (but only if not already set to avoid
  // roundoff)
  if(forceRecalculate)
    E = nu = lame = M = 0;

  if(E<=0)
    E = 9*K*G/(3.0*K+G);
  if(nu<=0)
    nu = (3*K-2*G)/(2.0*(3.0*K+G));
  if(lame<=0)
    lame = K - 2.0*G/3.0;
  if(M<=0)
    M = K + 4.0*G/3.0;

  return 0;
}



#endif /* UTILITIES_H_ */
