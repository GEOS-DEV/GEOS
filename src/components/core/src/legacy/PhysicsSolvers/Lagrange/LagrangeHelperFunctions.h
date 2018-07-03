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
 * LagrangeHelperFunctions.h
 *
 *  Created on: Nov 7, 2012
 *      Author: settgast1
 */

#ifndef LAGRANGEHELPERFUNCTIONS_H_
#define LAGRANGEHELPERFUNCTIONS_H_


#include "legacy/Common/Common.h"
//#include "legacy/ObjectManagers/DiscreteElementManagerT.h"

namespace LagrangeHelperFunctions
{

void LinearPointUpdatePart1( ObjectDataStructureBaseT& objectManager,
                             const realT& time,
                             const realT& dt,
                             const bool clearForces = true);

void LinearPointUpdatePart2( ObjectDataStructureBaseT& objectManager,
                             const realT& time,
                             const realT& dt,
                             const R1Tensor& gravity,
                             const realT& damping = 0.0,
                             const realT& stiffDamping = 0.0 );

//  void RotationalPointUpdatePart1(DiscreteElementManagerBaseT&
// discreteElementManager,
//                                  const realT& time, const realT& dt);
//
//  void RotationalPointUpdatePart1b(DiscreteElementManagerT&
// discreteElementManager);


void RotationalPointUpdatePart2(ObjectDataStructureBaseT& objectManager,
                                const realT& time, const realT& dt);

realT CalculateMaxStableExplicitTimestep( const realT& density,
                                          const realT& stiffness,
                                          const realT& BB );

inline realT BulkQ( const realT& rho,
                    const realT& soundSpeed,
                    const realT& qLinear,
                    const realT& qQuadratic,
                    const realT& trD,
                    const realT& L )
{
  realT Q = 0;
  if( trD < 0.0 )
  {
    Q = rho * ( soundSpeed * qLinear + qQuadratic  * ( trD * L ) )  * ( trD * L );
  }
  return Q;
}

}



#endif /* LAGRANGEHELPERFUNCTIONS_H_ */
