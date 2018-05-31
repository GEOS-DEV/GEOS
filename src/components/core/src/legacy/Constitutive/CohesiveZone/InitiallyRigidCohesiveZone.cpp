// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.


/*
 * InitiallyRigidCohesiveZone.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */


#include "Utilities/GeometryUtilities.h"
#include "DataStructures/VectorFields/ElementRegionT.h"
#include "Constitutive/Material/MaterialBase.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "InitiallyRigidCohesiveZone.h"
#include "Constitutive/CohesiveZone/CohesiveZoneFactory.h"
#include <typeinfo>
#include <assert.h>
#include "../../IO/ticpp/HierarchicalDataNode.h.old"
#include "../../IO/ticpp/HierarchicalDataNode.h.old"

InitiallyRigidCohesiveZone::InitiallyRigidCohesiveZone( ):
  CohesiveZoneBase( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}

InitiallyRigidCohesiveZone::~InitiallyRigidCohesiveZone()
{}

int
InitiallyRigidCohesiveZone::UpdateCohesiveZone( const localIndex index,
                                                const R1Tensor& gap,
                                                const R1Tensor& N,
                                                const std::pair< ElementRegionT*, localIndex >& elem0,
                                                const std::pair< ElementRegionT*, localIndex >& elem1,
                                                R1Tensor& traction,
                                                R2Tensor& stiffness )
{
  const InitiallyRigidCohesiveZoneParameterData& params = *ParameterData(index);
  InitiallyRigidCohesiveZoneStateData& state      = *StateData(index,0);

  int rval = 2;


  realT gapMag = gap.L2_Norm();

  realT k;

  const realT s = gapMag / params.failGap < 1.0 ?  gapMag / params.failGap : 1.0;
  bool mainline = false;

  state.maxGap.SetMax( gap );

  if( state.separationCoeff<s )
  {
    state.separationCoeff = s;
    mainline = true;
  }


  R1Tensor direction;

  traction = gap;
  traction.Normalize();

  if( isZero(gapMag,1e-14) )
  {
    R2SymTensor stress0;
    R2SymTensor stress1;


    ElementRegionT& er0 = *(elem0.first);
    ElementRegionT& er1 = *(elem1.first);
    const localIndex elemIndex0 = elem0.second;
    const localIndex elemIndex1 = elem1.second;

    realT pressure;
    er0.m_mat->MeanPressureDevStress(elemIndex0,pressure, stress0);
    stress0.PlusIdentity(pressure);

    er1.m_mat->MeanPressureDevStress(elemIndex1,pressure, stress1);
    stress1.PlusIdentity(pressure);

    stress0 += stress1;
    stress0 *= 0.5;

    direction.AijBj( stress0, N );
    direction.Normalize();
    gapMag = 1.0e-14;
  }
  else
  {
    direction = gap;
    direction /= gapMag;
  }

  direction[0] = 0;
  direction[2] = 0;
  //todo: why is this constrained this way??
  traction = direction;

  realT scale = 0.0;
  if( state.separationCoeff >=1.0 )
  {
    // fully separated
    state.maxTraction = 0.0;
    k = 0.0;
    traction = 0.0;
    rval = 3;
  }
  else if( isZero(state.separationCoeff) )
  {
    // initial
    traction *= params.failStress;
    k = -params.failStress/params.failGap;
  }
  else if( mainline || params.unloadFlag==0 )
  {
    // going down the curve
    scale = 1.0 - state.separationCoeff;
    state.maxTraction = params.failStress * scale;
    traction *= state.maxTraction;
    k = -params.failStress/params.failGap;
  }
  else
  {
    // unload/reload line
    scale = s/state.separationCoeff;
    traction *= scale * state.maxTraction;

    k = state.maxTraction/(state.separationCoeff*params.failGap);

  }


  realT dterm = traction.L2_Norm() / gapMag;

  R2Tensor kI;
  kI.PlusIdentity( dterm );

  stiffness.dyadic_ab( direction, direction );
  stiffness *= k - dterm;

  stiffness += kI;

  state.traction = traction;
  state.stiffness = stiffness;

  return rval;

}

/// Register class in the class factory
REGISTER_COHESIVEZONE( InitiallyRigidCohesiveZone )
