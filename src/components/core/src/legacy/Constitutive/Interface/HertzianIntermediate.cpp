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
 * HertzianIntermediate.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */


#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "HertzianIntermediate.h"
#include <typeinfo>
#include <assert.h>
#include "../../IO/ticpp/HierarchicalDataNode.h.old"

HertzianIntermediate::HertzianIntermediate( const int paramSize, const int stateSize ):
  InterfaceBase( paramSize, stateSize )
{
  // TODO Auto-generated constructor stub

}

HertzianIntermediate::~HertzianIntermediate()
{}

void
HertzianIntermediateStateData::Update(const realT curvature1,
                                      const realT curvature2)
{
  radius = EffectiveRadius(curvature1, curvature2);
  //set other derived values:
  //get the Hertzian coefficient, k, where f = k*a^1.5, k = Hertzian coefficient
  hertzCf = youngs * sqrt(radius) * 4.0 / 3.0;
}

void
HertzianIntermediateStateData::Initialize(const realT curvature1, const realT curvature2,
                                          const realT poissons1, const realT poissons2,
                                          const realT youngs1, const realT youngs2,
                                          const realT mass1, const realT mass2,
                                          const realT, const realT,
                                          const realT, const realT,
                                          const realT, const realT,
                                          const realT, const realT,
                                          const realT, const realT )
{
  youngs = EffectiveYoungsModulus(poissons1, poissons2, youngs1, youngs2);
  mass = EffectiveMass(mass1, mass2);

  //set mechanical properties of the contact
  Update(curvature1, curvature2);
}

realT
HertzianIntermediateStateData::EffectiveRadius(const realT curvature1,
                                               const realT curvature2) const
{
  realT r = curvature1 + curvature2;
  r = 1. / r;
  return r;
}

realT
HertzianIntermediateStateData::EffectiveYoungsModulus(const realT poissons1,
                                                      const realT poissons2,
                                                      const realT youngs1,
                                                      const realT youngs2) const
{
  realT Eeff = (1 - poissons1 * poissons1) / youngs1;
  Eeff += (1 - poissons2 * poissons2) / youngs2;
  Eeff = 1. / Eeff;
  return Eeff;
}

realT
HertzianIntermediateStateData::EffectiveMass(const realT mass1,
                                             const realT mass2) const
{
  realT em1 = 1/mass1 + 1/mass2;
  em1 = 1. / em1;
  return em1;
// EBH comment: if either mass1 or mass2 is large due to boundary conditions,
// etc. (e.g. 1.e100)
// EBH comment: then c++ will produce garbage for the line below.  The fix above
// accounts for this.
//  return mass1 * mass2 / (mass1 + mass2);
}
void
HertzianIntermediate::UpdateProperties(const localIndex index, std::map<std::string, realT>& p0, std::map<std::string, realT>& p1)
{
  if(index >= NumStateIndex0())
    throw GPException("index of out range: HertzianIntermediate::UpdateProperties\n");

  //radius(0), poissons(1), youngs(2), mass(3), rest(4),
  //yield(5), velHalf(6), surfaceEnergy(7), cement(8)

  HertzianIntermediateStateData& hstate = *StateData(index, 0);
  if (hstate.youngs > 0)
    hstate.Update     ( p0["radius"], p1["radius"]);
  else
    hstate.Initialize ( p0["radius"], p1["radius"],
                        p0["Nu"],     p1["Nu"],
                        p0["E"],      p1["E"],
                        p0["mass"],   p1["mass"],
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}
realT
HertzianIntermediate::NormalStiffness(const InterfaceBaseParameterData&,
                                      InterfaceBaseStateData& matStateBase,
                                      const realT normalApproach,
                                      const bool setForces) const
{
  HertzianIntermediateStateData& matState = static_cast<HertzianIntermediateStateData&> (matStateBase);
  const realT stiffness = matState.hertzCf * sqrt(normalApproach);
  if(setForces)
    matState.stress = stiffness * normalApproach;
  return 1.5 * stiffness;
}

void
HertzianIntermediate::UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                                      InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  const HertzianIntermediateParameterData& matParams = static_cast < const HertzianIntermediateParameterData& > ( matParamsBase );
  HertzianIntermediateStateData& matState = static_cast < HertzianIntermediateStateData& > ( matStateBase );

  matState.mu = matParams.mu0;
}

realT
HertzianIntermediate::ShearStrength(const InterfaceBaseParameterData&,
                                    InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  HertzianIntermediateStateData& matState = static_cast < HertzianIntermediateStateData& > ( matStateBase );
  return matState.stress * matState.mu;
}
