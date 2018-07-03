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
 * PenaltyCoulombIntermediate.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */


#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "PenaltyCoulombIntermediate.h"
#include <typeinfo>
#include <assert.h>
#include "../../IO/ticpp/HierarchicalDataNode.h.old"

PenaltyCoulombIntermediate::PenaltyCoulombIntermediate( const int paramSize, const int stateSize ):
  InterfaceBase( paramSize, stateSize )
{
  // TODO Auto-generated constructor stub

}

PenaltyCoulombIntermediate::~PenaltyCoulombIntermediate()
{}

realT
PenaltyCoulombIntermediateParameterData::Stiffness(const realT normalApproach) const
{
  realT arealStiffness = 0.0;

  if(normalApproach <= 0)
    arealStiffness = 0.0;
  else if(normalApproach < normalApproachYield)
    arealStiffness = ktildeAperture / (aperture - normalApproach);
  else if(normalApproach < normalApproachSoften)
    arealStiffness = kyield;
  else
    arealStiffness = ksoften;

  return arealStiffness;
}

realT
PenaltyCoulombIntermediateParameterData::Stress(const realT normalApproach) const
{
  realT normalStress = 0.0;

  if(normalApproach <= 0)
    normalStress = 0;
  else if(normalApproach < normalApproachYield)
    normalStress = ktildeAperture * log(aperture / (aperture - normalApproach));
  else if(normalApproach < normalApproachSoften)
    normalStress = stressYield + kyield * (normalApproach - normalApproachYield);
  else
    normalStress = stressSoften + ksoften * (normalApproach - normalApproachSoften);

  return normalStress;
}
void
PenaltyCoulombIntermediateParameterData::Initialize()
{
  //check user settings
  const realT tol = 1e-6;
  if(normalApproachYield > aperture)
    normalApproachYield = (1.0 - tol) * aperture;
  if(stressYield > stressSoften)
    stressYield = stressSoften;

  //set derived values
  const realT factor = 1.0/(aperture - normalApproachYield);
  ktildeAperture = stressYield / log(aperture * factor);
  kyield = ktildeAperture * factor;
  if(isEqual(stressSoften, std::numeric_limits<realT>::max()))
    normalApproachSoften = std::numeric_limits<realT>::max();
  else
    normalApproachSoften = normalApproachYield + (stressSoften - stressYield) / kyield;
}

void
PenaltyCoulombIntermediateParameterData::PostReadXML( TICPP::HierarchicalDataNode& )
{
  if(aperture <= 0)
    aperture = 1e-2;

  //normal approach at the onset of yielding
  if(normalApproachYield <= 0)
    throw GPException("You must set a normal approach for the onset of yielding > 0");
  if(normalApproachYield >= aperture)
    throw GPException("You must set a normal approach for the onset of yielding < aperture");
  if(stressYield <= 0)
    throw GPException("You must set a stress for the onset of yielding > 0");
  if(stressSoften <= 0)
    stressSoften = std::numeric_limits<realT>::max();
  if(stressSoften < stressYield)
    stressSoften = stressYield;

  //set all derived values
  Initialize();

  //set the rest ... arealStiffnessSoften is not calculated or used in
  // SimpleInitialize
  if(ksoften <= 0)
    ksoften = (1e-5) * kyield;
  if(kshear <= 0)
    kshear = 0.7 * kyield;
}

void
PenaltyCoulombIntermediateParameterData::PostSetValue()
{
  TICPP::HierarchicalDataNode node;
  PostReadXML(node);
}
void
PenaltyCoulombIntermediate::StrainDrivenUpdate( const localIndex index )
{
  //get temporary references
  const PenaltyCoulombIntermediateParameterData& matParams = *this->ParameterData(index);
  PenaltyCoulombIntermediateStateData& matState = *this->StateData(index, 0);

  //(1) evolve normal stress
  matState.stress = matParams.Stress(matState.normalApproach);

  //(2) evolve friction
  UpdateFriction(matParams, matState);

  //(3) evolve shear stress
  const realT dssdt = matParams.kshear * matState.dxsdt;
  matState.stressShear += dssdt * matState.dt;

  //(4) return shear stress to the failure surface if necessary
  ThresholdToFailureSurface(matParams, matState, matState.dxsdt * matState.dt);

  matState.stressShearVector.Normalize();
  matState.stressShearVector *= matState.stressShear;

  return;
}

realT
PenaltyCoulombIntermediate::StiffnessProjected(const localIndex index)
{
  const PenaltyCoulombIntermediateParameterData& matParams = *this->ParameterData(index);
  const PenaltyCoulombIntermediateStateData& matState = *this->StateData(index, 0);
  const realT normalApproachEstimate = matState.normalApproach +
                                       (matState.dxndt > 0 ? matState.dxndt * matState.dt : 0.0);
  return matParams.Stiffness(normalApproachEstimate);
}
