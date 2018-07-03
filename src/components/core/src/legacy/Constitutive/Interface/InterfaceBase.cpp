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
 * InterfaceBase.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */


#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "InterfaceBase.h"
#include "../../IO/ticpp/HierarchicalDataNode.h.old"

InterfaceBase::InterfaceBase( const int paramSize, const int stateSize ):
  ConstitutiveBase(),
  m_paramSize( paramSize ),
  m_stateSize( stateSize )
{}

InterfaceBase::~InterfaceBase()
{}

void
InterfaceBaseStateData::UpdateOrientation(const realT normalApproachIn,
                                          const realT dtIn,
                                          const R1Tensor& normal,
                                          const R1Tensor& velocity,
                                          R1Tensor& dShearSlip,
                                          R1Tensor& shearSlip)
{
  dt = dtIn;
  normalApproach = normalApproachIn;

  dxndt = Dot(velocity, normal);

  //transform the tangential stress and slip
  const realT shearSlipNormal = Dot(shearSlip, normal);
  if (!isZero(shearSlipNormal))
  {
    //project shear slip and shear stress to current coordinates
    const realT shearSlipMagnitude = shearSlip.L2_Norm();
    R1Tensor slipNormal(normal);
    slipNormal *= shearSlipNormal;
    shearSlip -= slipNormal;
    shearSlip.Normalize();
    stressShearVector = shearSlip;
    shearSlip *= shearSlipMagnitude;
  }
  else
  {
    stressShearVector = shearSlip;
    stressShearVector.Normalize();
  }
  stressShearVector *= -stressShear;

  //update the contact velocity
  GeometryUtilities::OrthogonalVectorComponent(velocity, normal, dShearSlip);

  //next statement assumes ever-increasing slip ... should handle reversal of
  // slip, too!
  dxsdt = dShearSlip.L2_Norm();
  dShearSlip *= dt;

  //update the slip magnitude
  xs = shearSlip.L2_Norm();
}
void
InterfaceBase::Initialize( const localIndex index,
                           const realT stressNormal,
                           const realT stressShear)
{
  throw GPException("Cannot call InterfaceBase::Initialize\n");
}

void
InterfaceBase::UpdateProperties(const localIndex, std::map<std::string, realT>&, std::map<std::string, realT>& )
{}

realT
InterfaceBase::SetPermeabilityTerm(const localIndex index)
{
  InterfaceBaseStateData& matState = *this->StateData(index, 0);
  const InterfaceBaseParameterData& matParams = *this->ParameterData(index);

  if(matState.normalApproach <= 0)
  {
    matState.kappa = matParams.kappa0 - (matState.normalApproach * matState.normalApproach * matState.normalApproach);
  }
  else
  {
    matState.kappa = matParams.kappa0 +
                     matParams.dkappadnFct *
                     (isEqual(matParams.dkappadnExp, 1.0) ? matState.normalApproach : pow(matState.normalApproach, matParams.dkappadnExp)) +
                     matParams.dkappadsFct * (isEqual(matParams.dkappadsExp, 1.0) ? matState.xs : pow(matState.xs, matParams.dkappadsExp));
  }
  return matState.kappa;
}

realT
InterfaceBase::StiffnessProjected( const localIndex index )
{
  InterfaceBaseStateData& matState = *this->StateData(index, 0);
  const InterfaceBaseParameterData& matParams = *this->ParameterData(index);

  const realT normalApproachEstimate = matState.normalApproach +
                                       (matState.dxndt > 0 ? matState.dxndt * matState.dt : 0.0);

  return NormalStiffness(matParams, matState, normalApproachEstimate, false);
}

void
InterfaceBase::StrainDrivenUpdate( const localIndex index )
{
  InterfaceBaseStateData& matState = *this->StateData(index, 0);
  const InterfaceBaseParameterData& matParams = *this->ParameterData(index);

  //evolve the normal force ("stress" ... bad terminology, but I like the
  // inheritance)
  const realT stiffness = NormalStiffness(matParams, matState, matState.normalApproach, true);

  //evolve the friction
  UpdateFriction(matParams, matState);

  //evolve the shear stress
  {
    matState.stressShear += -stiffness * matParams.ston * matState.dxsdt * matState.dt;
    matState.xs += matState.dxsdt * matState.dt;
  }

  //determine whether strength exceeded and return to failure surface if
  // necessary
  ThresholdToFailureSurface(matParams, matState, matState.dxsdt * matState.dt);

  //update vector representation of shear stress
  //NOTE: JointBaseStateData::UpdateOrientation sets the direction as -(shear
  // slip direction)
  matState.stressShearVector.Normalize();
  matState.stressShearVector *= matState.stressShear;

  return;
}
realT InterfaceBase::ShearStrength(const InterfaceBaseParameterData& matParams,
                                   InterfaceBaseStateData& matState) const
{
  return matState.stress * matParams.mu0;
}

void InterfaceBase::UpdateFriction(const InterfaceBaseParameterData&,
                                   InterfaceBaseStateData&) const
{}

realT
InterfaceBase::NormalStiffness(const InterfaceBaseParameterData&,
                               InterfaceBaseStateData&,
                               const realT,
                               const bool ) const
{
  return 0.0;
}

void
InterfaceBase::ThresholdToFailureSurface(const InterfaceBaseParameterData& matParams,
                                         InterfaceBaseStateData& matState,
                                         const realT dxs)
{
  //determine whether strength exceeded and return to failure surface if
  // necessary
  const realT strength = ShearStrength( matParams, matState);
  if(strength < fabs(matState.stressShear))
  {
    matState.DissipatedEnergy += fabs(dxs) * strength;
    matState.stressShear = matState.stressShear < 0 ? -strength : strength;
  }
}
