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

//FUNCTION_BEGIN_PARSE
virtual_realT InterfaceBase::ShearStrength(const InterfaceBaseParameterData& matParams,
                                           InterfaceBaseStateData& matState) const
{
  return matState.stress * matParams.mu0;
}

//FUNCTION_BEGIN_PARSE
virtual_void InterfaceBase::UpdateFriction(const InterfaceBaseParameterData&,
                                           InterfaceBaseStateData&) const
{}

//FUNCTION_BEGIN_PARSE
virtual_realT
InterfaceBase::NormalStiffness(const InterfaceBaseParameterData&,
                               InterfaceBaseStateData&,
                               const realT,
                               const bool ) const
{
  return 0.0;
}

//FUNCTION_BEGIN_PARSE
virtual_void
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
