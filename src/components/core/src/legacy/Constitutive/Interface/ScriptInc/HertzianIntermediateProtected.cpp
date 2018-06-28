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
virtual_realT
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

//FUNCTION_BEGIN_PARSE
virtual_void
HertzianIntermediate::UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                                      InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  const HertzianIntermediateParameterData& matParams = static_cast < const HertzianIntermediateParameterData& > ( matParamsBase );
  HertzianIntermediateStateData& matState = static_cast < HertzianIntermediateStateData& > ( matStateBase );

  matState.mu = matParams.mu0;
}

//FUNCTION_BEGIN_PARSE
virtual_realT
HertzianIntermediate::ShearStrength(const InterfaceBaseParameterData&,
                                    InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  HertzianIntermediateStateData& matState = static_cast < HertzianIntermediateStateData& > ( matStateBase );
  return matState.stress * matState.mu;
}
