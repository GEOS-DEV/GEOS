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

//FUNCTION_BEGIN_PARSE
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
