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

#ifndef PHYSICSSOLVERSTRINGS_H_
#define PHYSICSSOLVERSTRINGS_H_

#include <string>
/*
 * PhysicsSolverStrings.h
 *
 *  Created on: Feb 18, 2011
 *      Author: walsh24
 */

/// Field and Model Parameter Strings
namespace PS_STR
{

const static std::string HeadStr = "Head";
const static std::string PressureStr = "Pressure";
const static std::string DarcyFluxStr = "DarcyFlux";
const static std::string FluidVelocityStr = "FluidVelocity";

const static std::string PermeabilityStr = "Permeability";
const static std::string DiffusivityStr = "Diffusivity";

const static std::string ConcentrationStr = "Concentration";
const static std::string ConcentrationFluxStr = "ConcentrationFlux";
const static std::string ReactionRateStr = "ReactionRate";

const static std::string PorosityStr = "Porosity";
const static std::string ApertureStr = "Aperture";
const static std::string MinimumApertureStr = "MinimumAperture";
const static std::string MaximumApertureStr = "MaximumAperture";

const static std::string BulkModulusStr = "BulkModulus";

const static std::string VolumetricFluxStr =  "VolumetricFlux";

const static std::string ElementCenterStr =  "ElementCenter";

const static std::string FaceAreaStr =  "FaceArea";
const static std::string FaceCenterStr =  "FaceCenter";
const static std::string FaceNormalStr =  "FaceNormal";

const static std::string EdgeCenterStr =  "EdgeCenter";
const static std::string EdgeLengthStr =  "EdgeLength";

const static std::string VerboseStr = "Verbose";

const static std::string ViscosityStr = "Viscosity";

const static std::string ProppantVolumeFractionStr = "ProppantVolumeFraction";
const static std::string ProppantPackVolumeFractionStr = "ProppantPackVolumeFraction";
}

#endif /*PHYSICSSOLVERSTRINGS_H_*/
