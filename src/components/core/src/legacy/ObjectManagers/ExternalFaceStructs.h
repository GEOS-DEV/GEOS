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

/**
 * @file ExternalFaceStructs.h
 * @author Scott Johnson
 * @date created on December 21, 2011
 */

#ifndef EXTERNALFACESTRUCTS_H_
#define EXTERNALFACESTRUCTS_H_

#include "Common/Common.h"

struct ExternalFaceStruct
{
  R1Tensor xmin, xmax;
  realT area;
  array<R1Tensor> xs, dxs;
};

struct CommonPlaneStruct
{
  R1Tensor centerCommonPlane, normalCommonPlane;
  array<R1Tensor> pointsCommonPlane;
  realT areaCommonPlane, penetration;
};

struct ToleranceStruct
{
  realT spatial, area, cosMin, penetration, feParentSolution;
  realT searchRadiusFactor, searchRadiusVelocityFactor, maximumSeparation;
};

#endif /* EXTERNALFACESTRUCTS_H_ */
