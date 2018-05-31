// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/**
 * File: EllipsoidFileManagerDataT.cpp
 * Class provides file IO data structure
 * created : SJ (11/16/2012)
 * ported from work by RRS (10/2001)
 */

#include "EllipsoidFileManagerDataT.h"

namespace GPAC_IO
{
EllipsoidFileManagerDataT::EllipsoidFileManagerDataT(): FileManagerDataT()
{}

realT EllipsoidFileManagerDataT::ReadLine(const std::string& inputline,
                                          globalIndex& globalNodeNumber,
                                          R1Tensor& position)
{
  std::istringstream linestream(inputline);
  R1Tensor principalRadii;
  linestream >> globalNodeNumber >> position(0) >> position(1) >> position(2) >>
  principalRadii(0) >> principalRadii(1) >> principalRadii(2);
  return principalRadii.MaxVal();
}

void EllipsoidFileManagerDataT::ReadLine(const std::string& inputline,
                                         globalIndex& globalNodeNumber,
                                         R1Tensor& position,
                                         R1Tensor& principalRadii,
                                         R1TensorT<4>& rotation)
{
  std::istringstream linestream(inputline);
  linestream >> globalNodeNumber >> position(0) >> position(1) >> position(2) >>
  principalRadii(0) >> principalRadii(1) >> principalRadii(2) >>
  rotation(0) >> rotation(1) >> rotation(2) >> rotation(3);
}

bool EllipsoidFileManagerDataT::AdvanceLine(std::string& inputline)
{
  FileManagerDataT::AdvanceLine(inputline);
  replace(inputline.begin(), inputline.end(), ',', ' ');
  return (inputline.compare(0, 1, "*") == 0 || inputline.size() == 0);
}

void EllipsoidFileManagerDataT::AddNodalPositionLine(std::string& inputline, const realT geometryUnits)
{
  R1Tensor nodePosition;
  globalIndex globalNodeNumber;
  realT rad = 0.0;
  {
    rad = ReadLine(inputline, globalNodeNumber, nodePosition);
    nodePosition *= geometryUnits;   // rescale to problem units
  }
  if (globalNodeNumber > maxGlobalNodeID)
  {
    maxGlobalNodeID = globalNodeNumber;
    nodalPositions.resize(maxGlobalNodeID + 1);
    isNodeInDomain.resize(maxGlobalNodeID + 1);
  }

  isNodeInDomain[globalNodeNumber] = 0;
  nodalPositions[globalNodeNumber] = nodePosition;

  R1Tensor tpos(nodePosition);
  tpos += rad;
  spatialMax.SetMax(tpos);
  tpos -= 2*rad;
  spatialMin.SetMin(tpos);
}
}
