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
 * File: AbaqusFileManagerDataT.cpp
 * Class provides file IO data structure
 * created : SJ (11/16/2012)
 * ported from work by RRS (10/2001)
 */

#include "AbaqusFileManagerDataT.h"

namespace GPAC_IO
{
const std::string AbaqusFileManagerDataT::nodeKey = "*NODE";
const std::string AbaqusFileManagerDataT::elemKey = "*ELEMENT";
const std::string AbaqusFileManagerDataT::nodeSetKey = "*NSET";

const std::string AbaqusFileManagerDataT::typeKey = "TYPE";
const std::string AbaqusFileManagerDataT::elSetKey = "ELSET";

/**
 * @brief Constructor
 * @author Scott Johnson
 */
AbaqusFileManagerDataT::AbaqusFileManagerDataT():
  FileManagerDataT(),
  mode(UNDEF)
{}

bool AbaqusFileManagerDataT::AdvanceLine(std::string& inputline)
{
  FileManagerDataT::AdvanceLine(inputline);

  replace(inputline.begin(), inputline.end(), ',', ' ');

  if (inputline.compare(0, nodeKey.size(), nodeKey) == 0)
  {
    if (std::streamoff(beginNodes) == 0)
      beginNodes = this->current;
    mode = NODE;
    return true;
  }
  else if (inputline.compare(0, elemKey.size(), elemKey) == 0)
  {
    if (std::streamoff(beginElems) == 0)
      beginElems = this->current;
    mode = ELEM;
    return true;
  }
  else if (inputline.compare(0, nodeSetKey.size(), nodeSetKey) == 0)
  {
    if (std::streamoff(beginNodesets) == 0)
      beginNodesets = this->current;
    mode = NSET;
    return true;
  }
  else if (inputline.compare(0, 1, "*") == 0 || inputline.size() == 0)
  {
    bool rval = false;
    if( mode != UNDEF )
    {
      rval = true;
    }
    mode = UNDEF;
    return rval;
  }
  return false;
}

//R.W. get maxGlobalNodeID, nodal position
void AbaqusFileManagerDataT::AddNodeLine(std::string& inputline, const realT geometryUnits, std::map<globalIndex,R1Tensor> &tempNodalPositionsMap,
                                         std::map<globalIndex,globalIndex> &resortNodeMap)
{
  std::istringstream linestream(inputline);
  globalIndex globalNodeNumber, newglobalNodeNumber;
  R1Tensor nodePosition;
  {
    linestream >> globalNodeNumber >> nodePosition(0) >> nodePosition(1) >> nodePosition(2);
    newglobalNodeNumber = resortNodeMap.size();
    resortNodeMap[globalNodeNumber] = newglobalNodeNumber;
    //globalNodeNumber = newglobalNodeNumber + 1;
    nodePosition*=geometryUnits;   // rescale to problem units
  }
  if (newglobalNodeNumber > maxGlobalNodeID)
  {
    maxGlobalNodeID = newglobalNodeNumber;
  }

  tempNodalPositionsMap[newglobalNodeNumber] = nodePosition;

}

//without MPI: just to get the basic spatial partition information
void AbaqusFileManagerDataT::AddNodalPositionLine(std::string& inputline, const realT geometryUnits)
{
  std::istringstream linestream(inputline);
  globalIndex globalNodeNumber;
  R1Tensor nodePosition;
  {
    linestream >> globalNodeNumber >> nodePosition(0) >> nodePosition(1) >> nodePosition(2);
    nodePosition*=geometryUnits;   // rescale to problem units
  }
  if (globalNodeNumber > maxGlobalNodeID)
  {
    maxGlobalNodeID = globalNodeNumber;
  }

  spatialMin.SetMin(nodePosition);
  spatialMax.SetMax(nodePosition);
}

void AbaqusFileManagerDataT::AddNodalPositionLineDiscrete(std::string& inputline, const realT geometryUnits)
{
  std::istringstream linestream(inputline);
  globalIndex globalNodeNumber;
  R1Tensor nodePosition;
  {
    linestream >> globalNodeNumber >> nodePosition(0) >> nodePosition(1) >> nodePosition(2);
    nodePosition*=geometryUnits;   // rescale to problem units
  }
  if (globalNodeNumber > maxGlobalNodeID)
  {
    maxGlobalNodeID = globalNodeNumber;
    nodalPositions.resize(maxGlobalNodeID + 1);
    isNodeInDomain.resize(maxGlobalNodeID + 1);
  }

  isNodeInDomain[globalNodeNumber] = 0;
  nodalPositions[globalNodeNumber] = nodePosition;

  spatialMin.SetMin(nodePosition);
  spatialMax.SetMax(nodePosition);
}

//with MPI: read certain number of nodes into nodalPositions with the size <=
// mpiNodeLimit
void AbaqusFileManagerDataT::AddNodalPositionLine(std::string& inputline, const realT geometryUnits, std::map<globalIndex,R1Tensor> &tempNodalPositionsMap)
{
  std::istringstream linestream(inputline);
  globalIndex globalNodeNumber;
  R1Tensor nodePosition;
  {
    linestream >> globalNodeNumber >> nodePosition(0) >> nodePosition(1) >> nodePosition(2);
    nodePosition*=geometryUnits;   // rescale to problem units
  }

  tempNodalPositionsMap[globalNodeNumber] = nodePosition;
}

void AbaqusFileManagerDataT::AddElementRegionLine(std::string& inputline)
{
  ++numElementRegions;
  elementRegionTypes.push_back(FileManagerDataT::ExtractValue(typeKey, inputline));
  elementRegionNames.push_back(FileManagerDataT::ExtractValue(elSetKey, inputline));
  numElementsInRegion[elementRegionNames.back()];
  elemsInRegion[elementRegionNames.back()];
  elementRegionNameTypes[elementRegionNames.back()] = elementRegionTypes.back();
}

void AbaqusFileManagerDataT::Reset()
{
  FileManagerDataT::Reset();

  mode = UNDEF;
}

void AbaqusFileManagerDataT::GoToSection(READMODE mode_)
{
  geometry.clear();
  switch(mode_)
  {
  case NODE:
    geometry.seekg(beginNodes);
    break;
  case ELEM:
    geometry.seekg(beginElems);
    break;
  case NSET:
    geometry.seekg(beginNodesets);
    break;
  default:
    break;
  }
}
}
