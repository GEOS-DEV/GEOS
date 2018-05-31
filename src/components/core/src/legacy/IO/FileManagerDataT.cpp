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
 * File: FileManagerDataT.cpp
 * Class provides file IO data structure
 * created : SJ (11/16/2012)
 * ported from work by RRS (10/2001)
 */

#include "FileManagerDataT.h"

namespace GPAC_IO
{
/**
 * @brief Constructor
 * @author Scott Johnson
 */
FileManagerDataT::FileManagerDataT():
  exist(false),
  numNodes(0),
  maxGlobalNodeID(0),
  maxGlobalElemID(0),
  totalNumElem(0),
  numElementRegions(0),
  numNodeSets(0),
  numElementsInRegion(),
  elementRegionNameTypes(),
  elementRegionNames(),
  elementRegionTypes(),
  elemToNodes(),
  nodalPositions(),
  spatialMin(),
  spatialMax(),
  isElemInDomain(),
  isNodeInDomain(),
  elemsInRegion(),
  beginElems(0),
  beginNodes(0),
  beginNodesets(0),
  current(0),
  GlobalToLocalNodeMap(),
  geometry(),
  filename()
{
  spatialMin = 1e100;
  spatialMax = -1e100;
}

/**
 * @brief Brief
 * @author Scott Johnson
 * @param[in] geometryfilename File name for the Abaqus geometry file
 * @return Flag for successfully opening the file
 */
bool FileManagerDataT::OpenFile()
{
  if (filename.empty())
  {
    this->exist = false;
  }
  else
  {
    geometry.open(filename.c_str());
    this->exist = geometry.is_open();
    if(!this->exist)
      throw GPException("Could not open mesh file!");
  }
  return this->exist;
}

/**
 * @brief Write the current state to std out
 * @author Scott Johnson
 */
void FileManagerDataT::Write()
{
  std::cout << "elementRegionNameTypes " << elementRegionNameTypes.size() << "\n";
  for (std::map<std::string, std::string>::const_iterator it = elementRegionNameTypes.begin() ;
       it != elementRegionNameTypes.end() ; ++it)
    std::cout << " " << it->first << " " << it->second << "\n";
  std::cout << "nodalPositions " << nodalPositions.size() << "\n";
  for (array<R1Tensor>::size_type i = 0 ; i < nodalPositions.size() ; ++i)
    std::cout << " " << nodalPositions[i](0) << " " << nodalPositions[i](1) << " "
              << nodalPositions[i](2) << "\n";
  std::cout << "maxGlobalNodeID " << maxGlobalNodeID << "\n";
  std::cout << "numElementsInRegion " << numElementsInRegion.size() << "\n";
  for (std::map<std::string, localIndex>::const_iterator it = numElementsInRegion.begin() ;
       it != numElementsInRegion.end() ; ++it)
    std::cout << " " << it->first << " " << it->second << "\n";
  std::cout << "elemsInRegion " << elemsInRegion.size() << "\n";
  for (std::map<std::string, gArray1d>::const_iterator it = elemsInRegion.begin() ;
       it != elemsInRegion.end() ; ++it)
    for (gArray1d::size_type j = 0 ; j < it->second.size() ; ++j)
      std::cout << " " << it->first << " " << it->second[j] << "\n";
  std::cout << "maxGlobalElemID " << maxGlobalElemID << "\n";
  std::cout << "isNodeInDomain " << isNodeInDomain.size() << "\n";
  for (array<integer>::size_type i = 0 ; i < isNodeInDomain.size() ; ++i)
    std::cout << " " << isNodeInDomain[i] << "\n";
  std::cout << "isElemInDomain " << isElemInDomain.size() << "\n";
  for (array<integer>::size_type i = 0 ; i < isElemInDomain.size() ; ++i)
    std::cout << " " << isElemInDomain[i] << "\n";
  std::cout << "numNodes " << numNodes << "\n";
  std::cout << "elemToNodes " << elemToNodes.size() << "\n";
  //    for(int i = 0; i < elemToNodes.size(); ++i)
  //      std::cout << " " << elemToNodes[i,j] << "\n";
  std::cout << "numElementRegions " << numElementRegions << "\n";
  std::cout << "elementRegionNames " << elementRegionNames.size() << "\n";
  for (array<string>::size_type i = 0 ; i < elementRegionNames.size() ; ++i)
    std::cout << " " << elementRegionNames[i] << "\n";
  std::cout << "elementRegionTypes " << elementRegionTypes.size() << "\n";
  for (array<string>::size_type i = 0 ; i < elementRegionTypes.size() ; ++i)
    std::cout << " " << elementRegionTypes[i] << "\n";
  std::cout << "elementRegionNameTypes " << elementRegionNameTypes.size() << "\n";
  for (std::map<std::string, std::string>::const_iterator it = elementRegionNameTypes.begin() ;
       it != elementRegionNameTypes.end() ; ++it)
    std::cout << " " << it->first << " " << it->second << "\n";
  std::cout << "elemsInRegion " << elemsInRegion.size() << "\n";
  for (std::map<std::string, gArray1d>::const_iterator it = elemsInRegion.begin() ;
       it != elemsInRegion.end() ; ++it)
    for (gArray1d::size_type i = 0 ; i < it->second.size() ; ++i)
      std::cout << " " << it->first << " " << it->second[i] << "\n";
  std::cout << "GlobalToLocalNodeMap " << GlobalToLocalNodeMap.size() << "\n";
  for (lvector::size_type i = 0 ; i < GlobalToLocalNodeMap.size() ; ++i)
    std::cout << " " << GlobalToLocalNodeMap[i] << "\n";
}

void FileManagerDataT::Reset()
{
  numNodes = 0;   // Number of nodes
  maxGlobalNodeID = 0;
  nodalPositions.resize(1);
  isNodeInDomain.resize(1);

  maxGlobalElemID = 0;
  totalNumElem = 0;
  numElementRegions = 0;
}

bool FileManagerDataT::AdvanceLine(std::string& inputline)
{
  current = geometry.tellg();
  getline(geometry, inputline);
  return false;
}

//for parsing format: "key=value,"
std::string FileManagerDataT::ExtractValue(const std::string& key, const std::string& inputline)
{
  const int pos1 = inputline.find(key) + key.length();
  std::string first = inputline.substr(pos1);

  const int pos2 = first.find("=") + 1;
  std::string second = first.substr(pos2);

  const int pos3 = second.find(" ");
  return pos3 > -1 ? second.substr(0, pos3) : second;
}
}
