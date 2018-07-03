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
 * File: FileManagerDataT.h
 * Class provides file IO data structure
 * created : SJ (11/16/2013)
 * ported from work by RRS (10/2001)
 */

#ifndef _FILEMANAGERDATAT_H_
#define _FILEMANAGERDATAT_H_

#include <fstream>
#include <string>
#include <stdio.h>
#include "Common/Common.h"

namespace GPAC_IO
{
/**
 * @author Scott Johnson
 * @brief Class to manager the data associated with the Abaqus geometry file
 * reader
 */
class FileManagerDataT
{

public:

  FileManagerDataT();
  ~FileManagerDataT(){}

  //file operations
  bool OpenFile();
  bool OK() {return !geometry.eof();}
  void CloseFile() {geometry.close();}
  void Write();

  //manage state
  virtual void Reset();

  virtual bool AdvanceLine(std::string& inputline);
  virtual void AddNodalPositionLine(std::string&, const realT){}
  virtual void AddElementRegionLine(std::string&){}

  static std::string ExtractValue(const std::string& key, const std::string& inputline);

  bool exist;   // Does the manager exist - set ReadMesh
  localIndex numNodes;   // Number of nodes - init (A) set (B)
  globalIndex maxGlobalNodeID;   // maximum nodal ID - set (A)
  globalIndex maxGlobalElemID;   // maximum element ID
  globalIndex totalNumElem;   //total number of elements
  localIndex numElementRegions;
  localIndex numNodeSets;   // Number of nodes - init (A) set (B)
  std::map<std::string, localIndex> numElementsInRegion;   // Number of elements
                                                           // in a region
  std::map<std::string, std::string> elementRegionNameTypes;   // Map of element
                                                               // region names
                                                               // to their types
  array<string> elementRegionNames;
  array<string> elementRegionTypes;
  gArray2d elemToNodes;

  array<R1Tensor> nodalPositions;
  std::map<globalIndex,R1Tensor> nodalPositionsMap;
  std::map<globalIndex,R1Tensor> potentialNodalPositionsMap;
  std::map<globalIndex,gArray1d> elemToNodesMap;
  gSet neighborList;


  R1Tensor spatialMin;
  R1Tensor spatialMax;
  // these arrays indicate if the objects status in the computational domain.
  // The values are:
  //  0 = Not in the domain
  //  1 = In the domain
  array<integer> isElemInDomain;
  array<integer> isNodeInDomain;
  std::map<std::string, gArray1d> elemsInRegion;
  std::map<std::string, globalIndex> maxNumElementsInRegion;
  std::streampos beginElems;
  std::streampos beginNodes;
  std::streampos beginNodesets;
  std::streampos current;
  lvector GlobalToLocalNodeMap;
  std::ifstream geometry;        /** geometry file stream */
  std::string filename;
};
}
#endif //_FILEMANAGERDATAT_H_
