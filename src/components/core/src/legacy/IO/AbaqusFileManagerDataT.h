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
 * File: AbaqusFileManagerDataT.h
 * Class provides file IO data structure
 * created : SJ (11/16/2013)
 * ported from work by RRS (10/2001)
 */

#ifndef _ABAQUSFILEMANAGERDATAT_H_
#define _ABAQUSFILEMANAGERDATAT_H_

#include "FileManagerDataT.h"

namespace GPAC_IO
{
enum READMODE
{
  UNDEF, NODE, ELEM, NSET
};

/**
 * @author Scott Johnson
 * @brief Class to manager the data associated with the Abaqus geometry file
 * reader
 */
class AbaqusFileManagerDataT : public FileManagerDataT
{

public:

  AbaqusFileManagerDataT();

  void Reset();

  bool AdvanceLine(std::string& inputline);
  virtual void AddNodeLine(std::string& inputline, const realT geometryUnits, std::map<globalIndex,R1Tensor> &tempNodalPositionsMap, std::map<globalIndex,
                                                                                                                                              globalIndex> &resortNodeMap);
  void AddNodalPositionLineDiscrete(std::string& inputline, const realT geometryUnits);
  virtual void AddNodalPositionLine(std::string& inputline, const realT geometryUnits);
  virtual void AddNodalPositionLine(std::string& inputline, const realT geometryUnits, std::map<globalIndex,R1Tensor> &tempNodalPositionsMap);
  void AddElementRegionLine(std::string& inputline);
  void GoToSection(READMODE mode);

  const static std::string nodeKey;
  const static std::string elemKey;
  const static std::string nodeSetKey;

  const static std::string typeKey;
  const static std::string elSetKey;

  READMODE mode;   // The current read mode
};
}
#endif //_ABAQUSFILEMANAGERDATAT_H_
