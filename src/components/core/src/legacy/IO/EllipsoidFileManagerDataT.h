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
 * File: EllipsoidFileManagerDataT.h
 * Class provides file IO data structure
 * created : SJ (11/16/2013)
 * ported from work by RRS (10/2001)
 */

#ifndef _ELLIPSOIDFILEMANAGERDATAT_H_
#define _ELLIPSOIDFILEMANAGERDATAT_H_

#include "FileManagerDataT.h"

namespace GPAC_IO
{
/**
 * @author Scott Johnson
 * @brief Class to manager the data associated with the Abaqus geometry file
 * reader
 */
class EllipsoidFileManagerDataT : public FileManagerDataT
{
public:
  EllipsoidFileManagerDataT();

  void AddNodalPositionLine(std::string& inputline, const realT geometryUnits);
  void AddElementRegionLine(std::string&){}
  bool AdvanceLine(std::string& inputline);

  static void ReadLine(const std::string& inputline,
                       globalIndex& globalNodeNumber,
                       R1Tensor& position,
                       R1Tensor& principalRadii,
                       R1TensorT<4>& rotation);
  static realT ReadLine(const std::string& inputline,
                        globalIndex& globalNodeNumber,
                        R1Tensor& position);
};
}
#endif //_ELLIPSOIDFILEMANAGERDATAT_H_
