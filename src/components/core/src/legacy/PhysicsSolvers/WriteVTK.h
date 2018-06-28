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

/*
 * WriteVTK.h
 *
 *  Created on: May 27, 2015
 *      Author: stuartwalsh
 */

#ifndef SRC_PHYSICSSOLVERS_WRITEVTK_H_
#define SRC_PHYSICSSOLVERS_WRITEVTK_H_

#include "WriteFieldToFile.h"

/// Update a field using a function of other fields defined on the same object.
class WriteVTK : public WriteFieldToFile
{
public:
  WriteVTK(  const std::string& name,
             ProblemManagerT* const pm );
  virtual ~WriteVTK(){};

  void ReadXML( TICPP::HierarchicalDataNode* const hdn );


  double TimeStep( const realT& time,
                   const realT& dt,
                   const int cycleNumber,
                   PhysicalDomainT& domain,
                   const array<string>& namesOfSolverRegions,
                   SpatialPartition& partition,
                   FractunatorBase* const fractunator );


  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects
  static const char* SolverName(){return "WriteVTK";};

private:


};



#endif /* SRC_PHYSICSSOLVERS_WRITEVTK_H_ */
