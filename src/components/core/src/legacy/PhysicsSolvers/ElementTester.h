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

#ifndef ELEMENTTESTER_H_
#define ELEMENTTESTER_H_

/*
 * ElementTester.h
 *
 *  Created on: Jun 7, 2011
 *      Author: walsh24
 */

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "ObjectManagers/FunctionManager.h"
#include "Common/Common.h"

#include "ElementLibrary/FiniteElement.h"

#include "Utilities/TrilinosUtilities.h"
#include "Utilities/RCVSparse.h"

#include "PhysicsSolverStrings.h"

#include <map>
#include <string>
#include <vector>


/// FEM Element tests
class ElementTester : public SolverBase
{
public:
  ElementTester( const std::string& name,
                 ProblemManagerT* const pm );
  ~ElementTester();

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition  ) {}

  void InitializeCommunications( PartitionBase& partition  ) {}

  virtual double TimeStep( const realT& time,
                           const realT& dt,
                           const int cycleNumber,
                           PhysicalDomainT * domain,
                           const array<string>& namesOfSolverRegions,
                           SpatialPartition& partition,
                           FractunatorBase* const fractunator );

  virtual void RegisterFields( PhysicalDomainT& domain );

  static const char* SolverName(){return "ElementTester";};


private:

  void TestRegion(const realT& time, const realT& dt,
                  ElementRegionT& elementRegion,
                  PhysicalDomainT& domain);
  void ReportCoordinates();
  void CheckHandedness();
  void CheckElement(FiniteElement<3>& feElement);

  // Element data
  array<R1Tensor> nodeCoords;
  const localIndex* nodeList;
  int numNodesPerElem;


  #if GPAC_MPI
  const Epetra_MpiComm* m_CommPtr;
  #else
  const Epetra_SerialComm* m_CommPtr;
  #endif

};



#endif /*ELEMENTTESTER_H_*/
