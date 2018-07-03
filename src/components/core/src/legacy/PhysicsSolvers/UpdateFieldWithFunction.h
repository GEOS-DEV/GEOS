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
 * @file UpdateFieldWithFunction.h
 * @author walsh24
 * @date June 19, 2011
 */

#ifndef UPDATEFIELDWFUNCTION_H_
#define UPDATEFIELDWFUNCTION_H_

#include "SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"

#include "Common/Common.h"


/// Update a field using a function of other fields defined on the same object.
class UpdateFieldWithFunction : public SolverBase
{
public:
  UpdateFieldWithFunction(  const std::string& name,
                            ProblemManagerT* const pm );
  virtual ~UpdateFieldWithFunction();

  void ReadXML( TICPP::HierarchicalDataNode* const hdn );
  void RegisterFields( PhysicalDomainT& domain );

  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition  ) {}

  void InitializeCommunications( PartitionBase& partition  ) {}


  double TimeStep( const realT& time,
                   const realT& dt,
                   const int cycleNumber,
                   PhysicalDomainT& domain,
                   const array<string>& namesOfSolverRegions, SpatialPartition& partition,
                   FractunatorBase* const fractunator );

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects
  static const char* SolverName(){return "UpdateFieldWithFunction";};

private:

  std::string m_functionName;
  array<string> m_variables;
  array<FieldType> m_variable_types;

  array<string> m_setNames;

  PhysicalDomainT::ObjectDataStructureKeys m_objectKey;
  std::string m_regionName;  // only used if setting field in an element region
  std::string m_fieldName;
  FieldType m_fieldType;
  int m_component;

};


#endif /* UPDATEFIELDWFUNCTION_H_ */
