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

#ifndef WriteFieldToFile_H_
#define WriteFieldToFile_H_

#include "SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"

#include "Common/Common.h"

#include <deque>


/// Update a field using a function of other fields defined on the same object.
class WriteFieldToFile : public SolverBase
{
public:
  WriteFieldToFile(  const std::string& name,
                     ProblemManagerT* const pm );
  virtual ~WriteFieldToFile();

  void ReadXML( TICPP::HierarchicalDataNode* const hdn );
  void RegisterFields( PhysicalDomainT& domain );

  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition ) {}

  void InitializeCommunications( PartitionBase& partition  ) {}


  virtual double TimeStep( const realT& time,
                           const realT& dt,
                           const int cycleNumber,
                           PhysicalDomainT& domain,
                           const array<string>& namesOfSolverRegions,
                           SpatialPartition& partition,
                           FractunatorBase* const fractunator );

  void SetMaxStableTimeStep( const realT& time, PhysicalDomainT& domain, const array<string>& namesOfSolverRegions,
                             SpatialPartition& partition );

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects
  static const char* SolverName(){return "WriteFieldToFile";};

protected:

  std::string m_filePrefix;
  array<string> m_setNames;
  realT m_dt;   // used for regular output intervals
  std::deque<realT> m_outputTimes;  // preset output times

  PhysicalDomainT::ObjectDataStructureKeys m_objectType;
  std::string m_regionName;  // only used if setting field in an element region

  array<string> m_fieldNames;
  std::vector<FieldType> m_fieldTypes;

  std::string m_headerString;

  bool m_appendToFile;
  bool m_isFirstTime;

};


#endif /* UPDATEFIELDWFUNCTION_H_ */
