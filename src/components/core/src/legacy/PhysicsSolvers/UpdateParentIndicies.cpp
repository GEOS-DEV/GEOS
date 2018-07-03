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
 * @file UpdateParentIndicies.cpp
 * @author walsh24
 * @date December 3, 2013
 */

#include "SolverFactory.h"
#include "UpdateParentIndicies.h"
#include "Common/Common.h"
#include "Common/typedefs.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"



//////////////////////////////////////////////////////////////////////////////////////////

// Upate field with function

UpdateParentIndicies::UpdateParentIndicies(  const std::string& name,
                                             ProblemManagerT* const pm ):
  SolverBase(name,pm)
{}

UpdateParentIndicies::~UpdateParentIndicies()
{
  // TODO Auto-generated destructor stub
}

/*
 * <UpdateParentIndicies name="upi"            * name of solver
 *            object="Face"  />                * object to calculate the global
 * parent indicies on (Default=Node)
 */
void UpdateParentIndicies::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML( hdn );

  m_setNames = hdn->GetStringVector("setnames");
  if(m_setNames.empty())
    m_setNames = hdn->GetStringVector("setname");


  std::string objectTypeStr = hdn->GetAttributeStringOrDefault("object", PhysicalDomainT::FiniteElementNodeManagerStr() );
  {
    std::string oldStr = hdn->GetAttributeStringOrDefault("objecttype", "" ); // throw
                                                                              // error
                                                                              // if
                                                                              // using
                                                                              // old
                                                                              // syntax
    if(!oldStr.empty())
    {
      throw GPException("UpdateParentIndicies: Attempting to set objecttype - use 'object' instead.");
    }
  }
  m_objectKey = PhysicalDomainT::GetObjectDataStructureConditionKey(objectTypeStr);
//=======
//  std::string objectTypeStr = hdn->GetAttributeString("objecttype");
//  if(objectTypeStr.empty())
//    throw GPException("Must define an object type for updating a field with a
// function");
//  m_objectKey =
// PhysicalDomainT::GetObjectDataStructureConditionKey(objectTypeStr);
//>>>>>>> .r1185

  m_regionName = hdn->GetAttributeStringOrDefault("regionname","");



}

void UpdateParentIndicies::RegisterFields( PhysicalDomainT& domain )
{

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);

  objectManager.AddKeylessDataField( FieldInfo::globalIndexField,  "parentGlobalIndex", true, true );

}



/**
 *
 *
 **/

double UpdateParentIndicies::TimeStep( const realT& time,
                                       const realT& dt,
                                       const int cycleNumber,
                                       PhysicalDomainT& domain,
                                       const array<string>& namesOfSolverRegions,
                                       SpatialPartition& partition,
                                       FractunatorBase* const fractunator )
{
  realT dt_return = dt;

  m_stabledt.m_maxdt = 0.9*std::numeric_limits<double>::max();

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);
  gArray1d& parentIndexes = objectManager.GetFieldData<globalIndex>("parentGlobalIndex");

  if(m_setNames.empty())
  {
    for( array<string>::size_type i =0 ; i < objectManager.DataLengths() ; ++i)
    {
      localIndex localPI =  objectManager.GetParentIndex(i);
      parentIndexes[i] =  objectManager.m_localToGlobalMap[localPI];
    }
  }
  else
  {
    for( array<string>::size_type i =0 ; i < m_setNames.size() ; ++i)
    {
      lSet& set = objectManager.GetSet(m_setNames[i]);
      for( lSet::iterator si = set.begin() ; si != set.end() ; ++si)
      {
        localIndex localPI =  objectManager.GetParentIndex(*si);
        parentIndexes[*si] =  objectManager.m_localToGlobalMap[localPI];
      }

    }
  }
  return dt_return;
}


REGISTER_SOLVER( UpdateParentIndicies )
