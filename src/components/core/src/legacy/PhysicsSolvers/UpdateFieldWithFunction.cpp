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
 * @file UpdateFieldWithFunction.cpp
 * @author walsh24
 * @date June 19, 2011
 */

#include "SolverFactory.h"
#include "UpdateFieldWithFunction.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"



//////////////////////////////////////////////////////////////////////////////////////////

// Upate field with function

UpdateFieldWithFunction::UpdateFieldWithFunction(  const std::string& name,
                                                   ProblemManagerT* const pm ):
  SolverBase(name,pm)
{}

UpdateFieldWithFunction::~UpdateFieldWithFunction()
{
  // TODO Auto-generated destructor stub
}

/*
 * <UpdateFieldWithFunction name="CaUpdate"       * name of the solver
 *            object="Face"                    * location of field to set
 *(Default=Node)
 *            fieldtype="Scalar"                    * type of field to set
 *   (Default=Scalar)
 *            fieldname="Ca"                       * name of field to update
 *            function="aFunction"                * name of the function
 *            variables="Ca CO2"                  * function variables
 *            variableTypes="Scalar Scalar" />   * variableTypes (assumed scalar
 * for all if omitted)
 */
void UpdateFieldWithFunction::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML( hdn );

  m_functionName = hdn->GetAttributeString("function");
  m_fieldName = hdn->GetAttributeString("fieldname");
  m_component = hdn->GetAttributeOrDefault("component",0);

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
      throw GPException("UpdateFieldWithFunction: Attempting to set objecttype - use 'object' instead.");
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

  std::string fieldTypeStr = hdn->GetAttributeStringOrDefault("fieldtype",FieldInfo::RealStr);
  m_fieldType = fromString<FieldType>(fieldTypeStr);

  m_regionName = hdn->GetAttributeStringOrDefault("regionname","");

  std::string varsStr = hdn->GetAttributeString("variables");
  m_variables = Tokenize(varsStr," ");
  std::string varTypesStr = hdn->GetAttributeStringOrDefault("variabletypes","");

  if(varTypesStr == "")
  {
    m_variable_types = array<FieldType>(m_variables.size(),FieldInfo::realField);
  }
  else
  {
    array<string> vTypesVect = Tokenize(varTypesStr," ");
    m_variable_types.resize(vTypesVect.size());

    if(m_variable_types.size() != m_variables.size())
      throw GPException("UpdateFieldWithFunction: Number of variable types not equal to number of variables.");

    for(size_t i =0 ; i < vTypesVect.size() ; ++i)
      m_variable_types[i] = fromString<FieldType>(vTypesVect[i]);

  }


}

void UpdateFieldWithFunction::RegisterFields( PhysicalDomainT& domain )
{

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);

  objectManager.AddKeylessDataField( m_fieldType,  m_fieldName, true, true );
  // register variables
  for(size_t i =0 ; i < m_variables.size() ; ++i)
  {
    if(!streq(m_variables[i],"dt") && !streq(m_variables[i],"time") )
    {
      objectManager.AddKeylessDataField( m_variable_types[i],  m_variables[i], true, true );
    }
  }
}



/**
 *
 *
 **/

double UpdateFieldWithFunction::TimeStep( const realT& time,
                                          const realT& dt,
                                          const int cycleNumber,
                                          PhysicalDomainT& domain,
                                          const array<string>& namesOfSolverRegions,
                                          SpatialPartition& partition,
                                          FractunatorBase* const fractunator )
{
  realT dt_return = dt;

  m_stabledt.m_maxdt = std::numeric_limits<double>::max();

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);
  if(m_setNames.empty())
  {
    objectManager.SetFieldEqualToFunction(m_fieldType,  m_fieldName, m_functionName, m_variables, m_variable_types, m_component,time,dt);
  }
  else
  {
    for( array<string>::size_type i =0 ; i < m_setNames.size() ; ++i)
    {
      lSet& set = objectManager.GetSet(m_setNames[i]);
      objectManager.SetFieldEqualToFunction(m_fieldType,  m_fieldName, m_functionName, m_variables, m_variable_types,set, m_component,time,dt);
    }
  }
  return dt_return;
}


REGISTER_SOLVER( UpdateFieldWithFunction )
