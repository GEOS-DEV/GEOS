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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file BoundaryConditions.cpp
 * @author walsh24
 * @date December 5, 2011
 */

#include "BoundaryConditions.h"

#include "../../codingUtilities/Functions.hpp"
#include "ObjectManagers/FunctionManager.h"
#include "ObjectManagers/PhysicalDomainT.h"
//#include "ObjectManagers/ProblemManagerT.h"
#include "ObjectManagers/TableManager.h"
#include "Utilities/Utilities.h"
#include "Utilities/FieldTypeMultiPtr.h"
#include "PhysicsSolvers/PhysicsSolverStrings.h"

#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"


#include "MPI_Communications/SpatialPartition.h"


////////////////////////////
//
/// Boundary condition base class
/**
 * @author walsh24
 * @brief Boundary condition base class - should not be instantiated directly.
 *
 **/
BoundaryConditionBase::BoundaryConditionBase( TICPP::HierarchicalDataNode* hdn,
                                              const ProblemManagerT* const pm  )
{
  ReadXML(hdn);
}

BoundaryConditionBase::~BoundaryConditionBase()
{
  // TODO Auto-generated destructor stub
}


/**
 * @author walsh24
 * @brief Records object and field data.
 *
 * <GENERIC_BOUNDARY_CONDITION  ...
 *      object="Element"               * Object type, Node,Face,Element etc
 *(Default = Node)
 *      regionname="Region1"           * Apply to all elements in Region1 (for
 * Element fields only - no effect otherwise)
 *      setname="someFieldSubset"      * Subset to apply initial condition to -
 * leave out or empty to apply to whole field
 *      fieldname="Concentration"      * Name of field
 *      fieldtype="Scalar"/>           * Type of field (Default = Scalar)
 *
 *
 **/
void BoundaryConditionBase::ReadXML(TICPP::HierarchicalDataNode* hdn){

  m_fieldName = hdn->GetAttributeStringOrDefault("fieldname","");
  {
    array<string> tempSetName;
    tempSetName = hdn->GetStringVector("setname");
    if (!tempSetName.empty())
      throw GPException("ERROR!!! setname is not supported anymore in boundary condition.  Use setnames instead.");
  }
  m_setNames = hdn->GetStringVector("setnames");

  if (streq(m_fieldName, "combinedFlowRate") && m_setNames.size() > 1)
    throw GPException("ERROR!!! One combinedFlowRate BC is only supposed to be applied to one set.");

  m_component = hdn->GetAttributeOrDefault<int>("component",-1);

  //get object type
  {
    std::string objectTypeStr = hdn->GetAttributeStringOrDefault("object","Node");
    if(objectTypeStr.empty())
      throw GPException("Must define the object type to which to apply a boundary condition!");
    m_objectKey = PhysicalDomainT::GetObjectDataStructureConditionKey(objectTypeStr);
  }


  m_regionName = hdn->GetAttributeStringOrDefault("toregion","");

  // set direction
  {
    std::string temp = hdn->GetAttributeString("direction");
    if( !temp.empty() )
    {
      m_direction.StrVal( temp );
      m_direction /= m_direction.L2_Norm();

      m_component = -1;
      // check if direction aligned with coordinate system
      for( int i=0 ; i<m_direction.Length() ; ++i )
      {
        if( m_direction[i] >= ( 1.0 - 1.0e-14 ) )
          m_component = i;
      }
    }
  }

  m_isConstantInSpace = true;
  m_isConstantInTime = true;


  m_time = -1e64; // unlikely value

  m_timeFactor = hdn->GetAttributeOrDefault("timeFactor",-1.0);

  m_startTime = hdn->GetAttributeOrDefault("startTime",-std::numeric_limits<realT>::max());
  m_endTime = hdn->GetAttributeOrDefault("endTime",std::numeric_limits<realT>::max());

  m_option = hdn->GetAttributeOrDefault("option",0);
}


////////////////////////////
//
/// Simple Boundary condition
/**
 * @author walsh24
 * @brief Basic boundary condition.
 *
 **/

SimpleBoundaryCondition::SimpleBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  BoundaryConditionBase(hdn,pm)
{
  ReadXML(hdn);
}

void SimpleBoundaryCondition::ReadXML( TICPP::HierarchicalDataNode* BCNode)
{

  m_fieldName         = BCNode->GetAttributeString("fieldname");  // BC field
  m_timeTableName     = BCNode->GetAttributeString("timetable");
  m_functionName      = BCNode->GetAttributeString("function");
  m_time = -1e64; // unlikely value
  m_scale             = BCNode->GetAttributeOrDefault<realT>("scale",1.0);
  m_value = m_scale;

  // get function
  if(!(m_functionName.empty()))
  {
    FunctionManager& fm = FunctionManager::Instance();
    m_function = &(fm.GetFunction(m_functionName));
  }

  m_isConstantInSpace = true;
  m_isConstantInTime = (  !(m_timeTableName.empty())  && !(m_functionName.empty()) );

}

realT SimpleBoundaryCondition::GetValue(const ObjectDataStructureBaseT& object,
                                        const lSet::const_iterator& si,
                                        realT time )
{

  if (!isEqual(m_time, time))
  {
    if (!(m_timeTableName.empty()))
    {
      array<real64> t(1);
      t[0] = time;
      const realT tableval = TableManager::Instance().LookupTable<1>(m_timeTableName, t);
      m_value = m_scale * tableval;
    }
    else if (!(m_functionName.empty()))
    {
      m_value = m_scale * (*m_function)(time);
    }
    m_time = time;
  }
  return m_value;
}


REGISTER_BoundaryCondition( SimpleBoundaryCondition )

////////////////////////////////////////////
//
/// Boundary condition function
/**
 * @author walsh24
 * @brief Boundary condition as a function of the fields on the object where the
 * boundary condition is defined (and time).
 *
 **/

BoundaryConditionFunction::BoundaryConditionFunction( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  BoundaryConditionBase(hdn,pm)
{
  ReadXML(hdn);
}

void BoundaryConditionFunction::ReadXML( TICPP::HierarchicalDataNode* hdn)
{

  m_fieldName         = hdn->GetAttributeString("fieldname");  // BC field
  m_functionName      = hdn->GetAttributeString("function");

  // get function
  if(!(m_functionName.empty()))
  {
    FunctionManager& fm = FunctionManager::Instance();
    m_function = &(fm.GetFunction(m_functionName));
  }

  std::string varsStr = hdn->GetAttributeStringOrDefault("variables","");
  if(varsStr.empty())
  {
    m_variableNames.resize(0);
  }
  else
  {
    m_variableNames = Tokenize(varsStr," ");
  }

  std::string varTypesStr = hdn->GetAttributeStringOrDefault("variableTypes","");
  if( varTypesStr.empty() )
  {
    // assume all variables are scalars by default
    m_variableTypes = array<FieldType>(m_variableNames.size(),FieldInfo::realField);
  }
  else
  {
    array<string> vTypesVect = Tokenize(varTypesStr," ");
    m_variableTypes.resize(vTypesVect.size());

    if(m_variableTypes.size() != m_variableNames.size())
      throw GPException("Error InitialConditionFunction: Number of variable types not equal to number of variables.");

    for( array<string>::size_type i=0 ; i < vTypesVect.size() ; ++i )
      m_variableTypes[i] = fromString<FieldType>(vTypesVect[i]);

  }
  m_nVars = m_variableNames.size();
  m_fieldPtrs = std::vector<FieldTypeMultiPtr>(m_nVars );

  m_isConstantInTime = hdn->GetAttributeOrDefault<bool>("constantInTime",false);  // very
                                                                                  // likely
                                                                                  // not
  m_isConstantInSpace =  hdn->GetAttributeOrDefault<bool>("constantInSpace",false);  // very
                                                                                     // likely
                                                                                     // not
  m_scale =  hdn->GetAttributeOrDefault<realT>("scale", 1.0);
  m_offset =  hdn->GetAttributeOrDefault<realT>("offset", 0.0);

  // x vector for function
  int xLength = 0;
  for(unsigned j =0 ; j < m_nVars ; ++j)
  {
    if(streq(m_variableNames[j],"time"))
    {
      xLength += 1;
      //   m_fieldPtrs[j].SetFieldPtr(&m_time); // must be pointer to m_time
      // (otherwise not preserved) ***
    }
    else
    {
      xLength += FieldSize( m_variableTypes[j] );
      //    m_fieldPtrs[j].SetFieldPtr(object, m_variableTypes[j],
      // m_variableNames[j]); ***
    }
  }
  m_x = std::vector<realT>(xLength );
  m_time = -1e-64;

  // *** this may be possible if pointers are preserved - but suspect that
  // fields are not registered
}

realT BoundaryConditionFunction::GetValue(const ObjectDataStructureBaseT& object,
                                          const lSet::const_iterator& si,
                                          realT time )
{

  if (!isEqual(m_time, time))
  {
    m_time = time;   // will need to keep this even if pointers init moved to
                     // read xml

    // only need to do this once (in reality may only need to do this once per
    // run
    // but not sure if object pointer addresses are preserved or if objects have
    // been initiated.

    ObjectDataStructureBaseT* ObjPtr = const_cast< ObjectDataStructureBaseT * >( &object );   // forgive
                                                                                              // me
                                                                                              // father

    // initialize field pointers
    for(unsigned j =0 ; j < m_nVars ; ++j)
    {
      if(streq(m_variableNames[j],"time"))
      {
        m_fieldPtrs[j].SetFieldPtr(&m_time);   // must be pointer to m_time
                                               // (otherwise not preserved)
      }
      else
      {
        m_fieldPtrs[j].SetFieldPtr(*ObjPtr, m_variableTypes[j], m_variableNames[j]);
      }
    }

  }


/*
   int nVars = m_variableNames.size();
   std::vector<FieldTypeMultiPtr> fieldPtr(nVars );

   /FieldTypeMultiPtr theFieldPtr;
   theFieldPtr.SetFieldPtr(object , fieldType, m_fieldName);

   int xLength = 0;
   for(int j =0; j < nVars; ++j){
    if(streq(variables[j],"time")){
      xLength += 1;
      fieldPtr[j].SetFieldPtr(&time);
    } else {
      xLength += FieldSize( m_variable_types[j] );
      fieldPtr[j].SetFieldPtr(*this, m_variable_types[j], m_variableNames[j]);
    }
   }
 */


  // pack field values into x
  std::vector<realT>::iterator xItr = m_x.begin();
  for(unsigned j =0 ; j < m_nVars ; ++j)
    xItr = m_fieldPtrs[j].CopyValues(*si, xItr);

  // calculate value
  m_value =  (*m_function)(m_x[0]); // nb needs to be stored (in case constant
                                    // in space)
  m_value *= m_scale;
  m_value += m_offset;

  return m_value;
}


REGISTER_BoundaryCondition( BoundaryConditionFunction )


///////////////////////////////////////////

/**

   varname="Pressure Temperature"
   allocated="true"

 **/

MultiVarBoundaryConditionBase::MultiVarBoundaryConditionBase(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  BoundaryConditionBase(hdn,pm)
{
  ReadXML(hdn);

}

void MultiVarBoundaryConditionBase::ReadXML(TICPP::HierarchicalDataNode* hdn)
{

  m_varName = hdn->GetStringVector("varName");

  if(m_varName.empty())
    throw GPException("ERROR: varName is not set for MultiVarBoundaryCondition!");

  m_time = -1e64; // unlikely value

  m_varName.resize(m_varName.size());

}

/**

   isClamped="1"
   timetables="P1Table T1Table"

 **/


MultiVarDirichletBoundaryCondition::MultiVarDirichletBoundaryCondition(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  MultiVarBoundaryConditionBase(hdn,pm)
{
  ReadXML(hdn);

}

void MultiVarDirichletBoundaryCondition::ReadXML(TICPP::HierarchicalDataNode* hdn)
{
  m_isClamped = hdn->GetAttributeOrDefault<bool>("isClamped",true);

  if(!m_isClamped)
  {
    m_tables = hdn->GetStringVector("timeTables");
    if(m_tables.size() == 0)
      throw GPException("ERROR: timeTables is not set for MultiVarDirichletBoundaryCondition!");
  }
}

void MultiVarDirichletBoundaryCondition::CheckVars(const array<string> &varName)
{
  if(!m_isClamped)
  {
    array<string> tmp;
    for(array<string>::const_iterator it = varName.begin() ; it != varName.end() ; ++it)
    {
      bool notFound = true;
      array<string>::const_iterator itt = m_tables.begin();
      for(array<string>::const_iterator itm = m_varName.begin() ; itm != m_varName.end() ; ++itm, ++itt)
      {
        if(streq(*it, *itm))
        {
          tmp.push_back(*itt);
          notFound = false;
          break;
        }
      }
      if(notFound)
        throw GPException("ERROR: " + *it + "is not found for MultiVarDirichletBoundaryCondition!");
    }

    m_tables.resize(tmp.size());
    m_tables = tmp;
    m_value.resize(tmp.size());
  }
}

const array<real64>& MultiVarDirichletBoundaryCondition::GetValues(realT time)
{
  if(!m_isClamped && !isEqual(m_time,time) )
  {
    m_value.resize(m_tables.size());

    array<real64> t(1, time);
    array<real64>::iterator itv = m_value.begin();
    for(array<string>::const_iterator it = m_tables.begin() ; it != m_tables.end() ; ++it, ++itv)
      *itv = TableManager::Instance().LookupTable<1>(*it, t);
    m_time = time;
  }
  return m_value;
}

REGISTER_BoundaryCondition( MultiVarDirichletBoundaryCondition )


/**

   allocedByWeight="1"
   timetables="P1Table T1Table"

 **/


MultiVarSrcFluxBoundaryCondition::MultiVarSrcFluxBoundaryCondition(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  MultiVarBoundaryConditionBase(hdn,pm)
{
  ReadXML(hdn);

}

void MultiVarSrcFluxBoundaryCondition::ReadXML(TICPP::HierarchicalDataNode* hdn)
{

  m_allocedByWeight = hdn->GetAttributeOrDefault<bool>("allocedByWeight",true);

  m_tables = hdn->GetStringVector("timeTables");
  if(m_tables.size() == 0)
    throw GPException("ERROR: timeTables is not set for MultiVarSrcFluxBoundaryCondition!");
}


void MultiVarSrcFluxBoundaryCondition::CheckVars(const array<string> &varName)
{
  array<string> tmp;
  for(array<string>::const_iterator it = varName.begin() ; it != varName.end() ; ++it)
  {
    bool notFound = true;
    array<string>::const_iterator itt = m_tables.begin();
    for(array<string>::const_iterator itm = m_varName.begin() ; itm != m_varName.end() ; ++itm, ++itt)
    {
      if(streq(*it, *itm))
      {
        tmp.push_back(*itt);
        notFound = false;
        break;
      }
    }
    if(notFound)
      throw GPException("ERROR: " + *it + "is not found for MultiVarDirichletBoundaryCondition!");
  }

  m_tables.resize(tmp.size());
  m_tables = tmp;
  m_value.resize(tmp.size());
}

const array<real64>& MultiVarSrcFluxBoundaryCondition::GetValues(realT time)
{
  if(!isEqual(m_time,time) )
  {
    array<real64> t(1, time);
    array<real64>::iterator itv = m_value.begin();
    for(array<string>::const_iterator it = m_tables.begin() ; it != m_tables.end() ; ++it, ++itv)
      *itv = TableManager::Instance().LookupTable<1>(*it, t);
    m_time = time;
  }
  return m_value;
}

REGISTER_BoundaryCondition( MultiVarSrcFluxBoundaryCondition )


////////////////////////////
//
/// Traction Boundary condition
/**
 * @author walsh24
 * @brief Traction boundary condition
 *
 **/

TractionBoundaryCondition::TractionBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  SimpleBoundaryCondition(hdn,pm)
{
  ReadXML(hdn);
  m_objectKey = PhysicalDomainT::FiniteElementFaceManager;
  m_fieldName = "Traction";
  m_fieldType = FieldInfo::R2TensorField;
  m_useNormalFlag = hdn->GetAttributeOrDefault<bool>("applyNormalTraction",false);

  m_isConstantInSpace = true;
  if(m_useNormalFlag)
  {
    m_isConstantInSpace = false;
  }
}

// Individual solvers should use this function rather than calculate the
// traction directly in the solver.
// This allows the traction bc to be developed separately from the individual
// solvers
R1Tensor TractionBoundaryCondition::GetTractionOnFace(PhysicalDomainT& domain, const lSet::const_iterator& fc,  realT& time){
  realT value = this->GetValue(domain.m_feFaceManager,fc,time);

  R1Tensor traction;
  if( this->IsNormalTraction() )
  {
    traction = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *fc );
    traction *= -1;
  }
  else
  {
    traction = this->GetDirection(time);
  }
  traction *= value;
  return traction;
}

REGISTER_BoundaryCondition( TractionBoundaryCondition )

////////////////////////////
//
/// Traction Boundary condition
/**
 * @author walsh24 fu4
 * @brief Traction boundary condition
 *
 **/

TractionBoundaryConditionFunction::TractionBoundaryConditionFunction( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  TractionBoundaryCondition(hdn,pm)
{
  ReadXML(hdn);
  m_objectKey = PhysicalDomainT::FiniteElementFaceManager;
  m_fieldName = "Traction";
  m_fieldType = FieldInfo::R2TensorField;
  m_useNormalFlag = hdn->GetAttributeOrDefault<bool>("applyNormalTraction",false);

}



void TractionBoundaryConditionFunction::ReadXML( TICPP::HierarchicalDataNode* hdn)
{

  m_functionName      = hdn->GetAttributeString("function");

  // get function
  if(!(m_functionName.empty()))
  {
    FunctionManager& fm = FunctionManager::Instance();
    m_function = &(fm.GetFunction(m_functionName));
  }

  std::string varsStr = hdn->GetAttributeStringOrDefault("variables","");
  if(varsStr.empty())
  {
    m_variableNames.resize(0);
  }
  else
  {
    m_variableNames = Tokenize(varsStr," ");
  }

  std::string varTypesStr = hdn->GetAttributeStringOrDefault("variableTypes","");
  if( varTypesStr.empty() )
  {
    // assume all variables are scalars by default
    m_variableTypes = array<FieldType>(m_variableNames.size(),FieldInfo::realField);
  }
  else
  {
    array<string> vTypesVect = Tokenize(varTypesStr," ");
    m_variableTypes.resize(vTypesVect.size());

    if(m_variableTypes.size() != m_variableNames.size())
      throw GPException("Error TractionBoundaryConditionFunction: Number of variable types not equal to number of variables.");

    for( array<string>::size_type i=0 ; i < vTypesVect.size() ; ++i )
      m_variableTypes[i] = fromString<FieldType>(vTypesVect[i]);

  }
  m_nVars = m_variableNames.size();
  m_fieldPtrs = std::vector<FieldTypeMultiPtr>(m_nVars );

  m_isConstantInTime = hdn->GetAttributeOrDefault<bool>("constantInTime",false);  // very
                                                                                  // likely
                                                                                  // not
  m_isConstantInSpace =  hdn->GetAttributeOrDefault<bool>("constantInSpace",false);  // very
                                                                                     // likely
                                                                                     // not
  m_scale =  hdn->GetAttributeOrDefault<realT>("scale", 1.0);
  m_offset =  hdn->GetAttributeOrDefault<realT>("offset", 0.0);

  // x vector for function
  int xLength = 0;
  for(unsigned j =0 ; j < m_nVars ; ++j)
  {
    if(streq(m_variableNames[j],"time"))
    {
      xLength += 1;
      //   m_fieldPtrs[j].SetFieldPtr(&m_time); // must be pointer to m_time
      // (otherwise not preserved) ***
    }
    else
    {
      xLength += FieldSize( m_variableTypes[j] );
      //    m_fieldPtrs[j].SetFieldPtr(object, m_variableTypes[j],
      // m_variableNames[j]); ***
    }
  }
  m_x = std::vector<realT>(xLength );
  m_time = -1e-64;

  // *** this may be possible if pointers are preserved - but suspect that
  // fields are not registered
}



R1Tensor TractionBoundaryConditionFunction::GetTractionOnFace(PhysicalDomainT& domain, const lSet::const_iterator& fc,  realT& time){

  R1Tensor traction;
  if( this->IsNormalTraction() )
  {
    traction = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *fc );
    traction *= -1;
  }
  else
  {
    traction = this->GetDirection(time);
  }



  if (!isEqual(m_time, time))
  {
    m_time = time;   // will need to keep this even if pointers init moved to
                     // read xml

    // only need to do this once (in reality may only need to do this once per
    // run
    // but not sure if object pointer addresses are preserved or if objects have
    // been initiated.


    // initialize field pointers
    for(unsigned j =0 ; j < m_nVars ; ++j)
    {
      if(streq(m_variableNames[j],"time"))
      {
        m_fieldPtrs[j].SetFieldPtr(&m_time);   // must be pointer to m_time
                                               // (otherwise not preserved)
      }
      else
      {
        m_fieldPtrs[j].SetFieldPtr(domain.m_feFaceManager, m_variableTypes[j], m_variableNames[j]);
      }
    }

  }

  // pack field values into x
  std::vector<realT>::iterator xItr = m_x.begin();
  for(unsigned j =0 ; j < m_nVars ; ++j)
    xItr = m_fieldPtrs[j].CopyValues(*fc, xItr);

  // calculate value
  m_value =  (*m_function)(m_x[0]); // nb needs to be stored (in case constant
                                    // in space)
  m_value *= m_scale;
  m_value += m_offset;

  traction *= m_value;

  return traction;
}



REGISTER_BoundaryCondition( TractionBoundaryConditionFunction )



////////////////////////////
//
/// Hydraulic Pressure Boundary condition
/**
 * @author walsh24
 * @brief Traction boundary condition, returns traction calculated from pressure
 * defined on face and normal direction.
 *
 **/

HydraulicPressureBoundaryCondition::HydraulicPressureBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  TractionBoundaryCondition(hdn,pm)
{
  ReadXML(hdn);
  m_objectKey = PhysicalDomainT::FiniteElementFaceManager;
  m_fieldName = "Traction";
  m_fieldType = FieldInfo::R2TensorField;
  m_useNormalFlag = true;
  m_isConstantInSpace = false;
}

// Individual solvers should use this function rather than calculate the
// traction directly in the solver.
// This allows the traction bc to be developed separately from the individual
// solvers
R1Tensor HydraulicPressureBoundaryCondition::GetTractionOnFace(PhysicalDomainT& domain, const lSet::const_iterator& fc,  realT& time){


  R1Tensor traction = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *fc );
  array<real64>& facePressures = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure> ();

  traction *= -facePressures[*fc];
  return traction;
}

REGISTER_BoundaryCondition( HydraulicPressureBoundaryCondition )

////////////////////////////
//
/// Uniform Pressure Boundary condition
/**
 * @author walsh24
 * @brief Traction boundary condition - applies uniform pressure based on bc
 * value and normal direction for all external faces
 * that are 1. not on the domain boundary or 2. flow faces on the domain
 * boundary
 *
 **/

UniformPressureBoundaryCondition::UniformPressureBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  TractionBoundaryCondition(hdn,pm)
{
  ReadXML(hdn);
  m_objectKey = PhysicalDomainT::FiniteElementFaceManager;
  m_fieldName = "UniformPressureBoundaryCondition";
  m_fieldType = FieldInfo::R2TensorField;
  m_useNormalFlag = true;
  m_isConstantInSpace = false;
}

// Individual solvers should use this function rather than calculate the
// traction directly in the solver.
// This allows the traction bc to be developed separately from the individual
// solvers
R1Tensor UniformPressureBoundaryCondition::GetTractionOnFace(PhysicalDomainT& domain,
                                                             const lSet::const_iterator& fc,
                                                             realT& time)
{

  R1Tensor traction;

  if (domain.m_feFaceManager.m_isExternal[*fc])
  {
    localIndex parentNodeIndex = domain.m_feFaceManager.m_parentIndex[*fc];
    if (parentNodeIndex == LOCALINDEX_MAX)
      parentNodeIndex = *fc;

    array<integer>& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
    if (flowFaceType[parentNodeIndex] == 1)
    { // is flow face

      realT value = this->GetValue(domain.m_feFaceManager, fc, time);
      traction = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, *fc);
      traction *= -value;

      /*
         array<real64>& facePressures =
            domain.m_feFaceManager.GetFieldData<FieldInfo::pressure> ();
         facePressures[*fc] = value;

         const array<real64>& fluidVolume  =
            domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
         array<real64>& faceFluidMass =
            domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
         array<real64>& faceFluidDensity =
            domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
       */

    }
  }
  return traction;
}

REGISTER_BoundaryCondition( UniformPressureBoundaryCondition )

////////////////////////////
//
/// Wall Boundary condition
/**
 * @author walsh24
 * @brief Wall boundary condition
 *
 **/

WallBoundaryCondition::WallBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  SimpleBoundaryCondition(hdn,pm)
{
  ReadXML(hdn);

  std::string temp    = hdn->GetAttributeString("wallcoord");
  if( !temp.empty() )
  {
    m_position.StrVal( temp );
  }
  m_objectKey = PhysicalDomainT::FiniteElementNodeManager;
}



REGISTER_BoundaryCondition( WallBoundaryCondition )

////////////////////////////
//
/// Outflow Boundary condition
/**
 * @author walsh24
 * @brief Outflow boundary condition
 *
 **/


OutflowBoundaryCondition::OutflowBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  SimpleBoundaryCondition(hdn,pm)
{
  ReadXML(hdn);
}

REGISTER_BoundaryCondition( OutflowBoundaryCondition )

////////////////////////////
//
/// Outflow Boundary condition
/**
 * @author walsh24
 * @brief Outflow boundary condition
 *
 **/


InflowBoundaryCondition::InflowBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  SimpleBoundaryCondition(hdn,pm)
{
  ReadXML(hdn);
}

REGISTER_BoundaryCondition( InflowBoundaryCondition )

////////////////////////////
//
/// Non-penetrating Boundary condition
/**
 * @author walsh24
 * @brief Prevent two surfaces from interpenetrating
 *
 **/

NonPenetratingBoundaryCondition::NonPenetratingBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  BoundaryConditionBase(hdn,pm)
{
  ReadXML(hdn);
}


void NonPenetratingBoundaryCondition::ReadXML( TICPP::HierarchicalDataNode* hdn)
{

  if(m_setNames.size() != 2)
  {
    throw GPException("Error - non penetrating boundary condition requires two sets.");
  }
  ;

  m_fieldName = NonPenetratingBoundaryCondition::BoundaryConditionName();
  m_objectKey = PhysicalDomainT::FiniteElementNodeManager;
  m_updatePressure = hdn->GetAttributeOrDefault<bool>("UpdatePressure",false);
}
void NonPenetratingBoundaryCondition::RegisterFields(PhysicalDomainT& domain )
{
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::displacement>();
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("FaceCenter");
}
realT NonPenetratingBoundaryCondition::GetValue(const ObjectDataStructureBaseT& object,
                                                const lSet::const_iterator& si,
                                                realT time )
{

  throw GPException("Error - non penetrating boundary condition queried for value.");
  return 0.0;
}

void NonPenetratingBoundaryCondition::UpdateNearestNeighborMaps(PhysicalDomainT& domain ){

  const array<R1Tensor>& u = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement> ();
  const array<R1Tensor>& X = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition> ();

  array<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>("FaceCenter");

  //nodes
  {
    std::map< std::string, lSet >::const_iterator setMapA = domain.m_feNodeManager.m_Sets.find( m_setNames[0] );
    std::map< std::string, lSet >::const_iterator setMapB = domain.m_feNodeManager.m_Sets.find( m_setNames[1] );

    if( setMapA != domain.m_feNodeManager.m_Sets.end() &&  setMapB != domain.m_feNodeManager.m_Sets.end() )
    {

      const lSet& setA = setMapA->second;
      const lSet& setB = setMapB->second;

      // fixme - brute force search for nearest neighbors.
      for( lSet::const_iterator ndA=setA.begin() ; ndA!=setA.end() ; ++ndA )
      {
        realT minSqrdDist = std::numeric_limits<realT>::max();
        localIndex nbr = 0;
        R1Tensor posA =  u[*ndA] + X[*ndA];
        for( lSet::const_iterator ndB=setB.begin() ; ndB!=setB.end() ; ++ndB )
        {
          R1Tensor l =  u[*ndB] + X[*ndB] - posA;

          realT ll = Dot(l,l);
          if(ll < minSqrdDist)
          {
            minSqrdDist = ll;
            nbr = *ndB;
          }
        }

        if(minSqrdDist < std::numeric_limits<realT>::max())
        {
          m_nearestNodeNeighborMap[*ndA] = nbr;

          if(isMember(nbr,m_nearestNodeNeighborMap) )
          {
            localIndex oldNdA = m_nearestNodeNeighborMap[nbr];
            R1Tensor lb =  u[nbr] + X[nbr] - u[oldNdA]-X[oldNdA];
            if(minSqrdDist < Dot(lb,lb) )
            {
              m_nearestNodeNeighborMap[nbr] = *ndA;
            }
          }
          else
          {

            m_nearestNodeNeighborMap[nbr] = *ndA;
          }
        }

      }

    }
  }

  //faces
  {
    std::map< std::string, lSet >::const_iterator setMapA = domain.m_feFaceManager.m_Sets.find( m_setNames[0] );
    std::map< std::string, lSet >::const_iterator setMapB = domain.m_feFaceManager.m_Sets.find( m_setNames[1] );

    if( setMapA != domain.m_feFaceManager.m_Sets.end() &&  setMapB != domain.m_feFaceManager.m_Sets.end() )
    {

      const lSet& setA = setMapA->second;
      const lSet& setB = setMapB->second;

      // update face centers
      for( lSet::const_iterator fcA=setA.begin() ; fcA!=setA.end() ; ++fcA )
      {
        domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager,*fcA, faceCenter(*fcA));
      }
      for( lSet::const_iterator fcB=setB.begin() ; fcB!=setB.end() ; ++fcB )
      {
        domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager,*fcB, faceCenter(*fcB));
      }

      // brute force search for nearest neighbors.
      for( lSet::const_iterator fcA=setA.begin() ; fcA!=setA.end() ; ++fcA )
      {
        realT minSqrdDist = std::numeric_limits<realT>::max();
        localIndex nbr = 0;
        R1Tensor posA =  faceCenter(*fcA);
        for( lSet::const_iterator fcB=setB.begin() ; fcB!=setB.end() ; ++fcB )
        {
          R1Tensor l =  faceCenter(*fcB) - posA;

          realT ll = Dot(l,l);
          if(ll < minSqrdDist)
          {
            minSqrdDist = ll;
            nbr = *fcB;
          }
        }

        if(minSqrdDist < std::numeric_limits<realT>::max())
        {
          m_nearestFaceNeighborMap[*fcA] = nbr;
          m_nearestFaceNeighborMap[nbr] = *fcA;
        }

      }

    }
  }
}

REGISTER_BoundaryCondition( NonPenetratingBoundaryCondition )

////////////////////////////
//
/// SinglePartitionPeriodicBoundaryCondition
/**
 * @author walsh24
 * @brief Single partition periodic boundary condition
 *
 **/

SinglePartitionPeriodicBoundaryCondition::SinglePartitionPeriodicBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  BoundaryConditionBase(hdn,pm)
{
  ReadXML(hdn);
}


void SinglePartitionPeriodicBoundaryCondition::ReadXML( TICPP::HierarchicalDataNode* hdn)
{

  if(m_setNames.size() != 2)
  {
    throw GPException("Error - periodic boundary condition requires two sets.");
  }
  ;

  m_dimension = hdn->GetAttributeValue<int>("dimension");

  m_fieldName = SinglePartitionPeriodicBoundaryCondition::BoundaryConditionName();
  // m_objectKey = PhysicalDomainT::NodeManager;
}
void SinglePartitionPeriodicBoundaryCondition::RegisterFields(PhysicalDomainT& domain ){

  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::displacement>();
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("FaceCenter");
  domain.m_feEdgeManager.AddKeylessDataField<R1Tensor>("EdgeCenter");
}
realT SinglePartitionPeriodicBoundaryCondition::GetValue(const ObjectDataStructureBaseT& object,
                                                         const lSet::const_iterator& si,
                                                         realT time )
{

  throw GPException("Error - periodic boundary condition queried for value.");

  return 0.0;
}

void SinglePartitionPeriodicBoundaryCondition::SetNeighborMaps(PhysicalDomainT& domain ){

  const array<R1Tensor>& refPositions = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition> ();
  array<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>("FaceCenter");
  array<R1Tensor>& edgeCenters = domain.m_feEdgeManager.GetFieldData<R1Tensor>("EdgeCenter");

  PlanarSorter planarSorterNodes(refPositions,m_dimension);
  PlanarSorter planarSorterFaces(faceCenters,m_dimension);
  PlanarSorter planarSorterEdges(edgeCenters,m_dimension);


  //nodes
  {
    std::vector< std::vector<localIndex> > sortedIndexes(2);

    for(int a =0 ; a < 2 ; ++a)
    {

      std::map< std::string, lSet >::const_iterator setMap = domain.m_feNodeManager.m_Sets.find( m_setNames[a] );
      if( setMap != domain.m_feNodeManager.m_Sets.end() )
      {
        const lSet& set = setMap->second;
        sortedIndexes[a].assign(set.begin(),set.end() );
        std::sort(sortedIndexes[a].begin(),sortedIndexes[a].end(),planarSorterNodes);
      }

    }

    if(sortedIndexes[0].size() != sortedIndexes[1].size())
    {
      throw GPException("SinglePartitionPeriodicBoundaryCondition::SetNeighborMaps: Size of " + m_setNames[0]
                        + " does not match size of " + m_setNames[1] + "\n");
    }
    else
    {
      // make maps
      for(unsigned int i =0 ; i < sortedIndexes[0].size() ; ++i)
      {
        m_nodeNeighborMapA[sortedIndexes[0][i]] = sortedIndexes[1][i];
        m_nodeNeighborMapB[sortedIndexes[1][i]] = sortedIndexes[0][i];
      }
    }
  }

  //faces
  {
    std::vector< std::vector<localIndex> > sortedIndexes(2);

    for(int a =0 ; a < 2 ; ++a)
    {

      std::map< std::string, lSet >::const_iterator setMap = domain.m_feFaceManager.m_Sets.find( m_setNames[a] );
      if( setMap != domain.m_feFaceManager.m_Sets.end() )
      {
        const lSet& set = setMap->second;

        // update face centers
        for( lSet::const_iterator kf=set.begin() ; kf!=set.end() ; ++kf )
          domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, *kf, faceCenters[*kf]);


        sortedIndexes[a].assign(set.begin(),set.end() );
        std::sort(sortedIndexes[a].begin(),sortedIndexes[a].end(),planarSorterFaces);
      }

    }

    if(sortedIndexes[0].size() != sortedIndexes[1].size())
    {
      throw GPException("SinglePartitionPeriodicBoundaryCondition::SetNeighborMaps: Size of " + m_setNames[0]
                        + " face set does not match size of " + m_setNames[1] + "\n");
    }
    else
    {
      // make maps
      for(unsigned int i =0 ; i < sortedIndexes[0].size() ; ++i)
      {
        m_faceNeighborMapA[sortedIndexes[0][i]] = sortedIndexes[1][i];
        m_faceNeighborMapB[sortedIndexes[1][i]] = sortedIndexes[0][i];
      }
    }
  }

  // edges
  {
    std::vector< std::vector<localIndex> > sortedIndexes(2);

    for(int a =0 ; a < 2 ; ++a)
    {

      std::map< std::string, lSet >::const_iterator setMap = domain.m_feEdgeManager.m_Sets.find( m_setNames[a] );
      if( setMap != domain.m_feEdgeManager.m_Sets.end() )
      {
        const lSet& set = setMap->second;

        // update edge centers
        for( lSet::const_iterator ke=set.begin() ; ke!=set.end() ; ++ke )
          domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, *ke, edgeCenters[*ke]);

        sortedIndexes[a].assign(set.begin(),set.end() );
        std::sort(sortedIndexes[a].begin(),sortedIndexes[a].end(),planarSorterEdges);
      }

    }

    if(sortedIndexes[0].size() != sortedIndexes[1].size())
    {
      throw GPException("SinglePartitionPeriodicBoundaryCondition::SetNeighborMaps: Size of " + m_setNames[0]
                        + " edge set does not match size of " + m_setNames[1] + "\n");
    }
    else
    {
      // make maps
      for(unsigned int i =0 ; i < sortedIndexes[0].size() ; ++i)
      {
        m_edgeNeighborMapA[sortedIndexes[0][i]] = sortedIndexes[1][i];
        m_edgeNeighborMapB[sortedIndexes[1][i]] = sortedIndexes[0][i];
      }
    }
  }

}

REGISTER_BoundaryCondition( SinglePartitionPeriodicBoundaryCondition )

/////////////////////////////
//
// SwitchBoundaryConditions
//
/**
 * @author walsh24
 * @brief Switch between boundary conditions on the same boundary
 *
 **/

/// Alternate between multiple boundary conditions
SwitchBoundaryConditions::SwitchBoundaryConditions( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  BoundaryConditionBase(hdn,pm),
  m_pmPtr(pm)
{
  ReadXML(hdn);
}


/**
 *  <SwitchBoundaryConditions function="timeSwitch"   * 0 for pressure bc, 1 for
 * dummy
 *                            setname="Zmax">
 *      <BoundaryCondition    fieldname="Pressure"
 *                            scale="1 GPa" />
 *      <BoundaryCondition    fieldname="Dummy"
 *                            scale="0.0" />
 *  </SwitchBoundaryConditions>
 *
 **/
void SwitchBoundaryConditions::ReadXML( TICPP::HierarchicalDataNode* hdn)
{

  m_functionName = hdn->GetAttributeString("function");
  m_default_index = 0;

  // get function
  FunctionManager& fm = FunctionManager::Instance();
  m_function = &(fm.GetFunction(m_functionName));

  m_isConstantInSpace = true;
  for (TICPP::HierarchicalDataNode* bcNode = hdn->Next(true) ; bcNode ; bcNode = hdn->Next())
  {
    std::string bcType = bcNode->Heading();

    BoundaryConditionBase* bcPtr = newBoundaryCondition(bcType, bcNode, m_pmPtr);
    m_children.push_back(bcPtr);
    m_isConstantInSpace = m_isConstantInSpace && bcPtr->m_isConstantInSpace;
  }

  m_isConstantInTime = false;

}


realT SwitchBoundaryConditions::GetValue(const ObjectDataStructureBaseT& object,
                                         const lSet::const_iterator& si,realT time){
  unsigned indx = (*m_function)(time);
  if(indx >= m_children.size())
    indx = m_default_index;

  return m_children[indx]->GetValue(object,si,time);
}


// provide pointer to underlying boundary conditions
// enables upcasting of switch boundary condition to underlying boundary
// condition type
BoundaryConditionBase* SwitchBoundaryConditions::GetActiveBCPointer(realT time){
  unsigned indx = (*m_function)(time);
  if(indx >= m_children.size())
    indx = m_default_index;

  return m_children[indx]->GetActiveBCPointer(time);
}

const std::string& SwitchBoundaryConditions::GetFieldName(realT time)
{
  unsigned indx = (*m_function)(time);
  if(indx > m_children.size())
    indx = m_default_index;

  return m_children[indx]->GetFieldName(time);
}

const FieldType& SwitchBoundaryConditions::GetFieldType(realT time)
{
  unsigned indx = (*m_function)(time);
  if(indx > m_children.size())
    indx = m_default_index;

  return m_children[indx]->GetFieldType(time);
}


int SwitchBoundaryConditions::GetComponent(realT time){
  unsigned indx = (*m_function)(time);
  if(indx > m_children.size())
    indx = m_default_index;
  return m_children[indx]->GetComponent(time);
}

const R1Tensor& SwitchBoundaryConditions::GetDirection(realT time){
  unsigned indx = (*m_function)(time);
  if(indx > m_children.size())
    indx = m_default_index;
  return m_children[indx]->GetDirection(time);
}

REGISTER_BoundaryCondition( SwitchBoundaryConditions )


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////////
//
// Boundary condition factory

typedef std::map<std::string, BoundaryConditionInitializer*> BoundaryConditionCatalogueType;

BoundaryConditionCatalogueType & getBoundaryConditionCatalogue(){
  static BoundaryConditionCatalogueType theCatalogue;
  return theCatalogue;
}

void getBoundaryConditionNames( std::vector<std::string>& nameList){
  for(BoundaryConditionCatalogueType::const_iterator it = getBoundaryConditionCatalogue().begin() ;
      it != getBoundaryConditionCatalogue().end() ; ++it)
  {
    nameList.push_back(it->first);
  }
  ;
}

BoundaryConditionBase* newBoundaryCondition(const std::string& BoundaryConditionName, TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm)
{

  BoundaryConditionInitializer* BoundaryConditionInitializer = getBoundaryConditionCatalogue()[BoundaryConditionName];
  BoundaryConditionBase *theNewBoundaryCondition = NULL;

  if(!BoundaryConditionInitializer)
    throw GPException("Could not create unrecognized BoundaryCondition "+ BoundaryConditionName);

  theNewBoundaryCondition = BoundaryConditionInitializer->initializeBoundaryCondition( hdn,pm );

  return theNewBoundaryCondition;
}
