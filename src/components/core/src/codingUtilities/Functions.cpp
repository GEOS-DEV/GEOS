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
 * @file   Functions.cpp
 * @author walsh24
 */

#include "Functions.hpp"

#include "IO/ticpp/TinyXMLParser.h"
//#include "ObjectManagers/ProblemManagerT.h"
#include "ObjectManagers/TableManager.h"
#include "Utilities.h"
#include "StringUtilities.h"


#include <map>
#include <cmath>
#include <mpi.h>

/**
 * @author walsh24
 * Base function class
 *
 **/
Function::Function(TICPP::HierarchicalDataNode* hdn,
                   const ProblemManagerT* const problemManager  ):
  m_name(hdn->GetAttributeString("name")){}


Function::Function(const std::string& name ):
  m_name(name){
  /** empty **/
}

////////////////////
//
// Constant Function
//

/**
 * @author walsh24
 * @brief A function to represent a constant scalar
 *
 **/
ConstantFunction::ConstantFunction(TICPP::HierarchicalDataNode* hdn,
                                   const ProblemManagerT* const pm):
  Function(hdn,pm){
  ReadXML(hdn);
}

/// XML: <ConstantFunction name="half" f="0.5">
void ConstantFunction::ReadXML(TICPP::HierarchicalDataNode* hdn){
  m_value = hdn->GetAttributeValue<realT>("f");
}

REGISTER_Function( ConstantFunction )

////////////////////
//
// Polynomial Function
//

/**
 * @author walsh24
 * @brief A function to represent a polynomial
 *
 **/
PolynomialFunction::PolynomialFunction(TICPP::HierarchicalDataNode* hdn,
                                       const ProblemManagerT* const pm):
  Function(hdn,pm){
  ReadXML(hdn);
}

/// XML: <PolynomialFunction name="halfXSquared" coeffs="0.0 0.0 0.5">
void PolynomialFunction::ReadXML(TICPP::HierarchicalDataNode* hdn){
  std::string coeffsStr = hdn->GetAttributeString("coeffs");
  array<string> csv = Tokenize(coeffsStr," ");
  coeffs.resize(csv.size() );
  for( size_t i =0 ; i < csv.size() ; ++i )
    coeffs[i] = fromString<realT>(csv[i]);
}

realT PolynomialFunction::operator() (const realT& x) {
  realT rv = coeffs.back();

  for(int i = coeffs.size()-2 ; i >= 0 ; --i )
  {
    rv *= x;
    rv += coeffs[i];
  }

  return rv;
}

REGISTER_Function( PolynomialFunction )

////////////////////
//
// UniformRandomDistribution
//

/**
 * @author walsh24
 * @brief A function to represent numbers drawn from a uniform random
 * distribution in a specified range
 *
 *
 **/
UniformRandomDistribution::UniformRandomDistribution(TICPP::HierarchicalDataNode* hdn,
                                                     const ProblemManagerT* const pm):
  Function(hdn,pm){
  ReadXML(hdn);
}

/// XML: <UniformRandomDistribution name="Default" min="0.0" max="1.0">
void UniformRandomDistribution::ReadXML(TICPP::HierarchicalDataNode* hdn){
  min = hdn->GetAttributeOrDefault("min",0.0);
  realT max = hdn->GetAttributeOrDefault("max",1.0);
  df =  max-min;
  uniqueAcrossProcessors= hdn->GetAttributeOrDefault<bool>("uniqueAcrossProcessors",true);

  rank =0;
  #if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
}

realT UniformRandomDistribution::operator() (const realT& x  ) {
  unsigned int newSeed(0);
  realT rv;
  if(uniqueAcrossProcessors)
  {
    newSeed = rand(); // record a new seed to resync the random number
                      // generators (or not if out of sync)
    srand(rand() + rank); // "uncorrelated" random seed for next number
    rv = min + df*rand()/(RAND_MAX+1.0);
    srand(newSeed); // resync random seed
  }
  else
  {
    rv = min + df*rand()/(RAND_MAX+1.0);
  }

  return rv;
}

REGISTER_Function( UniformRandomDistribution )

////////////////////
//
// SymbolicFunction
//

/**
 * @author walsh24
 * @brief  Class to evaluate a function string using the warp function parser
 **/
SymbolicFunction::SymbolicFunction(TICPP::HierarchicalDataNode* hdn,
                                   const ProblemManagerT* const pm):
  Function(hdn,pm){
  ReadXML(hdn);
}

/**
 * @author walsh24
 * XML input
 * f = the function string;
 * constants = a comma separated list of pairs of constants and their values;
 * variables = a space separated list of variables used in the function, in the
 * order they will be passed to f when called
 * eg <SymbolicFunction name="someFunction"
 *                      f="1.e-5*d*k+(-b+sqrt(b^2 - 4*a*c))/(2*a)"
 *                      constants="a 1, b 2, k 2.345"
 *                      variables="c d"/>
 **/
void SymbolicFunction::ReadXML(TICPP::HierarchicalDataNode* hdn){

  std::string fStr = hdn->GetAttributeString("f");
  std::string Vars = hdn->GetAttributeStringOrDefault("variables","");
  std::string constStr = hdn->GetAttributeStringOrDefault("constants","");

  if(constStr !="")
  {
    std::vector<std::string> constants =  Tokenize(constStr,",");
    for(std::vector<std::string>::size_type i =0 ; i < constants.size() ; ++i)
    {
      Trim(constants[i]);
      std::vector<std::string> keyValue = Split(constants[i]," ");
      if(keyValue.size()>1)
      {
        std::string& constant = keyValue[0];
        double value = fromString<double>(keyValue[1]);
        m_fParser.AddConstant(constant,value);
      }
      else if(keyValue.size()==1)
      {
        std::string& constant = keyValue[0];
        throw GPException("SymbolicFunction: Undefined constant " + constant + ".");
      }
    }
  }

  int err  = m_fParser.Parse(fStr.c_str(), Vars);
  if(err >= 0)
    throw GPException("Error detected in SymbolicFunction expression: '"+ fStr+ "' at character " + toString(err+1) + "." );
  m_fParser.Optimize();
}


REGISTER_Function( SymbolicFunction )

///////////////////////
//
// 4D Table lookup
//

/**
 * @author walsh24
 */
Lookup4DTable::Lookup4DTable(TICPP::HierarchicalDataNode* hdn,
                             const ProblemManagerT* const pm):
  Function(hdn,pm){
  ReadXML(hdn);
}

/// XML:
///  <Functions>
///    <Lookup4DTable name="porosityFunc" table="porosityTable"/>
///  </Functions>
void Lookup4DTable::ReadXML(TICPP::HierarchicalDataNode* hdn){
  std::string tableName = hdn->GetAttributeString("table");

  const TableManager& tableManager = TableManager::Instance();
  const std::map<std::string,Table4D >& spatialTables = tableManager.Tables<4>();
  if(isMember(tableName,spatialTables))
  {
    std::map<std::string,Table4D >::const_iterator stItr = spatialTables.find(tableName);
    m_tablePtr = &(stItr->second);
  }
  else
  {
    throw GPException("Lookup4DTable: Table '"+ tableName+ "' was not found." );
  }
}

REGISTER_Function( Lookup4DTable )

///////////////////////
//
// 3D Table lookup
//

/**
 * @author walsh24
 */
Lookup3DTable::Lookup3DTable(TICPP::HierarchicalDataNode* hdn,
                             const ProblemManagerT* const pm):
  Function(hdn,pm){
  ReadXML(hdn);
}

/// XML:
///  <Functions>
///    <Lookup3DTable name="porosityFunc" table="porosityTable"/>
///  </Functions>
void Lookup3DTable::ReadXML(TICPP::HierarchicalDataNode* hdn){
  std::string tableName = hdn->GetAttributeString("table");

  const TableManager& tableManager = TableManager::Instance();
  const std::map<std::string,Table3D >& spatialTables = tableManager.Tables<3>();
  if(isMember(tableName,spatialTables))
  {
    std::map<std::string,Table3D >::const_iterator stItr = spatialTables.find(tableName);
    m_tablePtr = &(stItr->second);
  }
  else
  {
    throw GPException("Lookup3DTable: Table '"+ tableName+ "' was not found." );
  }
}

//realT Lookup3DTable::operator() (const realT& x) {
//const realT* xPtr = &x;
//R1Tensor X(xPtr[0],xPtr[1],xPtr[2]);
//	return m_tablePtr->lookup(&x);
//}

REGISTER_Function( Lookup3DTable )

///////////////////////
//
// 2D Table lookup
//

/**
 * @author walsh24
 */
Lookup2DTable::Lookup2DTable(TICPP::HierarchicalDataNode* hdn,
                             const ProblemManagerT* const pm):
  Function(hdn,pm){
  ReadXML(hdn);
}

/// XML:
///  <Functions>
///    <Lookup2DTable name="appertureFunc" table="appertureTable"/>
///  </Functions>
void Lookup2DTable::ReadXML(TICPP::HierarchicalDataNode* hdn){
  std::string tableName = hdn->GetAttributeString("table");

  const TableManager& tableManager = TableManager::Instance();
  const std::map<std::string,Table2D >& spatialTables = tableManager.Tables<2>();
  if(isMember(tableName,spatialTables))
  {
    std::map<std::string,Table2D >::const_iterator stItr = spatialTables.find(tableName);
    m_tablePtr = &(stItr->second);
  }
  else
  {
    throw GPException("Lookup2DTable: Table '"+ tableName+ "' was not found." );
  }
}

/*realT Lookup2DTable::operator() (const realT& x) {
   //const realT* xPtr = &x;
   //R1Tensor X(xPtr[0],xPtr[1],0);
   return m_tablePtr->lookup(&x);
   }
 */

REGISTER_Function( Lookup2DTable )

///////////////////////
//
// 1D Table lookup
//

/**
 * @author walsh24
 */
Lookup1DTable::Lookup1DTable(TICPP::HierarchicalDataNode* hdn,
                             const ProblemManagerT* const pm):
  Function(hdn,pm){
  ReadXML(hdn);
}

/// XML:
///  <Functions>
///    <Lookup1DTable name="signalLookup" table="signalTable"/>
///  </Functions>
void Lookup1DTable::ReadXML(TICPP::HierarchicalDataNode* hdn){
  std::string tableName = hdn->GetAttributeString("table");

  const TableManager& tableManager = TableManager::Instance();
  const std::map<std::string,Table1D >& timeTables = tableManager.Tables<1>();
  if(isMember(tableName,timeTables))
  {
    std::map<std::string,Table1D >::const_iterator ttItr = timeTables.find(tableName);
    m_tablePtr = &(ttItr->second);
  }
  else
  {
    throw GPException("Lookup1DTable: Table '"+ tableName+ "' was not found." );
  }
}

REGISTER_Function( Lookup1DTable )

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


/*
 * @brief One time evaluation of a real number expression (eg.
 *"1-sin(0.5+3.14159/6.0)")
 *
 */
realT EvaluateStringFunction(const std::string& fString){
  FunctionParser fParser;
  int err = fParser.Parse(fString.c_str(), "");
  if(err >=0)
  {
    throw GPException("Error detected in EvaluateStringFunction expression: '"+ fString+ "' at character " + toString(err+1) + "." );
  }
  const double Dummy=0;
  return fParser.Eval(&Dummy);
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////
//
// Function factory

typedef std::map<std::string, FunctionInitializer*> FunctionCatalogueType;

FunctionCatalogueType & getFunctionCatalogue(){
  static FunctionCatalogueType theCatalogue;
  return theCatalogue;
}

void getFunctionNames( std::vector<std::string>& nameList){
  for(FunctionCatalogueType::const_iterator it = getFunctionCatalogue().begin() ;
      it != getFunctionCatalogue().end() ; ++it)
  {
    nameList.push_back(it->first);
  }
}

Function* newFunction(const std::string& FunctionName, TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm)
{

  FunctionInitializer* FunctionInitializer = getFunctionCatalogue()[FunctionName];
  Function *theNewFunction = NULL;

  if(!FunctionInitializer)
    throw GPException("Could not create unrecognized Function"+ FunctionName);

  theNewFunction = FunctionInitializer->initializeFunction( hdn,pm );


  return theNewFunction;
}
