/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#ifndef FUNCTIONS_H
#define FUNCTIONS_H


#include "common/DataTypes.hpp"

#if USE_FPARSER==1
#include "fparser.hh"
#endif

#include <string>

class ProblemManagerT;
namespace TICPP
{
class HierarchicalDataNode;
}

/// Base class for all user-defined functions
class Function
{

public:
  Function(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  Function(const std::string& name);
  virtual ~Function()
  {}

  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn)=0;
  virtual realT operator()(const realT& x)
  {
    return 0.0;
  }

  std::string Name()
  {
    return m_name;
  }

  /// returns name of specific function instance
  /// derived classes must define "static const char* FunctionName()"
  /// which returns name of derived class.

private:
  std::string m_name;
};

/// Class representing a constant function
class ConstantFunction : public Function
{

public:
  ConstantFunction(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  ~ConstantFunction()
  {}

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    return m_value;
  }

  static const char* FunctionName()
  {
    return "ConstantFunction";
  }

private:
  realT m_value;
};

/// Class representing a polynomial function
class PolynomialFunction : public Function
{

public:
  PolynomialFunction(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  ~PolynomialFunction()
  {}

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x);

  static const char* FunctionName()
  {
    return "PolynomialFunction";
  }

private:
  array1d<real64> coeffs;
};

/// Class representing a uniform random distribution
class UniformRandomDistribution : public Function
{

public:
  UniformRandomDistribution(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  ~UniformRandomDistribution()
  {}

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x);

  static const char* FunctionName()
  {
    return "UniformRandomDistribution";
  }

private:
  realT df;
  realT min;
  bool uniqueAcrossProcessors;
  int rank;
};

#if USE_FPARSER==1
/// Function class for user-defined functions
class SymbolicFunction : public Function
{

public:
  SymbolicFunction(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  ~SymbolicFunction()
  {}

  void ReadXML(TICPP::HierarchicalDataNode* hdn);

  realT operator()(const realT& x)
  {
    return m_fParser.Eval(&x);
  }

  static const char* FunctionName()
  {
    return "SymbolicFunction";
  }

private:
  FunctionParser m_fParser;
};
#endif


#include "managers/Tables/Table.hpp"
//#include "managers/TableManager.hpp"
#include "managers/Tables/TableTypes.hpp"

/// Function wrapper for 4D tables
class Lookup4DTable : public Function
{

public:
  Lookup4DTable(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  Lookup4DTable(const std::string& name, const Table4D* tablePtr):
    Function(name),
    m_tablePtr(tablePtr)
  {/** Empty **/
  }

  ~Lookup4DTable()
  {/** Empty **/
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    array1d<real64> xx(3);
    xx[0] = (&x)[0];
    xx[1] = (&x)[1];
    xx[2] = (&x)[2];
    xx[3] = (&x)[3];
    return m_tablePtr->Lookup(xx);
  }

  static const char* FunctionName()
  {
    return "Lookup4DTable";
  }

private:
  const Table4D* m_tablePtr;
};

/// Function wrapper for 3D tables
class Lookup3DTable : public Function
{

public:
  Lookup3DTable(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  Lookup3DTable(const std::string& name, const Table3D* tablePtr):
    Function(name),
    m_tablePtr(tablePtr)
  {/** Empty **/
  }

  ~Lookup3DTable()
  {/** Empty **/
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    array1d<real64> xx(3);
    xx[0] = (&x)[0];
    xx[1] = (&x)[1];
    xx[2] = (&x)[2];
    return m_tablePtr->Lookup(xx);
  }

  static const char* FunctionName()
  {
    return "Lookup3DTable";
  }

private:
  const Table3D* m_tablePtr;
};

/// Function wrapper for 2D tables
class Lookup2DTable : public Function
{

public:
  Lookup2DTable(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  Lookup2DTable(const std::string& name, const Table2D* tablePtr):
    Function(name),
    m_tablePtr(tablePtr)
  {/** Empty **/
  }

  ~Lookup2DTable()
  {/** Empty **/
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    array1d<real64> xx(2);
    xx[0] = (&x)[0];
    xx[1] = (&x)[1];
    return m_tablePtr->Lookup(xx);
  }

  static const char* FunctionName()
  {
    return "Lookup2DTable";
  }

private:
  const Table2D* m_tablePtr;
};

/// Function wrapper for 1D tables
class Lookup1DTable : public Function
{

public:
  Lookup1DTable(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  Lookup1DTable(const std::string& name, const Table1D* tablePtr):
    Function(name),
    m_tablePtr(tablePtr)
  {/** Empty **/
  }

  ~Lookup1DTable()
  {/** Empty **/
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    array1d<real64> xx(1);
    xx[0] = x;
    return m_tablePtr->Lookup(xx);
  }

  static const char* FunctionName()
  {
    return "Lookup1DTable";
  }
private:
  const Table1D* m_tablePtr;
};

/// Evaluates a mathematical expression in a string with the warp function
// parser.
/// Should be used for one-time evaluation of real number expressions only.
/// If a user-defined function must be evaluated multiple times use the
/// SymbolicFunction class instead.
realT EvaluateStringFunction(const std::string& fString);

//////////////////////////
//
// Function Factory
//
// Consists of the following parts:
//   * The function to generate new pointers: "newFunction"
//   * A base class to derive the functions to generate Function pointers:
// "FunctionInitializer"
//   * A String-to-Function-Intializer map hidden behind the
// getFunctionCatalogue function
//   * A template to create Function initializers: "FunctionRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_Function"
//
// Most initial conditions will only need to use one or two of the parts:
//   * To register a new Function in the factory: REGISTER_Function(
// FunctionClassName )
//   * To load a Function pointer from the factory:       Function* aFunctionPtr
// = newFunction(FunctionString, args );

/// The Function Factory.
Function* newFunction(const std::string& FunctionName, TICPP::HierarchicalDataNode* hdn,
                      const ProblemManagerT* const pm);

/// Base class to generate new Function pointers
class FunctionInitializer
{
public:
  virtual Function* initializeFunction(TICPP::HierarchicalDataNode* hdn,
                                       const ProblemManagerT* const pm) = 0;
  virtual ~FunctionInitializer() {}
};

/// Interface to the Function name -> Function initializer map
std::map<std::string, FunctionInitializer*> & getFunctionCatalogue();

/// Return a list of supported Function names
void getFunctionNames(std::vector<std::string>& nameList);

/// Template for creating classes derived from FunctionInitializer
template<class FunctionType>
class FunctionRegistrator : public FunctionInitializer
{

public:
  FunctionRegistrator(void)
  {
    std::string FunctionName = std::string(FunctionType::FunctionName());
    getFunctionCatalogue()[FunctionName] = this;
  }

  Function* initializeFunction(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm)
  {
    return new FunctionType(hdn, pm);
  }
};

/// Compiler directive to simplify function autoregistration
#define REGISTER_Function( ClassName ) namespace { FunctionRegistrator<ClassName> reg_ ## ClassName; }
#endif
