/*
 * FunctionManagerJIT.hpp
 *
 *  Created on: Jun 16, 2017
 *      Author: sherman
 */

#ifndef FUNCTIONMANAGERJIT_HPP_
#define FUNCTIONMANAGERJIT_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{


class JIT_Function : public dataRepository::ManagedGroup
{
public:
  JIT_Function( const std::string& name,
                dataRepository::ManagedGroup * const parent );

  virtual ~JIT_Function();
  static string CatalogName() { return "JIT_Function"; }
  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;
  virtual void BuildDataStructure( dataRepository::ManagedGroup * const domain ) override;

  void Compile();
  double Evaluate(double* input) { return parserExpression.evaluate(input); };

  using CatalogInterface = cxx_utilities::CatalogInterface< JIT_Function, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

private:
  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;
};


class FunctionManagerJIT : public dataRepository::ManagedGroup
{
public:
  FunctionManagerJIT( const std::string& name,
                      dataRepository::ManagedGroup * const parent );
  virtual ~FunctionManagerJIT();
  
  static string CatalogName() { return "FunctionManagerJIT"; }
  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;

  using dataRepository::ManagedGroup::ReadXML;
  void ReadXML( dataRepository::ManagedGroup& domain, xmlWrapper::xmlNode const & problemNode );
  JIT_Function & CreateFunction( string const & functionCatalogKey, string const & functionName );

};


} /* namespace geosx */

#endif /* FUNCTIONMANAGERJIT_HPP_ */
