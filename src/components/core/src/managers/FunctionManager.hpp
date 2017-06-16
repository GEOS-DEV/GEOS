/*
 * FunctionManager.hpp
 *
 *  Created on: Jun 16, 2017
 *      Author: sherman
 */

#ifndef FUNCTIONMANAGER_HPP_
#define FUNCTIONMANAGER_HPP_

#include "PhysicsSolvers/ManagedGroup.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{
namespace dataRepository
{
class ManagedGroup;
}


class JIT_Function : public ManagedGroup
{
public:
  JIT_Function( const std::string& name,
                    ManagedGroup * const parent );

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


class FunctionManager : public ManagedGroup
{
public:
  FunctionManager( const std::string& name,
                   ManagedGroup * const parent );
  virtual ~FunctionManager();
  static string CatalogName() { return "FunctionManager"; }
  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;
  virtual void BuildDataStructure( dataRepository::ManagedGroup * const domain ) override;

  void ReadXML( dataRepository::ManagedGroup& domain, xmlWrapper::xmlNode const & problemNode );
  JIT_Function & CreateFunction( string const & functionCatalogKey, string const & functionName );

};


} /* namespace geosx */

#endif /* FUNCTIONMANAGER_HPP_ */
