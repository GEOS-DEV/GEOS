/*
 * SymbolicFunction.hpp
 *
 *  Created on: July 6, 2017
 *      Author: sherman
 */

#ifndef SYMBOLICFUNCTION_HPP_
#define SYMBOLICFUNCTION_HPP_

#include "FunctionBase.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{


class SymbolicFunction : public FunctionBase
{
public:
  SymbolicFunction( const std::string& name,
                    dataRepository::ManagedGroup * const parent );

  virtual ~SymbolicFunction();
  static string CatalogName() { return "SymbolicFunction"; }
  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;
  virtual void BuildDataStructure( dataRepository::ManagedGroup * const domain ) override;

  virtual void InitializeFunction() override;
  virtual double Evaluate(double* input) override { return parserExpression.evaluate(input); }
  
private:
  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;
};


} /* namespace geosx */

#endif /* SYMBOLICFUNCTION_HPP_ */
