/*
 * CompositeFunction.hpp
 *
 *  Created on: August 18, 2017
 *      Author: sherman
 */

#ifndef COMPOSITEFUNCTION_HPP_
#define COMPOSITEFUNCTION_HPP_

#include "FunctionBase.hpp"
#include <mathpresso/mathpresso.h>

namespace geosx
{


class CompositeFunction : public FunctionBase
{
public:
  CompositeFunction( const std::string& name,
                    dataRepository::ManagedGroup * const parent );

  virtual ~CompositeFunction();
  static string CatalogName() { return "CompositeFunction"; }
  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;
  virtual void BuildDataStructure( dataRepository::ManagedGroup * const domain ) override;

  virtual void InitializeFunction() override;

  virtual void Evaluate( dataRepository::ManagedGroup const * const group,
                         real64 const time,
                         lSet const & sets,
                         real64_array & result ) const override final;

  virtual real64 Evaluate( real64 const * const input) const override final;
  
private:
  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;

  localIndex m_numSubFunctions;
  std::vector<FunctionBase*> m_subFunctions;

};


} /* namespace geosx */

#endif /* COMPOSITEFUNCTION_HPP_ */
