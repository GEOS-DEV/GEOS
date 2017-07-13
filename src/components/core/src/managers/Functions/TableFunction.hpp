/*
 * TableFunction.hpp
 *
 *  Created on: July 6, 2017
 *      Author: sherman
 */

#ifndef TABLEFUNCTION_HPP_
#define TABLEFUNCTION_HPP_

#include "FunctionBase.hpp"

namespace geosx
{


class TableFunction : public FunctionBase
{
public:
  TableFunction( const std::string& name,
                    dataRepository::ManagedGroup * const parent );

  virtual ~TableFunction();
  static string CatalogName() { return "TableFunction"; }
  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;
  virtual void BuildDataStructure( dataRepository::ManagedGroup * const domain ) override;

  virtual void InitializeFunction() override;
  virtual double Evaluate(double* input) override;

private:
  Array1dT<real64_array> m_coordinates;
  real64_array m_values;
  static localIndex constexpr m_maxDimensions = 3;
  localIndex m_dimensions;
  lArray1d m_size;
  lArray1d m_indexIncrement;
  Array1dT<lArray1d> m_corners;
};


} /* namespace geosx */

#endif /* TABLEFUNCTION_HPP_ */
