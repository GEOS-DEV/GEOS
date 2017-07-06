/*
 * FunctionBase.hpp
 *
 *  Created on: Jun 6, 2017
 *      Author: sherman
 */

#ifndef FUNCTIONBASE_HPP_
#define FUNCTIONBASE_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

class FunctionBase : public dataRepository::ManagedGroup
{
public:
  FunctionBase( const std::string& name,
                dataRepository::ManagedGroup * const parent );

  virtual ~FunctionBase();
  static string CatalogName() { return "FunctionBase"; }
  virtual void InitializeFunction();
  virtual double Evaluate(double* input);

  // Setup catalog
  using CatalogInterface = cxx_utilities::CatalogInterface< FunctionBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  };
};

} /* namespace geosx */

#endif /* FUNCTIONBASE_HPP_ */
