/*
 * FunctionBase.hpp
 *
 *  Created on: July 6, 2017
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

  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;
  
  static string CatalogName() { return "FunctionBase"; }
  virtual void InitializeFunction(){}
  virtual real64 Evaluate( real64 const * const input ) const = 0;

  // Setup catalog
  using CatalogInterface = cxx_utilities::CatalogInterface< FunctionBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }
};

} /* namespace geosx */

#endif /* FUNCTIONBASE_HPP_ */
