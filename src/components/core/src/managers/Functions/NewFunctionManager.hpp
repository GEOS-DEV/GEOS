/*
 * NewFunctionManager.hpp
 *
 *  Created on: July 6, 2017
 *      Author: sherman
 */

#ifndef NEWFUNCTIONMANAGER_HPP_
#define NEWFUNCTIONMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "FunctionBase.hpp"

namespace geosx
{



class NewFunctionManager : public dataRepository::ManagedGroup
{
public:
  NewFunctionManager( const std::string& name,
                      dataRepository::ManagedGroup * const parent );
  virtual ~NewFunctionManager();
  
  static NewFunctionManager * Instance()
  {
      static NewFunctionManager theFunctionManager("LastFunctionManagerOnEarth", nullptr);

    return &theFunctionManager;
  }


  static string CatalogName() { return "NewFunctionManager"; }
  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;

  using dataRepository::ManagedGroup::ReadXML;
  void ReadXML( dataRepository::ManagedGroup * domain, xmlWrapper::xmlNode const & problemNode );

  FunctionBase * CreateFunction( string const & functionCatalogKey, string const & functionName );
};


} /* namespace geosx */

#endif /* NEWFUNCTIONMANAGER_HPP_ */
