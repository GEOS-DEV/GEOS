/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FunctionManager.hpp
 */

#ifndef GEOSX_MANAGERS_FUNCTIONS_FUNCTIONMANAGER_HPP_
#define GEOSX_MANAGERS_FUNCTIONS_FUNCTIONMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "FunctionBase.hpp"

namespace geosx
{



class FunctionManager : public dataRepository::Group
{
public:
  FunctionManager( const std::string& name,
                      dataRepository::Group * const parent );
  virtual ~FunctionManager() override;

  static FunctionManager & Instance()
  {
    static FunctionManager theFunctionManager( "Functions", nullptr ); 
    return theFunctionManager;
  }

  static string CatalogName() { return "FunctionManager"; }
  virtual Group * CreateChild( string const & functionCatalogKey, string const & functionName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;
  
};


} /* namespace geosx */

#endif /* GEOSX_MANAGERS_FUNCTIONS_FUNCTIONMANAGER_HPP_ */
