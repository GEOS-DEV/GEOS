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
 * @file FunctionManager.h
 * @date May 24, 2011
 */

#ifndef FUNCTIONMANAGER_H_
#define FUNCTIONMANAGER_H_

//#include "Common/typedefs.h"

#include <map>
#include "codingUtilities/Functions.hpp"

class FunctionManager
{
public:

  static FunctionManager& Instance()
  {
    static FunctionManager theFunctionManager;

    return theFunctionManager;
  }

  inline std::map< std::string, Function* >& Functions()
  { return m_functions; }

  inline Function& GetFunction( const std::string& functionName ) const
  {
    std::map<std::string,Function* >::const_iterator function = m_functions.find( functionName );
    GEOS_ERROR_IF(function == m_functions.end(), "Error FunctionManager: Function name `" << functionName << "' not found\n")

    return *(function->second);
  }


  inline void AddFunction( const std::string& functionName, Function* aFunction ){
    m_functions.insert(std::make_pair(functionName, aFunction));
  }

  inline realT Eval( const std::string& functionName, const realT& key ) const
  {
    Function& f = GetFunction(functionName);
    return f(key);
  }

private:
  std::map<std::string,Function* > m_functions;

  FunctionManager():
    m_functions()
  {}

  ~FunctionManager() {
    for( std::map<std::string,Function*>::iterator it_func = m_functions.begin() ;
         it_func != m_functions.end() ; ++it_func )
    {
      delete it_func->second;
    }
  }

  FunctionManager( const FunctionManager& );
  FunctionManager& operator=( const FunctionManager& );
};

#endif /* FUNCTIONMANAGER_H_ */
