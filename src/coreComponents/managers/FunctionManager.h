/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file FunctionManager.h
 * @author walsh24
 * @date May 24, 2011
 */

#ifndef FUNCTIONMANAGER_H_
#define FUNCTIONMANAGER_H_

//#include "Common/typedefs.h"

#include <map>
#include "codingUtilities/Functions.hpp"

#ifdef GEOSX_USE_ATK
#include <slic/slic.hpp>
#endif

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
    if( function == m_functions.end() )
    {
#ifdef GEOSX_USE_ATK
      SLIC_ERROR("Error FunctionManager: Function name `" + functionName + "' not found\n");
#endif
    }

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
