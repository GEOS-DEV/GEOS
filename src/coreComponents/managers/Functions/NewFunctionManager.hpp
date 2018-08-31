/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
  virtual ~NewFunctionManager() override;

  static NewFunctionManager * Instance()
  {
    static NewFunctionManager theFunctionManager("LastFunctionManagerOnEarth", nullptr);

    return &theFunctionManager;
  }

  static string CatalogName() { return "NewFunctionManager"; }
  virtual void FillDocumentationNode() override;
  virtual void CreateChild( string const & functionCatalogKey, string const & functionName ) override;
};


} /* namespace geosx */

#endif /* NEWFUNCTIONMANAGER_HPP_ */
