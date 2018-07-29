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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * OutputManager.hpp
 *
 *  Created on: Jan 26, 2016
 *      Author: sherman
 */

#ifndef SRC_COMPONENTS_CORE_SRC_OUTPUTMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_OUTPUTMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"


namespace geosx
{
namespace dataRepository
{
namespace keys
{
}
}

class OutputManager : public dataRepository::ManagedGroup
{
public:
  OutputManager( std::string const & name,
                 ManagedGroup * const parent );

  virtual ~OutputManager() override;

  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  struct viewKeyStruct
  {
    dataRepository::ViewKey time = { "time" };
  } viewKeys;

};


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_OUTPUTMANAGER_HPP_ */
