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
 * MeshBodyManager.hpp
 *
 *  Created on: Oct 4, 2018
 *      Author: Antoine Mazuyer
 */

#pragma once

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

 class MeshBody;

/*!
* @brief Abstract Manager class for the inputs of GEOSX
* @details This class handles meshes and properties that
* can be inputed to GEOSX
*/
class MeshBodyManager : public dataRepository::ManagedGroup
{
public:
  MeshBodyManager( const std::string& name,
                      dataRepository::ManagedGroup * const parent );


  virtual void CreateChild( string const & childKey, string const & childName ) override;
};

} //namespace
