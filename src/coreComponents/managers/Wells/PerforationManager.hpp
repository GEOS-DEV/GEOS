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
 * @file PerforationManager.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATIONMANAGER_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATIONMANAGER_HPP

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

class PerforationManager : public dataRepository::ManagedGroup
{
public:

  explicit PerforationManager( string const & name, dataRepository::ManagedGroup * const parent );
  ~PerforationManager() override;

  PerforationManager() = delete;
  PerforationManager( PerforationManager const &) = delete;
  PerforationManager( PerforationManager && ) = delete;

  virtual void FillDocumentationNode() override;
  virtual void CreateChild( string const & childKey, string const & childName ) override;

  struct groupKeyStruct
  {
    static constexpr auto perforationString = "Perforation";

  } groupKeysPerfManager;

private:

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATIONMANAGER_HPP
