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

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_BASISFUNCTIONMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_BASISFUNCTIONMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "fileIO/xmlWrapper.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
}
}


class BasisFunctionManager : public dataRepository::ManagedGroup
{
public:
  BasisFunctionManager() = delete;
  BasisFunctionManager(string const & name, ManagedGroup * const parent);
  virtual ~BasisFunctionManager() override final;

  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

private:
  virtual void ProcessInputFile( xmlWrapper::xmlNode const & targetNode ) override final;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_BASISFUNCTIONMANAGER_HPP_ */
