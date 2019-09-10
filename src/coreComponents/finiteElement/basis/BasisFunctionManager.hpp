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

#include "dataRepository/Group.hpp"
#include "dataRepository/xmlWrapper.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
}
}


class BasisFunctionManager : public dataRepository::Group
{
public:
  BasisFunctionManager() = delete;
  BasisFunctionManager(string const & name, Group * const parent);
  virtual ~BasisFunctionManager() override final;

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;
  
  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

private:
  // virtual void ProcessInputFile( xmlWrapper::xmlNode const & targetNode ) override final;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_BASISFUNCTIONMANAGER_HPP_ */
