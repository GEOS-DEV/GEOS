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
 * @file FiniteVolumeManager.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FINITEVOLUMEMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FINITEVOLUMEMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

class DomainPartition;
class FluxApproximationBase;

class FiniteVolumeManager : public dataRepository::ManagedGroup
{
public:

  FiniteVolumeManager() = delete;
  FiniteVolumeManager(string const & name, ManagedGroup * const parent);
  virtual ~FiniteVolumeManager() override;

  void IntermediateInitializationPreSubGroups(ManagedGroup * const rootGroup) override;

  virtual ManagedGroup * CreateChild(string const & childKey, string const & childName) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

  FluxApproximationBase const * getFluxApproximation(string const &name) const;

private:

  void precomputeFiniteVolumeData(DomainPartition * const domain);

};

} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FINITEVOLUMEMANAGER_HPP_
