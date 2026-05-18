/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ProblemManagerBase.hpp
 */

#ifndef GEOS_DATAREPOSITORY_PROBLEMMANAGERBASE_HPP_
#define GEOS_DATAREPOSITORY_PROBLEMMANAGERBASE_HPP_

#include "dataRepository/Group.hpp"

namespace geos
{

class DomainPartition;
class EventManager;
class ExternalDataSourceManager;
class FieldSpecificationManager;
class FunctionManager;
class GeometricObjectManager;
class MeshManager;
class NumericalMethodsManager;
class OutputManager;
class PhysicsSolverManager;
class TasksManager;
namespace constitutive
{
class ConstitutiveManager;
}

namespace dataRepository
{


/**
 * @class ProblemManagerBase
 */
class ProblemManagerBase : public Group
{
public:

  using Group::Group;

  virtual DomainPartition & getDomainPartition() = 0;
  virtual DomainPartition const & getDomainPartition() const = 0;

  virtual FunctionManager & getFunctionManager() = 0;
  virtual FunctionManager const & getFunctionManager() const = 0;

  virtual FieldSpecificationManager & getFieldSpecificationManager() = 0;
  virtual FieldSpecificationManager const & getFieldSpecificationManager() const = 0;

  virtual EventManager & getEventManager() = 0;
  virtual TasksManager & getTasksManager() = 0;
  virtual PhysicsSolverManager & getPhysicsSolverManager() = 0;
  virtual NumericalMethodsManager & getNumericalMethodsManager() = 0;
  virtual MeshManager & getMeshManager() = 0;
  virtual OutputManager & getOutputManager() = 0;
  virtual GeometricObjectManager & getGeometricObjectManager() = 0;
  virtual constitutive::ConstitutiveManager & getConstitutiveManager() = 0;

  virtual string const & getProblemName() const = 0;
  virtual string const & getInputFileName() const = 0;
  virtual string const & getRestartFileName() const = 0;
  virtual string const & getSchemaFileName() const = 0;

};

} /* namespace dataRepository */
} /* namespace geos */


#endif /* GEOS_DATAREPOSITORY_PROBLEMMANAGERBASE_HPP_ */
