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
 * @file ProblemRepository.hpp
 */

#ifndef GEOS_DATAREPOSITORY_PROBLEMREPOSITORY_HPP_
#define GEOS_DATAREPOSITORY_PROBLEMREPOSITORY_HPP_

#include "dataRepository/ProblemRepositoryABC.hpp"
#include "dataRepository/Group.hpp"

namespace geos
{
namespace dataRepository
{

/**
 * @brief Base class for the problem data-repository repository, gives access to all the roots Groups.
 *         Usage examples:
 *        - to consult other managers:
 *            * from a FS instance:      ProblemRepository::getManager< FunctionManager >( myFieldSpec )
 *            * form a solver instance:  ProblemRepository::getManager< FieldSpecificationManager >( mySolver )
 *        - to consult data from the ProblemManager (mainInterface / high-level testing):
 *            * get the physics solvers manager:  problemManager.getManager< PhysicsSolverManager >()
 *            * get the CommandLine Group:        problemManager.getManager< CommandLine >()
 *        - to make a mutable problem-unique manager available (DomainPartition here):
 *            // EventManager is available through the ProblemRepository as a mutable problem-unique manager.
 *            template<> inline EventManager & dataRepository::ProblemRepository::getManager()
 *            { return getRootGroup().getGroup< EventManager >( m_gks.domain ); }
 *        - to make a const problem-unique manager available (DomainPartition here):
 *            // EventManager is available through the ProblemRepository as a const problem-unique manager.
 *            template<> inline EventManager const & dataRepository::ProblemRepository::getManager() const
 *            { return getRootGroup().getGroup< EventManager >( m_gks.domain ); }
 *
 *        Manager that are unique per-problem should provides a getManager() template specialization in its header.
 */
class ProblemRepository : public ProblemRepositoryABC
{
public:

  /**
   * @brief Gives the problem data repository interface from the given Group
   * @param group The current Group in the Problem tree
   * @return A reference to the problem data repository interface
   */
  static ProblemRepository & get( Group & group );

  /**
   * @copydoc get( Group & )
   */
  static ProblemRepository const & get( Group const & group );

  /**
   * @brief Get a root data-repository object (a manager) of a given problem. For each manager:
   *        - The consumer(s) need to include the type definition,
   *        - The manager implementation needs to implement a specialization.
   * @tparam ManagerType the type of the root data-repository object we want to get (DomainPartition, EventManager...)
   * @return ManagerType& the root data-repository object instance reference
   */
  template< typename ManagerType >
  ManagerType & getManager();

  /**
   * @copydoc getManager()
   */
  template< typename ManagerType >
  ManagerType const & getManager() const;

  /**
   * @copydoc getManager()
   */
  template< typename ManagerType >
  static ManagerType & getManager( Group & group )
  { return ProblemRepository::get( group ).getManager< ManagerType >(); }

  /**
   * @copydoc getManager()
   */
  template< typename ManagerType >
  static ManagerType const & getManager( Group const & group )
  { return ProblemRepository::get( group ).getManager< ManagerType >(); }

  /**
   * @return the root Group which contain all the problem data-repository
   */
  Group & getRootGroup()
  { return m_rootGroup; }

  /**
   * @return the root Group which contain all the problem data-repository
   */
  Group const & getRootGroup() const
  { return m_rootGroup; }

  // if an abstract Group getting method is absolutely needed, we can add:
  //
  //     Group & getManager( string_view managerKey ) = 0;
  //     Group const & getManager( string_view managerKey ) const = 0;
  //
  // ... but ideally, we don't want to propose these to remove any "invisible" circular dependency practice.

protected:

  /**
   * @brief Standard GEOS data-managers Group keys for efficient getManager() lookup.
   *        Note that this list is not a constraint, it can be extended with any type by adding a new
   *        getManager() specialization.
   *        This struct remains internal to avoid "invisible dependancies": being dependent of the
   *        existence of a Group through a specific data-repository path without of mention its type,
   *        complicating dependencies visibility.
   *        When adding a new key, const, and optionnally mutable specializations of getManager must be
   *        added in the Group type header.
   */
  struct ProblemGroupKeys : ProblemRepositoryABC::ProblemGroupKeys
  {
    dataRepository::GroupKey commandLine               = { "commandLine" };
    dataRepository::GroupKey domain                    = { "domain" };
    dataRepository::GroupKey constitutiveManager       = { "Constitutive" };
    dataRepository::GroupKey eventManager              = { "Events" };
    dataRepository::GroupKey externalDataSourceManager = { "ExternalDataSource" };
    dataRepository::GroupKey fieldSpecificationManager = { "FieldSpecifications" };
    dataRepository::GroupKey functionManager           = { "Functions" };
    dataRepository::GroupKey geometricObjectManager    = { "Geometry" };
    dataRepository::GroupKey meshManager               = { "Mesh" };
    dataRepository::GroupKey numericalMethodsManager   = { "NumericalMethods" };
    dataRepository::GroupKey outputManager             = { "Outputs" };
    dataRepository::GroupKey physicsSolverManager      = { "Solvers" };
    dataRepository::GroupKey tasksManager              = { "Tasks" };
  } m_gks;

  /**
   * @brief Construct a new ProblemRepository object
   * @param rootGroup The Problem root group reference.
   */
  ProblemRepository( Group & rootGroup )
    : ProblemRepositoryABC()
    , m_rootGroup( rootGroup )
  {}

  /**
   * @brief Deleted copy constructor for m_rootGroup reference validity.
   */
  ProblemRepository( ProblemRepository const & ) = delete;

  /**
   * @brief Deleted move constructor for m_rootGroup reference validity.
   */
  ProblemRepository( ProblemRepository && ) = delete;

  /**
   * @brief Deleted move operator for m_rootGroup reference validity.
   */
  ProblemRepository & operator=( ProblemRepository const & ) = delete;

  /**
   * @brief Deleted move operator for m_rootGroup reference validity.
   */
  ProblemRepository & operator=( ProblemRepository && ) = delete;

private:

  /**
   * @brief The Problem root group.
   */
  Group & m_rootGroup;

};

} /* namespace dataRepository */
} /* namespace geos */

#endif /* GEOS_DATAREPOSITORY_PROBLEMREPOSITORY_HPP_ */
