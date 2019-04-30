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

#ifndef GEOSX_CORECOMPONENTS_WELLS_PERFORATIONMANAGER_HPP
#define GEOSX_CORECOMPONENTS_WELLS_PERFORATIONMANAGER_HPP

#include "dataRepository/ManagedGroup.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto perforations = "Perforations";
}
}

class Perforation;

/**
 * @class PerforationManager
 *
 * This class processes the data from the XML and partitions the perforations
 */  
class PerforationManager : public ObjectManagerBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  explicit PerforationManager( string const & name, 
                               dataRepository::ManagedGroup * const parent );

  /**
   * @brief default destructor
   */
  ~PerforationManager() override;

  /// deleted default constructor
  PerforationManager() = delete;

  /// deleted copy constructor
  PerforationManager( PerforationManager const &) = delete;

  /// deleted move constructor
  PerforationManager( PerforationManager && ) = delete;

  /// deleted assignment operator
  PerforationManager & operator=( PerforationManager const & ) = delete;

  /// deleted move operator
  PerforationManager & operator=( PerforationManager && ) = delete;

  virtual const string getCatalogName() const override;
  
  dataRepository::ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Getter for the total number of perforations
   * @return the number of perforations on this rank
   */
  globalIndex numPerforationsGlobal()  const
  { return integer_conversion<globalIndex>(m_globalPerforationList.size()); }

  /**
   * @brief Getter for perforation iperf
   * @param the global index of the perforation
   * @return a pointer to the Perforation object
   */
  Perforation * getPerforation( globalIndex iperf );

  /**
   * @brief Const getter for perforation iperf
   * @param the global index of the perforation
   * @return a pointer to the const Perforation object
   */
  Perforation const * getPerforation( globalIndex iperf ) const;
  
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
  } viewKeysPerforationManager;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto perforationString = "Perforation";

    dataRepository::GroupKey perforation = { perforationString };

  } groupKeysPerforationManager;

protected:

  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override;

private:

  // global list of perforations for this well
  string_array m_globalPerforationList;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_WELLS_PERFORATIONMANAGER_HPP
