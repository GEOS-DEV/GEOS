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
 * @file WellElementManager.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_WELLELEMENTMANAGER_HPP
#define GEOSX_CORECOMPONENTS_WELLS_WELLELEMENTMANAGER_HPP

#include "dataRepository/ManagedGroup.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto wellElements = "Segments";
}
}

class WellElement;

/**
 * @class WellElementManager
 *
 * This class processes the data from the XML and partitions the wellElements
 */  
class WellElementManager : public ObjectManagerBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  explicit WellElementManager( string const & name, dataRepository::ManagedGroup * const parent );

  /**
   * @brief default destructor
   */
  ~WellElementManager() override;

  /// deleted default constructor
  WellElementManager() = delete;

  /// deleted copy constructor
  WellElementManager( WellElementManager const &) = delete;

  /// deleted move constructor
  WellElementManager( WellElementManager && ) = delete;

  /// deleted assignment operator
  WellElementManager & operator=( WellElementManager const & ) = delete;

  /// deleted move operator
  WellElementManager & operator=( WellElementManager && ) = delete;

  virtual const string getCatalogName() const override;
  
  dataRepository::ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Getter for the total number of well elements
   * @return the number of well elements on this rank
   */
  globalIndex numWellElementsGlobal()  const
  { return integer_conversion<globalIndex>(m_globalWellElementList.size()); }

  /**
   * @brief Getter for well element iwelem
   * @param the global index of the well element
   * @return a pointer to the WellElement object
   */
  WellElement * getWellElement( globalIndex iwelem );

  /**
   * @brief Const Getter for well element iwelem
   * @param the global index of the well element
   * @return a pointer to the const WellElement object
   */
  WellElement const * getWellElement( globalIndex iwelem ) const;


  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
  } viewKeysWellElementManager;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto wellElementString = "Segment";

    dataRepository::GroupKey wellElement = { wellElementString };

  } groupKeysWellElementManager;

protected:

  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override;

private:

  // global list of well elements for this well
  string_array m_globalWellElementList;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_WELLS_WELLELEMENTMANAGER_HPP
