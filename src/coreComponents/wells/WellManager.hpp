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
 * @file WellManager.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_WELLMANAGER_HPP_
#define GEOSX_CORECOMPONENTS_WELLS_WELLMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto wellManager = "Wells";
}
}

class Well;

/**
 * @class WellManager
 *
 * This class keeps tracks of all the wells
 */  
class WellManager : public dataRepository::ManagedGroup
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  explicit WellManager( string const & name, 
                        dataRepository::ManagedGroup * const parent );

  /**
   * @brief default destructor
   */
  virtual ~WellManager() override;

  /// deleted default constructor
  WellManager() = delete;

  /// deleted copy constructor
  WellManager( WellManager const & ) = delete;

  /// deleted move constructor
  WellManager( WellManager && ) = delete;

  /// deleted assignment operator
  WellManager & operator=( WellManager const & ) = delete;

  /// deleted move operator
  WellManager & operator=( WellManager && ) = delete;


  dataRepository::ManagedGroup * CreateChild( string const & childKey, 
                                              string const & childName ) override;

  /**
   * @brief Getter for a well by name
   * @param name the name of the well
   * @return a pointer to the well
   */
  Well * getWell( string const & name );

  /**
   * @brief Setter for the gravity vector and the gravity flag
   * @param gravity the R1 tensor contained the gravity vector
   * @param gravityFlag a flag to specify whether gravity is taken into account or not
   */
  void setGravityVector( R1Tensor const & gravity, 
                         bool gravityFlag = true );

  /**
   * @brief Getter for the gravity vector
   * @return an R1Tensor containing the gravity vector
   */
  R1Tensor const & getGravityVector() const { return m_gravityVector; }

  /**
   * @brief Getter for the gravity flag
   * @return a boolean containing the gravity flag
   */
  bool getGravityFlag() const { return m_gravityFlag; }

  /**
   * @brief Getter for the material list
   * @return a array of strings containing the list of materials
   */  
  string_array const & getMaterialList() const { return m_materialList; }

  struct viewKeyStruct 
  {

    static constexpr auto materialListString = "materialList";

  } m_ViewKeysWellManager;

  
private:

  // the gravity vector
  R1Tensor m_gravityVector;

  // a flag specifying whether gravity is accounted for or not
  bool m_gravityFlag;

  // the list of materials for the wells
  string_array m_materialList;
};

}

#endif //CORECOMPONENTS_WELLS_WELLMANAGER_HPP_
