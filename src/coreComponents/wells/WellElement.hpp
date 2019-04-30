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
 * @file WellElement.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_WELLELEMENT_HPP
#define GEOSX_CORECOMPONENTS_WELLS_WELLELEMENT_HPP

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

/**
 * @class WellElement
 *
 * This class describes a well element with its location and the next well element
 */  
class WellElement : public dataRepository::ManagedGroup
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  explicit WellElement( string const & name, 
                        dataRepository::ManagedGroup * const parent );

  /**
   * @brief default destructor
   */
  ~WellElement() override;

  /// deleted default constructor
  WellElement() = delete;

  /// deleted copy constructor
  WellElement( WellElement const &) = delete;

  /// deleted move constructor
  WellElement( WellElement && ) = delete;

  /// deleted assignment operator
  WellElement & operator=( WellElement const & ) = delete;

  /// deleted move operator
  WellElement & operator=( WellElement && ) = delete;

  /**
   * @brief Getter for the physical location of the well element
   * @return an R1Tensor containing the coordinates of the well element
   */
  R1Tensor const & getLocation() const
  { return m_location; }
  
  /**
   * @brief Getter for the name of the well element
   * @return a string containing the name of the well element
   */
  string const & getNextWellElementName() const 
  { return m_nextWellElementName; }

  struct viewKeyStruct
  {
    static constexpr auto locationString = "location";
    static constexpr auto nextWellElementNameString = "nextSegmentName";

    dataRepository::ViewKey location = { locationString };
    dataRepository::ViewKey nextWellElementName = { nextWellElementNameString };    

  } viewKeysWellElement;

private:

  // geometry
  R1Tensor m_location;

  // connectivity
  string m_nextWellElementName;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_WELLELEMENT_HPP
