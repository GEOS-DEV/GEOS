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
 * @file Perforation.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_PERFORATION_HPP
#define GEOSX_CORECOMPONENTS_WELLS_PERFORATION_HPP

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

/**
 * @class Perforation
 *
 * This class describes a perforation with its location, transmissibility and corresponding well element
 */  
class Perforation : public dataRepository::ManagedGroup
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  explicit Perforation( string const & name, dataRepository::ManagedGroup * const parent );
  
  /**
   * @brief default destructor
   */
  ~Perforation() override;

  /// deleted default constructor
  Perforation() = delete;

  /// deleted copy constructor
  Perforation( Perforation const &) = delete;

  /// deleted move constructor
  Perforation( Perforation && ) = delete;

  /// deleted assignment operator
  Perforation & operator=( Perforation const & ) = delete;

  /// deleted move operator
  Perforation & operator=( Perforation && ) = delete;

  /**
   * @brief Getter for the physical location of the perforation
   * @return an R1Tensor containing the coordinates of the perforation
   */
  R1Tensor const & getLocation() const { return m_location; }

  /**
   * @brief Getter for the transmissibility at the perforation
   * @return the transmissibility
   */
  real64 getTransmissibility() const { return m_transmissibility; }

  /**
   * @brief Getter for the name of the well element this perforation is attached to
   * @return a string for the name of the well element
   */
  string const & getWellElementName() const { return m_wellElementName; }

  struct viewKeyStruct
  {
    static constexpr auto locationString         = "location";
    static constexpr auto transmissibilityString = "transmissibility";
    static constexpr auto wellElementNameString  = "segmentName";

    dataRepository::ViewKey location         = { locationString         };
    dataRepository::ViewKey transmissibility = { transmissibilityString };
    dataRepository::ViewKey wellElementName  = { wellElementNameString };

  } viewKeysPerforation;

protected: 

  virtual void PostProcessInput() override;

private:
  
  // geometry
  R1Tensor   m_location;
  real64     m_transmissibility;

  // connectivity
  string m_wellElementName;
  
};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATION_HPP
