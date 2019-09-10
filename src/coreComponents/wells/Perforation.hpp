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

#include "dataRepository/Group.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto perforation = "Perforation";
}
}


/**
 * @class Perforation
 *
 * This class describes a perforation with its location, transmissibility and corresponding well element
 */  
class Perforation : public dataRepository::Group
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  explicit Perforation( string const & name, dataRepository::Group * const parent );
  
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
   * @brief Getter for the linear distance between the well head and the perforation
   * @return the distance between the well head and the perforation
   */
  real64 const & GetDistanceFromWellHead() const { return m_distanceFromHead; }

  /**
   * @brief Getter for the transmissibility at the perforation
   * @return the transmissibility
   */
  real64 GetTransmissibility() const { return m_transmissibility; }


  struct viewKeyStruct
  {
    static constexpr auto distanceFromHeadString = "distanceFromHead";
    static constexpr auto transmissibilityString = "transmissibility";

    dataRepository::ViewKey distanceFromHead = { distanceFromHeadString };
    dataRepository::ViewKey transmissibility = { transmissibilityString };

  } viewKeysPerforation;

protected:

  void PostProcessInput() override;

private:
  
  // linear distance from well head
  real64 m_distanceFromHead;

  // transmissibility (well index) at this perforation
  real64 m_transmissibility;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATION_HPP
