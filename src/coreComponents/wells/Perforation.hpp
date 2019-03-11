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

class Perforation : public dataRepository::ManagedGroup
{
public:

  explicit Perforation( string const & name, dataRepository::ManagedGroup * const parent );
  ~Perforation() override;

  Perforation() = delete;
  Perforation( Perforation const &) = delete;
  Perforation( Perforation && ) = delete;

  R1Tensor const & getLocation() const
  { return m_location; }

  real64 getTransmissibility() const
  { return m_transmissibility; }

  string const & getWellElementName() const
  { return m_wellElementName; }

  struct viewKeyStruct
  {
    static constexpr auto locationString         = "location";
    static constexpr auto transmissibilityString = "transmissibility";
    static constexpr auto wellElementNameString  = "segmentName";

    dataRepository::ViewKey location         = { locationString         };
    dataRepository::ViewKey transmissibility = { transmissibilityString };
    dataRepository::ViewKey wellElementName  = { wellElementNameString };

  } viewKeysPerforation;

private:
  
  // geometry
  R1Tensor   m_location;
  real64     m_transmissibility;

  // connectivity
  string m_wellElementName;

  // depending on whether we keep the WellStencil class or not,
  // we may need additional member variables here (i.e., cell id)
  
};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATION_HPP
