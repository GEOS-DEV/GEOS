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

class WellElement : public dataRepository::ManagedGroup
{
public:

  explicit WellElement( string const & name, dataRepository::ManagedGroup * const parent );
  ~WellElement() override;

  WellElement() = delete;
  WellElement( WellElement const &) = delete;
  WellElement( WellElement && ) = delete;

  R1Tensor const & getLocation() const
  { return m_location; }
  
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
