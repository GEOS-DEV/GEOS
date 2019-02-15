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
 * @file Connection.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_CONNECTION_HPP
#define GEOSX_CORECOMPONENTS_WELLS_CONNECTION_HPP

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

class Connection : public dataRepository::ManagedGroup
{
public:

  explicit Connection( string const & name, dataRepository::ManagedGroup * const parent );
  ~Connection() override;

  Connection() = delete;
  Connection( Connection const &) = delete;
  Connection( Connection && ) = delete;

  string const & getNextWellElementName() const
  { return m_nextWellElementName; }
  
  string const & getPreviousWellElementName() const
  { return m_prevWellElementName; }  
  
  // check if the connection is an exit
  bool isExitConnection() const
  { return (m_nextWellElementName.empty() || m_prevWellElementName.empty() ); }
  
  struct viewKeyStruct
  {
    static constexpr auto nextWellElementNameString = "nextSegmentName";
    static constexpr auto prevWellElementNameString = "prevSegmentName";

    using ViewKey = dataRepository::ViewKey;

    ViewKey nextWellElementName = { nextWellElementNameString };
    ViewKey prevWellElementName = { prevWellElementNameString };
    
  } viewKeysConnection;

private:

  // connectivity info
  string m_nextWellElementName;
  string m_prevWellElementName;
  
};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_WELLS_CONNECTION_HPP
