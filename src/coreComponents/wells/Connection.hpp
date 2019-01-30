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

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_CONNECTION_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_CONNECTION_HPP

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

  string const & getConnectionName() const    { return m_connectionName; }
  void setConnectionName(string const & name) { m_connectionName = name; }

  struct viewKeyStruct
  {

    static constexpr auto connectionNameString   = "connectionName";
    static constexpr auto nextSegmentIndexString = "nextSegmentIndex";
    static constexpr auto prevSegmentIndexString = "prevSegmentIndex";

    using ViewKey = dataRepository::ViewKey;
    
    ViewKey connectionName   = { connectionNameString };
    ViewKey nextSegmentIndex = { nextSegmentIndexString };
    ViewKey prevSegmentIndex = { prevSegmentIndexString };
    
  } viewKeysConnection;

private:

  string     m_connectionName;
  localIndex m_nextSegmentIndex;
  localIndex m_prevSegmentIndex;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_CONNECTION_HPP
