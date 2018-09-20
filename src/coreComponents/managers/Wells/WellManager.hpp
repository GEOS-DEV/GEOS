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

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_WELLMANAGER_HPP_
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_WELLMANAGER_HPP_

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

class WellBase;

class WellManager : public dataRepository::ManagedGroup
{
public:

  WellManager( string const & name, ManagedGroup * const parent );
  virtual ~WellManager() override;

  WellManager() = delete;
  WellManager( WellManager const &) = delete;
  WellManager( WellManager && ) = delete;

  void CreateChild( string const & childKey, string const & childName ) override;

  WellBase * getWell( string const & name );

  void setGravityVector(R1Tensor const & gravity) { m_gravityVector = gravity; }
  R1Tensor const & getGravityVector() const { return m_gravityVector; }

private:

  R1Tensor m_gravityVector;

};

}

#endif //CORECOMPONENTS_MANAGERS_WELLS_WELLMANAGER_HPP_
