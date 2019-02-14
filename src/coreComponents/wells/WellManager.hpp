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
  
class WellManager : public dataRepository::ManagedGroup
{
public:

  explicit WellManager( string const & name, dataRepository::ManagedGroup * const parent );
  virtual ~WellManager() override;

  WellManager() = delete;
  WellManager( WellManager const &) = delete;
  WellManager( WellManager && ) = delete;

  dataRepository::ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;
  
  Well * getWell( string const & name );

  void setGravityVector(R1Tensor const & gravity, bool gravityFlag = true);

  R1Tensor const & getGravityVector() const { return m_gravityVector; }
  bool getGravityFlag() const { return m_gravityFlag; }
  
private:

  R1Tensor m_gravityVector;
  bool m_gravityFlag;

};

}

#endif //CORECOMPONENTS_WELLS_WELLMANAGER_HPP_
