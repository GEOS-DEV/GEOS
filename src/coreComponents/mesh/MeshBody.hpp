/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * MeshBody.hpp
 *
 *  Created on: Sep 13, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHBODY_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHBODY_HPP_

#include "MeshLevel.hpp"


namespace geosx
{



class MeshLevel;

class MeshBody : public dataRepository::ManagedGroup
{
public:
  MeshBody( string const & name,
            ManagedGroup * const parent );
  virtual ~MeshBody();

  MeshLevel * CreateMeshLevel( integer const newLevel );

  MeshLevel * getMeshLevel( integer const level ) { return this->GetGroup<MeshLevel>(level); }
  MeshLevel const * getMeshLevel( integer const level ) const { return this->GetGroup<MeshLevel>(level); }

  void SetTolerance( real64 tolerance )
  {
    m_tolerance = tolerance;
  }

  real64 GetTolerance() const
  {
    return m_tolerance;
  }

  struct viewKeysStruct
  {


    dataRepository::ViewKey meshLevels                = { "meshLevels" };
  } viewKeys;

  struct groupStructKeys
  {} groupKeys;

private:
  /// By default, an absolute tolerance
  /// Can be set to another value
  real64 m_tolerance { std::numeric_limits<real64>::epsilon() };


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_ */
