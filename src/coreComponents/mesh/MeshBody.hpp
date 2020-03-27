/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MeshBody.hpp
 */

#ifndef GEOSX_MESH_MESHBODY_HPP_
#define GEOSX_MESH_MESHBODY_HPP_

#include "MeshLevel.hpp"


namespace geosx
{



class MeshLevel;

class MeshBody : public dataRepository::Group
{
public:
  MeshBody( string const & name,
            Group * const parent );
  virtual ~MeshBody();

  MeshLevel * CreateMeshLevel( localIndex const newLevel );

  MeshLevel * getMeshLevel( localIndex const level ) { return this->GetGroup< MeshLevel >( level ); }
  MeshLevel const * getMeshLevel( localIndex const level ) const { return this->GetGroup< MeshLevel >( level ); }

  void setGlobalLengthScale( real64 scale );

  real64 getGlobalLengthScale() const
  {
    return m_globalLengthScale;
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
  real64 m_globalLengthScale { 0. };


};

} /* namespace geosx */

#endif /* GEOSX_MESH_MESHBODY_HPP_ */
