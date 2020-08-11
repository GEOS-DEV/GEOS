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
 * @file Vertice.hpp
 */

#ifndef GEOSX_MESH_VERTICE_HPP_
#define GEOSX_MESH_VERTICE_HPP_

namespace geosx
{

/**
 * @class Vertice
 *
 * An event type for periodic events (using either time or cycle as a basis).
 */
class Vertice
{
public:

    Vertice( const int index);

  /// Destructor
  virtual ~Vertice();

  
private:
  int ind;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_VERTICE_HPP_ */
