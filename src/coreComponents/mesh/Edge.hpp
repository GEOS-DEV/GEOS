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
 * @file Edge.hpp
 */

#ifndef GEOSX_MESH_EDGE_HPP_
#define GEOSX_MESH_EDGE_HPP_

#include "common/Path.hpp"
#include "mesh/GraphBase.hpp"

namespace geosx
{

/**
 * @class Edge
 *
 * An event type for periodic events (using either time or cycle as a basis).
 */
class Edge : public GraphBase
{
public:

  /// @copydoc geosx::dataRepository::Group::Group( std::string const & name, Group * const parent )
  Edge( const std::string & name,
                 Group * const parent );

  /// Destructor
  virtual ~Edge() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string CatalogName() { return "Edge"; }

  virtual void GenerateGraph() override;

 
  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr auto fileString = "file";

     } viewKeys;
  /// @endcond

private:
  Path m_file;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_EDGE_HPP_ */
