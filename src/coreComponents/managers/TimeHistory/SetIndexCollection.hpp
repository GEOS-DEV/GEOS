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
 * @file SetIndexCollection.hpp
 */

#ifndef GEOSX_SetIndexCollection_HPP_
#define GEOSX_SetIndexCollection_HPP_

#include "TimeHistoryCollection.hpp"

namespace geosx
{

/**
 * @class SetIndexCollection
 *
 * A task class for serializing set collection-index information (indices into HistoryOutput related to a Set,
 *   not the simulation indices of the Set).
 */
class SetIndexCollection : public HistoryCollection
{
public:
/**
 * @brief Constructor
 * @param object_path The data repository path of the object to collect set index information about
 * @param set_name The name of the set to collect indices from
 * @param set_index_offset An offset constituting all 'prior' rank-sets + local-set index counts to correctly provide history output set
 * indices.
 */
  SetIndexCollection( string const & object_path, string const & set_name, globalIndex set_index_offset );

/// @copydoc geosx::HistoryCollection::GetMetadata
  virtual HistoryMetadata GetMetadata( Group * problem_group ) override;

/// @copydoc geosx::HistoryCollection::Collect
  virtual void Collect( Group * domain,
                        real64 const time_n,
                        real64 const dt,
                        buffer_unit_type * & buffer ) override;

protected:
  /// The path of the object to collect set indices from
  string m_object_path;
  /// The name of the set to collect indices about
  string m_set_name;
  /// The offset of this MPI rank for just this index set's indices (assuming densely ordered indices)
  globalIndex m_set_index_offset;
};
}
#endif
