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
 * @file MsrsbLevelBuilder.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBSTRATEGY_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBSTRATEGY_HPP

#include "linearAlgebra/common/PreconditionerBase.hpp"
#include "linearAlgebra/multiscale/mesh/DofManager.hpp"
#include "linearAlgebra/multiscale/mesh/MeshLevel.hpp"
#include "linearAlgebra/multiscale/msrsb/MsrsbLevelBuilderBase.hpp"

namespace geosx
{
namespace multiscale
{

template< typename LAI >
class MsrsbLevelBuilder : public MsrsbLevelBuilderBase< LAI >
{
public:

  /// Alias for base type
  using Base = MsrsbLevelBuilderBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  explicit MsrsbLevelBuilder( string name, LinearSolverParameters params );

  virtual void initializeFineLevel( DomainPartition & domain,
                                    geosx::DofManager const & dofManager,
                                    MPI_Comm const & comm ) override;

  virtual void initializeCoarseLevel( LevelBuilderBase< LAI > & fine,
                                      Matrix const & fineMat ) override;

  multiscale::MeshLevel       & mesh()       { return m_mesh; }
  multiscale::MeshLevel const & mesh() const { return m_mesh; }

  multiscale::MeshObjectManager & manager()
  {
    return m_location == FieldLocation::Node ? m_mesh.nodeManager() : m_mesh.cellManager();
  }

  multiscale::MeshObjectManager const & manager() const
  {
    return m_location == FieldLocation::Node ? m_mesh.nodeManager() : m_mesh.cellManager();
  }

  integer numComp() const { return m_dofManager.numComponents(); }

  virtual bool updateProlongation( Matrix const & fineMatrix ) override;

  virtual std::unique_ptr< PreconditionerBase< LAI > > makeCoarseSolver() const override;

private:

  void createSmoothers();

  void writeProlongationForDebug() const;

  using Base::m_params;
  using Base::m_name;
  using Base::m_prolongation;
  using Base::m_restriction;
  using Base::m_matrix;
  using Base::m_dofManager;
  using Base::m_preSmoother;
  using Base::m_postSmoother;
  using Base::m_fineLevel;

  /// Dof location (cell or node)
  FieldLocation m_location = FieldLocation::Node;

  /// Mesh description at current level
  multiscale::MeshLevel m_mesh;

  /// List of nodes on global boundary
  array1d< globalIndex > m_boundaryDof;

  /// List of nodes that are interior
  array1d< globalIndex > m_interiorDof;

  /// Previous number of smoothing iterations
  integer m_lastNumIter = std::numeric_limits< integer >::max();

  /// Cumulative number of smoothing iterations since last RAP
  integer m_updateLag = 0;
};

} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBSTRATEGY_HPP
