/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MatrixFreeSolidMechanicsFEMOperator.hpp
 */

// Source includes
#define SELECTED_FE_TYPES H1_Hexahedron_Lagrange1_GaussLegendre2
#include "MatrixFreeSolidMechanicsFEMOperator.hpp"
#include "kernels/SmallStrainResidual.hpp"
#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "mesh/MeshBody.hpp"
#include "linearAlgebra/solvers/CgSolver.hpp"
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"
#include "linearAlgebra/common/LinearOperatorWithBC.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

#include "SolidMechanicsFields.hpp"
#include "LvArray/src/output.hpp"

namespace geos
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;
using namespace constitutive;

MatrixFreeSolidMechanicsFEMOperator::
  MatrixFreeSolidMechanicsFEMOperator( DomainPartition & domain,
                                       map< std::pair< string, string >, array1d< string > > const & meshTargets,
                                       DofManager & dofManager,
                                       string const & finiteElementName,
                                       int const kernelOptimizationOption ):
  m_meshBodies( domain.getMeshBodies() ),
  m_meshTargets( meshTargets ),
  m_dofManager( dofManager ),
  m_finiteElementName( finiteElementName ),
  m_kernelOptimizationOption( kernelOptimizationOption )
{ }

MatrixFreeSolidMechanicsFEMOperator::
  MatrixFreeSolidMechanicsFEMOperator( dataRepository::Group & meshBodies,
                                       map< std::pair< string, string >, array1d< string > > const & meshTargets,
                                       DofManager & dofManager,
                                       string const & finiteElementName,
                                       int const kernelOptimizationOption ):
  m_meshBodies( meshBodies ),
  m_meshTargets( meshTargets ),
  m_dofManager( dofManager ),
  m_finiteElementName( finiteElementName ),
  m_kernelOptimizationOption( kernelOptimizationOption )
{ }

void MatrixFreeSolidMechanicsFEMOperator::apply( ParallelVector const & src, ParallelVector & dst ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const localSrc = src.values();
  arrayView1d< real64 > const localDst = dst.open();
  // We do it by hand to avoid hypre call
  using POLICY = parallelDeviceAsyncPolicy< 1024 >;
  forAll< POLICY >( localDst.size(), [localDst] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localDst[ i ] = 0.0;
  } );

  // {
  // std::cout<<"MatrixFreeSolidMechanicsFEMOperator::apply - bp1"<<std::endl;
  // LvArray::print< parallelDevicePolicy< 32 > >( localSrc );
  // std::cout<<"MatrixFreeSolidMechanicsFEMOperator::apply - bp2"<<std::endl;
  // }



  for( auto const & target: m_meshTargets )
  {
    string const meshBodyName = target.first.first;
    string const meshLevelName = target.first.second;
    arrayView1d< string const > const & regionNames = target.second.toViewConst();
    MeshBody & meshBody = m_meshBodies.getGroup< MeshBody >( meshBodyName );

    MeshLevel * meshLevelPtr = meshBody.getMeshLevels().getGroupPointer< MeshLevel >( meshLevelName );
    if( meshLevelPtr==nullptr )
    {
      meshLevelPtr = meshBody.getMeshLevels().getGroupPointer< MeshLevel >( MeshBody::groupStructKeys::baseDiscretizationString() );
    }
    MeshLevel & mesh = *meshLevelPtr;

    auto const & totalDisplacement = mesh.getNodeManager().getField< fields::solidMechanics::totalDisplacement >();
    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > localSrc2d( totalDisplacement.dimsArray(), 
                                                                           totalDisplacement.stridesArray(),
                                                                           localSrc.dataBuffer() );
    arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > localDst2d( totalDisplacement.dimsArray(), 
                                                                     totalDisplacement.stridesArray(), 
                                                                     localDst.dataBuffer() );

    auto kernelFactory = solidMechanicsLagrangianFEMKernels::SmallStrainResidualFactory( localSrc2d,
                                                                                         localDst2d,
                                                                                         0,
                                                                                         "",
                                                                                         m_kernelOptimizationOption );

    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                    constitutive::ElasticIsotropic,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            m_finiteElementName,
                                                            "solidMaterialNames",
                                                            kernelFactory );
    //parallelDeviceSync();

  }

  // {
  //   std::cout<<"MatrixFreeSolidMechanicsFEMOperator::apply - bp3"<<std::endl;
  //   LvArray::print< parallelDevicePolicy< 32 > >( localDst );
  //   std::cout<<"MatrixFreeSolidMechanicsFEMOperator::apply - bp4"<<std::endl;
  // }

  dst.close();
}

void MatrixFreeSolidMechanicsFEMOperator::computeDiagonal( ParallelVector & GEOS_UNUSED_PARAM( diagonal ) ) const
{
  GEOS_ERROR( "computeDiagonal: operation not yet implemented" );
}

globalIndex MatrixFreeSolidMechanicsFEMOperator::numGlobalRows() const
{
  return m_dofManager.numGlobalDofs();
}

globalIndex MatrixFreeSolidMechanicsFEMOperator::numGlobalCols() const
{
  return m_dofManager.numGlobalDofs();
}

localIndex MatrixFreeSolidMechanicsFEMOperator::numLocalRows() const
{
  return m_dofManager.numLocalDofs();
}

localIndex MatrixFreeSolidMechanicsFEMOperator::numLocalCols() const
{
  return m_dofManager.numLocalDofs();
}

MPI_Comm MatrixFreeSolidMechanicsFEMOperator::comm() const
{
  return MPI_COMM_GEOSX;
}

} /* namespace geos */
#undef SELECTED_FE_TYPES
