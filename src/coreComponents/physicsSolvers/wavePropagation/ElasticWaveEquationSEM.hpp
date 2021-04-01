/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file ElasticWaveEquationSEM.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEM_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEM_HPP_

#include "mesh/ExtrinsicMeshData.hpp"
#include "physicsSolvers/SolverBase.hpp"


namespace geosx
{

class ElasticWaveEquationSEM : public SolverBase
{
public:
  ElasticWaveEquationSEM( const std::string & name,
                           Group * const parent );

  virtual ~ElasticWaveEquationSEM() override;

  ElasticWaveEquationSEM() = delete;
  ElasticWaveEquationSEM( ElasticWaveEquationSEM const & ) = delete;
  ElasticWaveEquationSEM( ElasticWaveEquationSEM && ) = default;

  ElasticWaveEquationSEM & operator=( ElasticWaveEquationSEM const & ) = delete;
  ElasticWaveEquationSEM & operator=( ElasticWaveEquationSEM && ) = delete;


  static string catalogName() { return "ElasticSEM"; }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override final;


  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual
  real64 solverStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain ) override;

  virtual
  real64 explicitStep( real64 const & time_n,
                       real64 const & dt,
                       integer const cycleNumber,
                       DomainPartition & domain ) override;

  /// Returns the value of a Ricker at time t0 with central frequency f0
  virtual
  real64 evaluateRicker( real64 const & t0, real64 const & f0 );

  /// Returns the value of the first derivative of a Ricker at time t0 with central frequency f0
  virtual
  real64 evaluateRickerOrder1( real64 const & t, real64 const & f0 );

  /// Returns the value of the second derivative of a Ricker at time t0 with central frequency f0
  virtual
  real64 evaluateRickerOrder2( real64 const & t, real64 const & f0 );

  /**@}*/


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr char const * sourceCoordinatesString() { return "sourceCoordinates"; }
    static constexpr char const * sourceNodeIdsString() { return "sourceNodeIds"; }
    static constexpr char const * sourceConstantsString() { return "sourceConstants"; }
    static constexpr char const * sourceIsLocalString() { return "sourceIsLocal"; }

    static constexpr char const * timeSourceFrequencyString() { return "timeSourceFrequency"; }

    static constexpr char const * receiverCoordinatesString() { return "receiverCoordinates"; }
    static constexpr char const * receiverNodeIdsString() { return "receiverNodeIds"; }
    static constexpr char const * receiverConstantsString() {return "receiverConstants"; }
    static constexpr char const * receiverIsLocalString() { return "receiverIsLocal"; }

    static constexpr char const * rickerOrderString() { return "rickerOrder"; }
    
    static constexpr char const * displacementNp1AtReceiversString() { return "displacementNp1AtReceivers"; }

    static constexpr char const * outputSismoTraceString() { return "outputSismoTrace";}

  } waveEquationViewKeys;


protected:

  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  /// Locates the source term and precomputes the constant part of the source term
  /// And locate receivers and pre_evaluate the basis functions at each receiver coordinate
  void precomputeSourceAndReceiverTerm( MeshLevel & mesh );

  /// Multiply the precomputed term by the ricker and add to the right-hand side
  void addSourceToRightHandSide( real64 const & time, arrayView2d< real64 > const rhs );

  /// Apply free surface condition to the face define in the geometry box from the xml
  void applyFreeSurfaceBC( real64 const time, DomainPartition & domain );

  /// Compute the pressure at each receiver coordinate in one time step
  void computeSismoTrace( localIndex const num_timestep, arrayView2d< real64 > const pressure_np1 );

  /// save the sismo trace in file
  void saveSismo( localIndex isismo, real64 val_pressure, char *filename );

  /// Coordinates of the sources in the mesh
  array2d< real64 > m_sourceCoordinates;

  /// Indices of the nodes (in the right order) for each source point
  array2d< localIndex > m_sourceNodeIds;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds
  array2d< real64 > m_sourceConstants;

  /// Flag that indicates whether the source is local or not to the MPI rank
  array1d< localIndex > m_sourceIsLocal;

  /// Central frequency for the Ricker time source
  real64 m_timeSourceFrequency;

  /// Coordinates of the receivers in the mesh
  array2d< real64 > m_receiverCoordinates;

  /// Indices of the element nodes (in the right order) for each receiver point
  array2d< localIndex > m_receiverNodeIds;

  /// Basis function evaluated at the receiver for the nodes listed in m_receiverNodeIds
  array2d< real64 > m_receiverConstants;

  /// Flag that indicates whether the receiver is local or not to the MPI rank
  array1d< localIndex > m_receiverIsLocal;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real64 > m_displacementNp1AtReceivers;

    /// Flag that indicates the order of the Ricker to be used, order 2 by default
  localIndex m_rickerOrder;

  /// Flag that indicates if we write the sismo trace in a file .txt, 0 no output, 1 otherwise
  localIndex m_outputSismoTrace;

};

namespace extrinsicMeshData
{

EXTRINSIC_MESH_DATA_TRAIT( Displacement_nm1,
                           "displacement_nm1",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Vectorial displacement at time n-1." );

EXTRINSIC_MESH_DATA_TRAIT( Displacement_n,
                           "displacement_n",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Vectorial displacement at time n." );

EXTRINSIC_MESH_DATA_TRAIT( Displacement_np1,
                           "displacement_np1",
                           array2d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Vectorial displacement at time n+1." );

EXTRINSIC_MESH_DATA_TRAIT( ForcingRHS,
                           "rhs",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "RHS" );

EXTRINSIC_MESH_DATA_TRAIT( MassVector,
                           "massVector",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Diagonal Mass Matrix." );

EXTRINSIC_MESH_DATA_TRAIT( DampingVector,
                           "dampingVector",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Diagonal Damping Matrix." );

EXTRINSIC_MESH_DATA_TRAIT( MediumVelocityVp,
                           "mediumVelocityVp",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "P-waves speed in the cell" );

EXTRINSIC_MESH_DATA_TRAIT( MediumVelocityVs,
                           "mediumVelocityVs",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "S-waves speed in the cell" );

EXTRINSIC_MESH_DATA_TRAIT( MediumDensity,
                           "mediumDensity",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Medium density of the cell" );

EXTRINSIC_MESH_DATA_TRAIT( StiffnessVector,
                           "stiffnessVector",
                           array2d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Stiffness vector contains R_h*Pressure_n." );

EXTRINSIC_MESH_DATA_TRAIT( FreeSurfaceFaceIndicator,
                           "freeSurfaceFaceIndicator",
                           array1d< localIndex >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

EXTRINSIC_MESH_DATA_TRAIT( FreeSurfaceNodeIndicator,
                           "freeSurfaceNodeIndicator",
                           array1d< localIndex >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Free surface indicator, 1 if a node is on free surface 0 otherwise." );


}


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASSTICWAVEEQUATIONSEM_HPP_ */
