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
 * @file ElasticFirstOrderWaveEquationSEM.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICFIRSTORDERWAVEEQUATIONSEM_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICFIRSTORDERWAVEEQUATIONSEM_HPP_

#include "mesh/MeshFields.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "WaveSolverBase.hpp"



namespace geosx
{

class ElasticFirstOrderWaveEquationSEM : public WaveSolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy< 32 >;
  using ATOMIC_POLICY = parallelDeviceAtomic;

  static constexpr real64 epsilonLoc = 1e-8;

  ElasticFirstOrderWaveEquationSEM( const std::string & name,
                                    Group * const parent );

  virtual ~ElasticFirstOrderWaveEquationSEM() override;

  ElasticFirstOrderWaveEquationSEM() = delete;
  ElasticFirstOrderWaveEquationSEM( ElasticFirstOrderWaveEquationSEM const & ) = delete;
  ElasticFirstOrderWaveEquationSEM( ElasticFirstOrderWaveEquationSEM && ) = default;

  ElasticFirstOrderWaveEquationSEM & operator=( ElasticFirstOrderWaveEquationSEM const & ) = delete;
  ElasticFirstOrderWaveEquationSEM & operator=( ElasticFirstOrderWaveEquationSEM && ) = delete;


  static string catalogName() { return "ElasticFirstOrderSEM"; }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override final;


  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual real64 explicitStepForward( real64 const & time_n,
                                      real64 const & dt,
                                      integer const cycleNumber,
                                      DomainPartition & domain,
                                      bool const computeGradient ) override;

  virtual real64 explicitStepBackward( real64 const & time_n,
                                       real64 const & dt,
                                       integer const cycleNumber,
                                       DomainPartition & domain,
                                       bool const computeGradient ) override;

  /**@}*/

  /**
   * @brief Multiply the precomputed term by the Ricker and add to the right-hand side
   * @param cycleNumber the cycle number/step number of evaluation of the source
   * @param rhs the right hand side vector to be computed
   */
  void addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs );

  /**
   * TODO: move implementation into WaveSolverBase
   * @brief Compute the sesimic traces for a given variable at each receiver coordinate at a given time, using the field values at the
   * last two timesteps.
   * @param time_n the time corresponding to the field values pressure_n
   * @param dt the simulation timestep
   * @param timeSeismo the time at which the seismogram is computed
   * @param iSeismo the index of the seismogram time in the seismogram array
   * @param var_np1 the field values at time_n + dt
   * @param var_n the field values at time_n
   * @param varAtreceivers the array holding the trace values, where the output is written
   */
  virtual void computeSeismoTrace( real64 const time_n,
                                   real64 const dt,
                                   real64 const timeSeismo,
                                   localIndex const iSeismo,
                                   arrayView1d< real32 const > const var_np1,
                                   arrayView1d< real32 const > const var_n,
                                   arrayView2d< real32 > varAtReceivers ) override;

  /**
   * TODO: move implementation into WaveSolverBase
   * @brief Computes the traces on all receivers (see @computeSeismoTraces) up to time_n+dt
   * @param time_n the time corresponding to the field values pressure_n
   * @param dt the simulation timestep
   * @param var_np1 the field values at time_n + dt
   * @param var_n the field values at time_n
   * @param varAtreceivers the array holding the trace values, where the output is written
   */
  virtual void computeAllSeismoTraces( real64 const time_n,
                                       real64 const dt,
                                       arrayView1d< real32 const > const var_np1,
                                       arrayView1d< real32 const > const var_n,
                                       arrayView2d< real32 > varAtReceivers );

  /**
   * @brief Initialize Perfectly Matched Layer (PML) information
   */
  virtual void initializePML() override;


  /**
   * @brief Overridden from ExecutableGroup. Used to write last seismogram if needed.
   */
  virtual void cleanup( real64 const time_n, integer const cycleNumber, integer const eventCounter, real64 const eventProgress, DomainPartition & domain ) override;


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr char const * sourceNodeIdsString() { return "sourceNodeIds"; }
    static constexpr char const * sourceConstantsString() { return "sourceConstants"; }
    static constexpr char const * sourceIsAccessibleString() { return "sourceIsAccessible"; }

    static constexpr char const * receiverNodeIdsString() { return "receiverNodeIds"; }
    static constexpr char const * receiverConstantsString() {return "receiverConstants"; }
    static constexpr char const * receiverIsLocalString() { return "receiverIsLocal"; }

    static constexpr char const * displacementxNp1AtReceiversString() { return "displacementxNp1AtReceivers"; }
    static constexpr char const * displacementyNp1AtReceiversString() { return "displacementyNp1AtReceivers"; }
    static constexpr char const * displacementzNp1AtReceiversString() { return "displacementzNp1AtReceivers"; }

  } waveEquationViewKeys;

  /** internal function to the class to compute explicitStep either for backward or forward.
   * (requires not to be private because it is called from GEOSX_HOST_DEVICE method)
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @return return the timestep that was achieved during the step.
   */
  real64 explicitStepInternal( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain );

protected:

  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  /**
   * @brief Locate sources and receivers position in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param mesh mesh of the computational domain
   * @param regionNames name of the region you are currently on
   */
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) override;

  /**
   * @brief Apply free surface condition to the face define in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyFreeSurfaceBC( real64 const time, DomainPartition & domain ) override;

  /**
   * @brief Apply Perfectly Matched Layer (PML) to the regions defined in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyPML( real64 const time, DomainPartition & domain ) override;

  /// save the sismo trace in file
  void saveSeismo( localIndex iseismo, real32 valDisplacement, string const & filename ) override;

  localIndex getNumNodesPerElem();

  /// Indices of the nodes (in the right order) for each source point
  array2d< localIndex > m_sourceNodeIds;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds
  array2d< real64 > m_sourceConstants;

  /// Flag that indicates whether the source is local or not to the MPI rank
  array1d< localIndex > m_sourceIsAccessible;

  /// Indices of the element nodes (in the right order) for each receiver point
  array2d< localIndex > m_receiverNodeIds;

  /// Basis function evaluated at the receiver for the nodes listed in m_receiverNodeIds
  array2d< real64 > m_receiverConstants;

  /// Flag that indicates whether the receiver is local or not to the MPI rank
  array1d< localIndex > m_receiverIsLocal;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_displacementxNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_displacementyNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_displacementzNp1AtReceivers;

  /// Array containing the elements which contain a source
  array1d< localIndex > m_sourceElem;


};


namespace fields
{

DECLARE_FIELD( Displacementx_np1,
               "displacementx_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "x-component of displacement at time n+1." );

DECLARE_FIELD( Displacementy_np1,
               "displacementy_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "y-component of displacement at time n+1." );

DECLARE_FIELD( Displacementz_np1,
               "displacementz_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "z-component of displacement at time n+1." );

DECLARE_FIELD( Stresstensorxx,
               "stresstensorxx",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "xx-components of the stress tensor." );

DECLARE_FIELD( Stresstensoryy,
               "stresstensoryy",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "yy-components of the stress tensor." );

DECLARE_FIELD( Stresstensorzz,
               "stresstensorzz",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "zz-components of the stress tensor." );

DECLARE_FIELD( Stresstensorxy,
               "stresstensorxy",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "xy-components of the stress tensor (symetric of yx-component)." );

DECLARE_FIELD( Stresstensorxz,
               "stresstensorxz",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "xz-components of the stress tensor (symetric of zx-component)." );

DECLARE_FIELD( Stresstensoryz,
               "stresstensoryz",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "yz-components of the stress tensor (symetric of zy-component)." );


DECLARE_FIELD( ForcingRHS,
               "rhs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS" );

DECLARE_FIELD( MassVector,
               "massVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal Mass Matrix." );

DECLARE_FIELD( DampingVectorx,
               "dampingVectorx",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal Damping Matrix in x-direction." );

DECLARE_FIELD( DampingVectory,
               "dampingVectory",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal Damping Matrix in y-direction." );

DECLARE_FIELD( DampingVectorz,
               "dampingVectorz",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal Damping Matrix in z-direction." );

DECLARE_FIELD( MediumVelocityVp,
               "mediumVelocityVp",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "P-waves speed in the cell" );

DECLARE_FIELD( MediumVelocityVs,
               "mediumVelocityVs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "S-waves speed in the cell" );

DECLARE_FIELD( MediumDensity,
               "mediumDensity",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium density of the cell" );

DECLARE_FIELD( Lambda,
               "lambda",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "First Lame parameter: lambda" );

DECLARE_FIELD( Mu,
               "mu",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Second Lame parameter: mu" );

DECLARE_FIELD( FreeSurfaceFaceIndicator,
               "freeSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

DECLARE_FIELD( FreeSurfaceNodeIndicator,
               "freeSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a node is on free surface 0 otherwise." );

}


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASSTICWAVEEQUATIONSEM_HPP_ */
