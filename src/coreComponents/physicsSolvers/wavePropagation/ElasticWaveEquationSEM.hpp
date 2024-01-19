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
 * @file ElasticWaveEquationSEM.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEM_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEM_HPP_

#include "WaveSolverBase.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "ElasticFields.hpp"

namespace geos
{

class ElasticWaveEquationSEM : public WaveSolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy<  >;
  using ATOMIC_POLICY = parallelDeviceAtomic;

  /**
   * @brief Safeguard for timeStep. Used to avoid memory issue due to too small value.
   */
  static constexpr real64 epsilonLoc = 1e-8;

  ElasticWaveEquationSEM( const std::string & name,
                          Group * const parent );

  virtual ~ElasticWaveEquationSEM() override;

  ElasticWaveEquationSEM() = delete;
  ElasticWaveEquationSEM( ElasticWaveEquationSEM const & ) = delete;
  ElasticWaveEquationSEM( ElasticWaveEquationSEM && ) = default;

  ElasticWaveEquationSEM & operator=( ElasticWaveEquationSEM const & ) = delete;
  ElasticWaveEquationSEM & operator=( ElasticWaveEquationSEM && ) = delete;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "elastic"; }

  static string catalogName() { return "ElasticSEM"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

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
   * @param rhsx the right hand side vector to be computed (x-component)
   * @param rhsy the right hand side vector to be computed (y-component)
   * @param rhsz the right hand side vector to be computed (z-component)
   */
  void addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhsx, arrayView1d< real32 > const rhsy, arrayView1d< real32 > const rhsz );

  /**
   * TODO: move implementation into WaveSolverBase once 'm_receiverIsLocal' is also moved
   * @brief Compute DAS data from the appropriate three-component receiver pairs
   * @param xCompRcv the array holding the x-component of pairs of receivers
   * @param yCompRcv the array holding the y-component of pairs of receivers
   * @param zCompRcv the array holding the z-component of pairs of receivers
   */
  void computeDAS( arrayView2d< real32 > const xCompRcv,
                   arrayView2d< real32 > const yCompRcv,
                   arrayView2d< real32 > const zCompRcv );


  /**
   * @brief Overridden from ExecutableGroup. Used to write last seismogram if needed.
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  struct viewKeyStruct : WaveSolverBase::viewKeyStruct
  {
    static constexpr char const * displacementXNp1AtReceiversString() { return "displacementXNp1AtReceivers"; }
    static constexpr char const * displacementYNp1AtReceiversString() { return "displacementYNp1AtReceivers"; }
    static constexpr char const * displacementZNp1AtReceiversString() { return "displacementZNp1AtReceivers"; }

    static constexpr char const * sourceForceString() { return "sourceForce"; }
    static constexpr char const * sourceMomentString() { return "sourceMoment"; }

  } waveEquationViewKeys;


  /** internal function to the class to compute explicitStep either for backward or forward.
   * (requires not to be private because it is called from GEOS_HOST_DEVICE method)
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

  void computeUnknowns( real64 const & time_n,
                        real64 const & dt,
                        integer const cycleNumber,
                        DomainPartition & domain,
                        MeshLevel & mesh,
                        arrayView1d< string const > const & regionNames );

  void synchronizeUnknowns( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain,
                            MeshLevel & mesh,
                            arrayView1d< string const > const & regionNames );

  void prepareNextTimestep( MeshLevel & mesh );
protected:

  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  /**
   * @brief Locate sources and receivers position in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param mesh mesh of the computational domain
   * @param regionNames the names of the region you loop on
   */
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) override;

  /**
   * @brief Apply free surface condition to the face define in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyFreeSurfaceBC( real64 const time, DomainPartition & domain ) override;

  /**
   * @brief Initialize Perfectly Matched Layer (PML) information
   */
  virtual void initializePML() override;

  /**
   * @brief Apply Perfectly Matched Layer (PML) to the regions defined in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyPML( real64 const time, DomainPartition & domain ) override;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds in x-direction
  array2d< real64 > m_sourceConstantsx;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds in y-direction
  array2d< real64 > m_sourceConstantsy;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds in z-direction
  array2d< real64 > m_sourceConstantsz;

  /// Displacement_np1 at the receiver location for each time step for each receiver (x-component)
  array2d< real32 > m_displacementXNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver (y-component)
  array2d< real32 > m_displacementYNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver (z-component)
  array2d< real32 > m_displacementZNp1AtReceivers;

  /// Vector describing the force of the source
  R1Tensor m_sourceForce;

  /// Symmetric tensor describing the moment of the source
  R2SymTensor m_sourceMoment;

};


// namespace fields
// {

// DECLARE_FIELD( Displacementx_nm1,
//                "displacementx_nm1",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "x-component of displacement at time n-1." );

// DECLARE_FIELD( Displacementy_nm1,
//                "displacementy_nm1",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "y-component of displacement at time n-1." );

// DECLARE_FIELD( Displacementz_nm1,
//                "displacementz_nm1",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "z-component of displacement at time n-1." );

// DECLARE_FIELD( Displacementx_n,
//                "displacementx_n",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "x-component of displacement at time n." );

// DECLARE_FIELD( Displacementy_n,
//                "displacementy_n",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "y-component of displacement at time n." );

// DECLARE_FIELD( Displacementz_n,
//                "displacementz_n",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "z-component of displacement at time n." );

// DECLARE_FIELD( Displacementx_np1,
//                "displacementx_np1",
//                array1d< real32 >,
//                0,
//                LEVEL_0,
//                WRITE_AND_READ,
//                "x-component of displacement at time n+1." );

// DECLARE_FIELD( Displacementy_np1,
//                "displacementy_np1",
//                array1d< real32 >,
//                0,
//                LEVEL_0,
//                WRITE_AND_READ,
//                "y-component of displacement at time n+1." );

// DECLARE_FIELD( Displacementz_np1,
//                "displacementz_np1",
//                array1d< real32 >,
//                0,
//                LEVEL_0,
//                WRITE_AND_READ,
//                "z-component of displacement at time n+1." );

// DECLARE_FIELD( ForcingRHSx,
//                "rhsx",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "RHS for x-direction" );

// DECLARE_FIELD( ForcingRHSy,
//                "rhsy",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "RHS for y-direction" );

// DECLARE_FIELD( ForcingRHSz,
//                "rhsz",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "RHS for z-direction" );

// DECLARE_FIELD( ElasticMassVector,
//                "elasticMassVector",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "Diagonal Mass Matrix." );

// DECLARE_FIELD( DampingVectorx,
//                "dampingVectorx",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "Diagonal Damping Matrix in x-direction." );

// DECLARE_FIELD( DampingVectory,
//                "dampingVectory",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "Diagonal Damping Matrix in y-direction." );

// DECLARE_FIELD( DampingVectorz,
//                "dampingVectorz",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "Diagonal Damping Matrix in z-direction." );

// DECLARE_FIELD( StiffnessVectorx,
//                "stiffnessVectorx",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "x-component of stiffness vector." );

// DECLARE_FIELD( StiffnessVectory,
//                "stiffnessVectory",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "y-component of stiffness vector." );

// DECLARE_FIELD( StiffnessVectorz,
//                "stiffnessVectorz",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "z-component of stiffness vector." );

// DECLARE_FIELD( ElasticVelocityVp,
//                "elasticVelocityVp",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "P-waves speed in the cell" );

// DECLARE_FIELD( ElasticVelocityVs,
//                "elasticVelocityVs",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "S-waves speed in the cell" );

// DECLARE_FIELD( ElasticDensity,
//                "elasticDensity",
//                array1d< real32 >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "Medium density of the cell" );

// DECLARE_FIELD( ElasticFreeSurfaceFaceIndicator,
//                "elasticFreeSurfaceFaceIndicator",
//                array1d< localIndex >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

// DECLARE_FIELD( ElasticFreeSurfaceNodeIndicator,
//                "elasticFreeSurfaceNodeIndicator",
//                array1d< localIndex >,
//                0,
//                NOPLOT,
//                WRITE_AND_READ,
//                "Free surface indicator, 1 if a node is on free surface 0 otherwise." );

// }


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASSTICWAVEEQUATIONSEM_HPP_ */
