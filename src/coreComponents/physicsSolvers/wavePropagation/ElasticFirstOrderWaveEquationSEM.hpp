/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOS Contributors
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
#include "WaveSolverUtils.hpp"
#include "WaveSolverBaseFields.hpp"



namespace geos
{

class ElasticFirstOrderWaveEquationSEM : public WaveSolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy< >;
  using ATOMIC_POLICY = parallelDeviceAtomic;

  static constexpr real64 epsilonLoc = 1e-8;

  ElasticFirstOrderWaveEquationSEM( const std::string & name,
                                    Group * const parent );

  virtual ~ElasticFirstOrderWaveEquationSEM() override;

  static string catalogName() { return "ElasticFirstOrderSEM"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override  { return catalogName(); }

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
   * TODO: move implementation into WaveSolverBase
   * @brief Computes the traces on all receivers (see @computeSeismoTraces) up to time_n+dt
   * @param time_n the time corresponding to the field values pressure_n
   * @param dt the simulation timestep
   * @param var_np1 the field values at time_n + dt
   * @param var_n the field values at time_n
   * @param varAtreceivers the array holding the trace values, where the output is written
   */
  virtual void compute2dVariableAllSeismoTraces( localIndex const regionIndex,
                                                 real64 const time_n,
                                                 real64 const dt,
                                                 arrayView2d< real32 const > const var_np1,
                                                 arrayView2d< real32 const > const var_n,
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
    static constexpr char const * displacementxNp1AtReceiversString() { return "displacementxNp1AtReceivers"; }
    static constexpr char const * displacementyNp1AtReceiversString() { return "displacementyNp1AtReceivers"; }
    static constexpr char const * displacementzNp1AtReceiversString() { return "displacementzNp1AtReceivers"; }

    static constexpr char const * sigmaxxNp1AtReceiversString() { return "sigmaxxNp1AtReceivers"; }
    static constexpr char const * sigmayyNp1AtReceiversString() { return "sigmayyNp1AtReceivers"; }
    static constexpr char const * sigmazzNp1AtReceiversString() { return "sigmazzNp1AtReceivers"; }
    static constexpr char const * sigmaxyNp1AtReceiversString() { return "sigmaxyNp1AtReceivers"; }
    static constexpr char const * sigmaxzNp1AtReceiversString() { return "sigmaxzNp1AtReceivers"; }
    static constexpr char const * sigmayzNp1AtReceiversString() { return "sigmayzNp1AtReceivers"; }

    static constexpr char const * sourceElemString() { return "sourceElem"; }
    static constexpr char const * sourceRegionString() { return "sourceRegion"; }
    static constexpr char const * receiverElemString() { return "rcvElem"; }
    static constexpr char const * receiverRegionString() { return "receiverRegion"; }

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

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_displacementxNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_displacementyNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_displacementzNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_sigmaxxNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_sigmayyNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_sigmazzNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_sigmaxyNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_sigmaxzNp1AtReceivers;

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_sigmayzNp1AtReceivers;

  /// Array containing the elements which contain a source
  array1d< localIndex > m_sourceElem;

  /// Array containing the elements which contain the region which the source belongs
  array1d< localIndex > m_sourceRegion;

  /// Array containing the elements which contain a receiver
  array1d< localIndex > m_rcvElem;

  /// Array containing the elements which contain the region which the receiver belongs
  array1d< localIndex > m_receiverRegion;

};

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASSTICWAVEEQUATIONSEM_HPP_ */
