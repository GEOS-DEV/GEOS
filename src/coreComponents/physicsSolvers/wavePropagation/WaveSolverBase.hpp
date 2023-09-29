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
 * @file WaveSolverBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERBASE_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERBASE_HPP_


#include "mesh/MeshFields.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "common/LifoStorage.hpp"
#if !defined( GEOS_USE_HIP )
#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#endif
#include "WaveSolverUtils.hpp"

#if !defined( GEOS_USE_HIP )
#define SEM_FE_TYPES \
  finiteElement::Q1_Hexahedron_Lagrange_GaussLobatto, \
  finiteElement::Q2_Hexahedron_Lagrange_GaussLobatto, \
  finiteElement::Q3_Hexahedron_Lagrange_GaussLobatto, \
  finiteElement::Q4_Hexahedron_Lagrange_GaussLobatto, \
  finiteElement::Q5_Hexahedron_Lagrange_GaussLobatto
#else
#define SEM_FE_TYPES
#endif

#define SELECTED_FE_TYPES SEM_FE_TYPES

namespace geos
{

class WaveSolverBase : public SolverBase
{
public:

  static constexpr real64 epsilonLoc = WaveSolverUtils::epsilonLoc;
  using EXEC_POLICY = WaveSolverUtils::EXEC_POLICY;
  using wsCoordType = WaveSolverUtils::wsCoordType;

  WaveSolverBase( const std::string & name,
                  Group * const parent );

  virtual ~WaveSolverBase() override;

  WaveSolverBase() = delete;
  WaveSolverBase( WaveSolverBase const & ) = delete;
  WaveSolverBase( WaveSolverBase && ) = default;

  WaveSolverBase & operator=( WaveSolverBase const & ) = delete;
  WaveSolverBase & operator=( WaveSolverBase && ) = delete;

  virtual void initializePreSubGroups() override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;


  virtual real64 explicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr char const * sourceCoordinatesString() { return "sourceCoordinates"; }
    static constexpr char const * sourceValueString() { return "sourceValue"; }

    static constexpr char const * timeSourceFrequencyString() { return "timeSourceFrequency"; }

    static constexpr char const * receiverCoordinatesString() { return "receiverCoordinates"; }

    static constexpr char const * sourceNodeIdsString() { return "sourceNodeIds"; }
    static constexpr char const * sourceConstantsString() { return "sourceConstants"; }
    static constexpr char const * sourceIsAccessibleString() { return "sourceIsAccessible"; }

    static constexpr char const * receiverNodeIdsString() { return "receiverNodeIds"; }
    static constexpr char const * receiverConstantsString() {return "receiverConstants"; }
    static constexpr char const * receiverIsLocalString() { return "receiverIsLocal"; }

    static constexpr char const * rickerOrderString() { return "rickerOrder"; }
    static constexpr char const * outputSeismoTraceString() { return "outputSeismoTrace"; }
    static constexpr char const * dtSeismoTraceString() { return "dtSeismoTrace"; }
    static constexpr char const * indexSeismoTraceString() { return "indexSeismoTrace"; }
    static constexpr char const * forwardString() { return "forward"; }
    static constexpr char const * saveFieldsString() { return "saveFields"; }
    static constexpr char const * shotIndexString() { return "shotIndex"; }
    static constexpr char const * enableLifoString() { return "enableLifo"; }
    static constexpr char const * lifoSizeString() { return "lifoSize"; }
    static constexpr char const * lifoOnDeviceString() { return "lifoOnDevice"; }
    static constexpr char const * lifoOnHostString() { return "lifoOnHost"; }

    static constexpr char const * useDASString() { return "useDAS"; }
    static constexpr char const * linearDASGeometryString() { return "linearDASGeometry"; }

    static constexpr char const * usePMLString() { return "usePML"; }
    static constexpr char const * parametersPMLString() { return "parametersPML"; }

    static constexpr char const * receiverElemString() { return "rcvElem"; }
    static constexpr char const * receiverRegionString() { return "receiverRegion"; }
  };

  /**
   * @brief Re-initialize source and receivers positions in the mesh, and resize the pressureNp1_at_receivers array
   */
  void reinit() override final;

protected:

  virtual void postProcessInput() override;

  /**
   * @brief Apply free surface condition to the face defined in the geometry box of the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyFreeSurfaceBC( real64 const time, DomainPartition & domain ) = 0;

  /**
   * @brief Initialize DAS fiber geometry. This will duplicate the number of point receivers to be modeled
   */
  virtual void initializeDAS();

  /**
   * @brief Initialize Perfectly Matched Layer (PML) information
   */
  virtual void initializePML() = 0;

  virtual void incrementIndexSeismoTrace( real64 const time_n );

  /**
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
   * @brief Apply Perfectly Matched Layer (PML) to the regions defined in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyPML( real64 const time, DomainPartition & domain ) = 0;

  /**
   * @brief Locate sources and receivers positions in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param mesh mesh of the computational domain
   */
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) = 0;

  /**
   * @brief Perform forward explicit step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @param computeGradient Indicates if we want to compute gradient at this step
   * @return return the timestep that was achieved during the step.
   */
  virtual real64 explicitStepForward( real64 const & time_n,
                                      real64 const & dt,
                                      integer const cycleNumber,
                                      DomainPartition & domain,
                                      bool const computeGradient ) = 0;
  /**
   * @brief Perform backward explicit step
   * @param time_n time at the beginning of the step
   * @param dt the perscribed timestep
   * @param cycleNumber the current cycle number
   * @param domain the domain object
   * @param computeGradient Indicates if we want to compute gradient at this step
   * @return return the timestep that was achieved during the step.
   */
  virtual real64 explicitStepBackward( real64 const & time_n,
                                       real64 const & dt,
                                       integer const cycleNumber,
                                       DomainPartition & domain,
                                       bool const computeGradient ) = 0;


  virtual void registerDataOnMesh( Group & meshBodies ) override;

  localIndex getNumNodesPerElem();

  /// Coordinates of the sources in the mesh
  array2d< real64 > m_sourceCoordinates;

  /// Precomputed value of the source terms
  array2d< real32 > m_sourceValue;

  /// Central frequency for the Ricker time source
  real32 m_timeSourceFrequency;

  /// Coordinates of the receivers in the mesh
  array2d< real64 > m_receiverCoordinates;

  /// Flag that indicates the order of the Ricker to be used, order 2 by default
  localIndex m_rickerOrder;

  /// Flag that indicates if we write the seismo trace in a file .txt, 0 no output, 1 otherwise
  integer m_outputSeismoTrace;

  /// Time step for seismoTrace output
  real64 m_dtSeismoTrace;

  /// Cycle number for output SeismoTrace
  localIndex m_indexSeismoTrace;

  /// Amount of seismoTrace that will be recorded for each receiver
  localIndex m_nsamplesSeismoTrace;

  /// Flag to indicate if DAS type of data will be modeled
  integer m_useDAS;

  /// Geometry parameters for a linear DAS fiber (dip, azimuth, gauge length)
  array2d< real64 > m_linearDASGeometry;

  /// Indicate if we want to compute forward ou backward
  localIndex m_forward;

  /// Indicate if we want to save fields to restore them during backward
  localIndex m_saveFields;

  // Indicate the current shot computed for naming saved temporary data
  integer m_shotIndex;

  /// Flag to apply PML
  integer m_usePML;

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

  /// Array containing the elements which contain a receiver
  array1d< localIndex > m_rcvElem;

  /// Array containing the elements which contain the region which the receiver belongs
  array1d< localIndex > m_receiverRegion;

  /// Flag to enable LIFO
  localIndex m_enableLifo;

  /// lifo size (should be the total number of buffer to save in LIFO)
  localIndex m_lifoSize;

  /// Number of buffers to store on device by LIFO  (if negative, opposite of percentage of remaining memory)
  localIndex m_lifoOnDevice;

  /// Number of buffers to store on host by LIFO  (if negative, opposite of percentage of remaining memory)
  localIndex m_lifoOnHost;

  /// LIFO to store p_dt2
  std::unique_ptr< LifoStorage< real32, localIndex > > m_lifo;

  struct parametersPML
  {
    /// Mininum (x,y,z) coordinates of inner PML boundaries
    R1Tensor32 xMinPML;

    /// Maximum (x,y,z) coordinates of inner PML boundaries
    R1Tensor32 xMaxPML;

    /// Desired reflectivity of the PML region, used to compute the damping profile
    real32 reflectivityPML;

    /// Thickness of the PML region, used to compute the damping profile
    R1Tensor32 thicknessMinXYZPML;
    R1Tensor32 thicknessMaxXYZPML;

    /// Wave speed in the PML region, used to compute the damping profile
    R1Tensor32 waveSpeedMinXYZPML;
    R1Tensor32 waveSpeedMaxXYZPML;
  };

};

namespace fields
{
using reference32Type = array2d< WaveSolverUtils::wsCoordType, nodes::REFERENCE_POSITION_PERM >;
DECLARE_FIELD( referencePosition32,
               "referencePosition32",
               reference32Type,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Copy of the referencePosition from NodeManager in 32 bits integer" );
}
} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERBASE_HPP_ */
