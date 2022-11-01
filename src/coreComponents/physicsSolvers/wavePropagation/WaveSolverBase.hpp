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

#ifndef GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERBASE_HPP_

#include "mesh/MeshFields.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "finiteElement/elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp"
#include "finiteElement/elementFormulations/Q3_Hexahedron_Lagrange_GaussLobatto.hpp"

#define SEM_FE_TYPES \
  finiteElement::H1_Hexahedron_Lagrange1_GaussLegendre2, \
  finiteElement::Q3_Hexahedron_Lagrange_GaussLobatto

namespace geosx
{

class WaveSolverBase : public SolverBase
{
public:


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

    static constexpr char const * rickerOrderString() { return "rickerOrder"; }
    static constexpr char const * outputSeismoTraceString() { return "outputSeismoTrace"; }
    static constexpr char const * dtSeismoTraceString() { return "dtSeismoTrace"; }
    static constexpr char const * indexSeismoTraceString() { return "indexSeismoTrace"; }
    static constexpr char const * forwardString() { return "forward"; }
    static constexpr char const * saveFieldsString() { return "saveFields"; }
    static constexpr char const * shotIndexString() { return "shotIndex"; }

    static constexpr char const * useDASString() { return "useDAS"; }
    static constexpr char const * linearDASGeometryString() { return "linearDASGeometry"; }

    static constexpr char const * usePMLString() { return "usePML"; }
    static constexpr char const * parametersPMLString() { return "parametersPML"; }

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

  /**
   * @brief Apply Perfectly Matched Layer (PML) to the regions defined in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyPML( real64 const time, DomainPartition & domain ) = 0;


  /**
   * @brief Compute the value of a Ricker (a Gaussian function)
   * @param time_n time to evaluate the Ricker
   * @param f0 central frequency of the Ricker
   * @param order order of the ricker
   * @return the value of a Ricker evaluated a time_n with f0
   */
  virtual real32 evaluateRicker( real64 const & time_n, real32 const & f0, localIndex order );

  /**
   * @brief Locate sources and receivers positions in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param mesh mesh of the computational domain
   */
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) = 0;

  /**
   * @brief Compute the sesimic traces for a given variable at each receiver coordinate at a given time, using the field values at the
   * last two timesteps.
   * @param time_n the time corresponding to the field values pressure_n
   * @param dt the simulation timestep
   * @param timeSeismo the time at which the seismogram is computed
   * @param iSeismo the index of the seismogram time in the seismogram array
   * @param var_at_np1 the field values at time_n + dt
   * @param var_at_n the field values at time_n
   * @param var_at_receivers the array holding the trace values, where the output is written
   */
  virtual void computeSeismoTrace( real64 const time_n,
                                   real64 const dt,
                                   real64 const timeSeismo,
                                   localIndex iSeismo,
                                   arrayView1d< real32 const > const var_np1,
                                   arrayView1d< real32 const > const var_n,
                                   arrayView2d< real32 > varAtReceivers ) = 0;

  /**
   * @brief Temporary debug function. Saves the sismo trace to a file.
   * @param iSeismo index number of the seismo trace
   * @param val value to be written in seismo
   * @param filename name of the output file
   */
  virtual void saveSeismo( localIndex const iSeismo, real32 val, string const & filename ) = 0;



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
  localIndex m_outputSeismoTrace;

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

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERBASE_HPP_ */
