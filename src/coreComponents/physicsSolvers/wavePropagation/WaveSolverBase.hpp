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

#include "mesh/ExtrinsicMeshData.hpp"
#include "physicsSolvers/SolverBase.hpp"


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

    static constexpr char const * xMinPMLString() { return "xMinPML"; }
    static constexpr char const * xMaxPMLString() { return "xMaxPML"; }
    static constexpr char const * yMinPMLString() { return "yMinPML"; }
    static constexpr char const * yMaxPMLString() { return "yMaxPML"; }
    static constexpr char const * zMinPMLString() { return "zMinPML"; }
    static constexpr char const * zMaxPMLString() { return "zMaxPML"; }
    static constexpr char const * maxThicknessPMLString() { return "maxThicknessPML"; }
    static constexpr char const * reflectivityPMLString() { return "reflectivityPML"; }
    static constexpr char const * flagPMLString() { return "flagPML"; }

  };

  /**
   * @brief Re-initialize source and receivers positions in the mesh, and resize the pressureNp1_at_receivers array
   */
  void reinit() override final;

protected:

  /**
   * @brief Apply free surface condition to the face defined in the geometry box of the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyFreeSurfaceBC( real64 const time, DomainPartition & domain ) = 0;


  /**
   * @brief Apply Perfectly Matched Layer (PML) to the regions defined in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyPML( real64 const time, DomainPartition & domain ) = 0;

  /**
   * @brief Compute the damping profile for the Perfectly Matched Layer (PML)
   * @param xLocal a given x-y-z coordinates (3-components array)
   * @param c wave speed at the location xLocal
   * @param xMin etc... coordinate limits used to compute the damping profile
   * @param dPML maximum thickness used to compute the damping profile
   * @param rPML desired reflectivity of the PMLs
   * @param sigma 3-components array to hold the damping profile in each direction
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void computeDampingProfilePML( real64 const (&xLocal)[3],
                                        real64 const c,
                                        real64 const xMin,
                                        real64 const xMax,
                                        real64 const yMin,
                                        real64 const yMax,
                                        real64 const zMin,
                                        real64 const zMax,
                                        real64 const dPML,
                                        real64 const rPML,
                                        real64 (&sigma)[3])
  {
    real64 const factor =  -3.0/2.0*c*log(rPML)/(dPML*dPML*dPML);

    sigma[0] = 0;
    sigma[1] = 0;
    sigma[2] = 0;
    if (xLocal[0] < xMin)
    {
      sigma[0] = factor*(xLocal[0]-xMin)*(xLocal[0]-xMin);
    }
    else if (xLocal[0] > xMax)
    {
      sigma[0] = factor*(xLocal[0]-xMax)*(xLocal[0]-xMax);
    }
    if (xLocal[1] < yMin)
    {
      sigma[1] = factor*(xLocal[1]-yMin)*(xLocal[1]-yMin);
    }
    else if (xLocal[1] > yMax)
    {
      sigma[1] = factor*(xLocal[1]-yMax)*(xLocal[1]-yMax);
    }
    if (xLocal[2] < zMin)
    {
      sigma[2] = factor*(xLocal[2]-zMin)*(xLocal[2]-zMin);
    }
    else if (xLocal[2] > zMax)
    {
      sigma[2] = factor*(xLocal[2]-zMax)*(xLocal[2]-zMax);
    }
  }

  template<int numNodesPerElem>
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void computeGradientPML( real64 const (&gradN)[numNodesPerElem][3],
                                  real64 const (&var)[numNodesPerElem],
                                  real64 (&varGrad)[3])
  {
    for (int j=0; j<3; ++j)
    {
      varGrad[j]=0;
      for (int i=0; i<numNodesPerElem; ++i)
      {
        varGrad[j] += var[i]*gradN[i][j];
      }
    }
  }


  /**
   * @brief Compute the value of a Ricker (a Gaussian function)
   * @param time_n time to evaluate the Ricker
   * @param f0 central frequency of the Ricker
   * @param order order of the ricker
   * @return the value of a Ricker evaluated a time_n with f0
   */
  virtual
  real64 evaluateRicker( real64 const & time_n, real64 const & f0, localIndex order );

  /**
   * @brief Locate sources and receivers positions in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param mesh mesh of the computational domain
   */
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) = 0;

  /**
   * @brief Multiply the precomputed term by the Ricker and add to the right-hand side
   * @param time_n the time of evaluation of the source
   * @param rhs the right hand side vector to be computed
   */
  virtual void addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real64 > const rhs ) = 0;

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
                                   arrayView1d< real64 const > const var_np1,
                                   arrayView1d< real64 const > const var_n,
                                   arrayView2d< real64 > varAtReceivers ) = 0;

  /**
   * @brief Temporary debug function. Saves the sismo trace to a file.
   * @param iSeismo index number of the seismo trace
   * @param val value to be written in seismo
   * @param filename name of the output file
   */
  virtual void saveSeismo( localIndex const iSeismo, real64 val, string const & filename ) = 0;



  /// Coordinates of the sources in the mesh
  array2d< real64 > m_sourceCoordinates;

  array2d< real64 > m_sourceValue;

  /// Central frequency for the Ricker time source
  real64 m_timeSourceFrequency;

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

  /// Limits to be used to compute the damping profile for the PMLs
  real64 m_xMinPML;
  real64 m_xMaxPML;
  real64 m_yMinPML;
  real64 m_yMaxPML;
  real64 m_zMinPML;
  real64 m_zMaxPML;

  /// Desired reflectivity and maximum thickness of the PMLs
  real64 m_maxThicknessPML;
  real64 m_reflectivityPML;

  /// Temporary flag for PML testing
  int m_flagPML;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERBASE_HPP_ */
