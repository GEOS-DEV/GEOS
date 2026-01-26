/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A / TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERANISOTROPICWAVEEQUATIONSEM_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERANISOTROPICWAVEEQUATIONSEM_HPP_

#include "physicsSolvers/wavePropagation/shared/WaveSolverBase.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticFields.hpp"
#include "AcousticVTIFields.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverTypeDefSEM.hpp"


namespace geos
{
// Forward declarations for types used in method signatures
namespace finiteElement
{
class FiniteElementBase;
}
class CellElementSubRegion;

/**
 * @brief Abstract base class for second-order SEM-based acoustic anisotropic wave equation solvers.
 *
 * This class implements all the common logic for VTI/TTI and Fletcher/Zhang solvers.
 * Derived classes are required to provide a specific implementation for calculating the stiffness contribution.
 */
class AcousticSecondOrderAnisotropicWaveEquationSEM : public WaveSolverBase
{
public:
  using EXEC_POLICY = parallelDevicePolicy<>;
  using ATOMIC_POLICY = AtomicPolicy< EXEC_POLICY >;

  AcousticSecondOrderAnisotropicWaveEquationSEM( const std::string & name, Group * const parent );

  virtual ~AcousticSecondOrderAnisotropicWaveEquationSEM() override = default;

  // Prevent copy/move to avoid slicing and ensure unique resource ownership
  AcousticSecondOrderAnisotropicWaveEquationSEM( const AcousticSecondOrderAnisotropicWaveEquationSEM & ) = delete;
  AcousticSecondOrderAnisotropicWaveEquationSEM( AcousticSecondOrderAnisotropicWaveEquationSEM && ) = delete;
  AcousticSecondOrderAnisotropicWaveEquationSEM & operator=( const AcousticSecondOrderAnisotropicWaveEquationSEM & ) = delete;
  AcousticSecondOrderAnisotropicWaveEquationSEM & operator=( AcousticSecondOrderAnisotropicWaveEquationSEM && ) = delete;

  // Common Public Interface Methods
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual real64 explicitStepForward( real64 const & time_n, real64 const & dt, integer const cycleNumber, DomainPartition & domain, integer const computeGradient ) override;
  virtual real64 explicitStepBackward( real64 const & time_n, real64 const & dt, integer const cycleNumber, DomainPartition & domain, integer const computeGradient ) override;
  virtual real32 getGlobalMinWavespeed( MeshLevel & mesh, string_array const & regionNames ) override;
  /**@}*/

  /**
   * @brief Multiply the precomputed term by the Ricker and add to the right-hand side
   * @param cycleNumber the cycle number/step number of evaluation of the source
   * @param rhs the right hand side vector to be computed
   */
  virtual void addSourceToRightHandSide( real64 const & time_n, arrayView1d< real32 > const rhs );


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
                               DomainPartition & domain,
                               bool const isForward );

  void computeUnknowns( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain,
                        MeshLevel & mesh,
                        string_array const & regionNames,
                        bool const isForward );

  void synchronizeUnknowns( real64 const & time_n,
                            real64 const & dt,
                            DomainPartition & domain,
                            MeshLevel & mesh,
                            string_array const & regionNames );

  void prepareNextTimestep( MeshLevel & mesh );


  /**
   * @brief Initialize Perfectly Matched Layer (PML) information
   */
  virtual void initializePML() override;

  /**
   */
  virtual real64 computeTimeStep( real64 & dtOut ) override;

  struct viewKeyStruct : WaveSolverBase::viewKeyStruct
  {
    static constexpr char const * pressureNp1AtReceiversString() { return "pressureNp1AtReceivers"; }
  } waveEquationViewKeys;

protected:
  virtual void postInputInitialization() override;
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & baseMesh, MeshLevel & mesh, string_array const & regionNames ) override;
  virtual void initializePostInitialConditionsPreSubGroups() override;

  /**
   * @brief Virtual method for initializing mass and damping matrices.
   * Different approximation methods (Fletcher vs Zhang) have different implementations.
   */
  virtual void initializeMatrices( MeshLevel & mesh, string_array const & regionNames ) = 0;


  /**
   * @brief Pure virtual function to be implemented by derived classes.
   *
   * This function is the primary customization point. It is responsible for applying the
   * specific spectral element kernel (e.g., Zhang, Fletcher) to compute the
   * stiffness matrix contribution for the wave equations.
   */
  virtual void applyStiffnessKernels( const real64 & dt, MeshLevel & mesh, const string_array & regionNames ) = 0;


  /**
   * @brief Template method for common matrix initialization (mass + VTI DOF arrays + damping)
   * @tparam DampingComputer The specific damping computer to use (Zhang or Fletcher)
   * @param mesh The mesh level
   * @param regionNames The region names to process
   */
  template< typename DampingComputer >
  void initializeMatricesTemplate( MeshLevel & mesh, string_array const & regionNames );

private:

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
  /**
   * @brief Precompute surface field indicators for boundary conditions.
   */
  virtual void precomputeSurfaceFieldIndicator( DomainPartition & domain );

  /// Pressure_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_pressureNp1AtReceivers;

  /// Array of size the number of receivers, used for calculating seismograms
  array1d< real32 > m_seismoCoeff;
};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERANISOTROPICWAVEEQUATIONSEM_HPP_
