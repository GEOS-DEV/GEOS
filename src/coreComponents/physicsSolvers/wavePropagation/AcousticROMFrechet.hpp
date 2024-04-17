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
 * @file AcousticROMFrechet.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICROMFRECHET_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICROMFRECHET_HPP_

#include "WaveSolverBase.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/SolverBase.hpp"

namespace geos
{

class AcousticROMFrechet : public WaveSolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy<  >;
  using ATOMIC_POLICY = AtomicPolicy< EXEC_POLICY >;

  AcousticROMFrechet( const std::string & name,
                           Group * const parent );

  virtual ~AcousticROMFrechet() override;

  AcousticROMFrechet() = delete;
  AcousticROMFrechet( AcousticROMFrechet const & ) = delete;
  AcousticROMFrechet( AcousticROMFrechet && ) = default;

  AcousticROMFrechet & operator=( AcousticROMFrechet const & ) = delete;
  AcousticROMFrechet & operator=( AcousticROMFrechet && ) = delete;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "acoustic"; }

  static string catalogName() { return "AcousticFrechet"; }
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
   * @param rhs the right hand side vector to be computed
   */
  virtual void addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs );


  /**
   * @brief Initialize Perfectly Matched Layer (PML) information
   */
  virtual void initializePML() override;


  /**
   * @brief Overridden from ExecutableGroup. Used to write last seismogram if needed.
   */
  virtual void cleanup( real64 const time_n, integer const cycleNumber, integer const eventCounter, real64 const eventProgress, DomainPartition & domain ) override;

  struct viewKeyStruct : WaveSolverBase::viewKeyStruct
  {
    static constexpr char const * pressureNp1AtReceiversString() { return "pressureNp1AtReceivers"; }
    static constexpr char const * orderFrechetString() { return "orderFrechet"; }
    static constexpr char const * orderGSString() { return "orderGS"; }
    static constexpr char const * epsilonGSString() { return "epsilonGS"; }
    static constexpr char const * count_qString() { return "count_q"; }
    static constexpr char const * totcount_qString() { return "totcount_q"; }
    static constexpr char const * selectionOrderString() { return "selectionOrder"; }
    static constexpr char const * cycleOrderString() { return "cycleOrder"; }
    
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

  bool gramSchmidtROMStiffness(finiteElement::FiniteElementBase const & fe,
                               arrayView1d< real32 const > const Ku,
                               arrayView1d< real32 const > const u,
                               arrayView1d< integer const > const nodeghostrank,
                               localIndex const elemRegionSize,
                               arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                               arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X,
                               localIndex const ordF);

  void gramSchmidtROMStiffnessFinal(finiteElement::FiniteElementBase const & fe,
				    arrayView1d< integer const > const nodeghostrank,
				    localIndex const elemRegionSize,
				    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
				    arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X);
  
  void reorthogonalization(finiteElement::FiniteElementBase const & fe,
			   arrayView1d< integer const > const nodeghostrank,
			   localIndex const elemRegionSize,
			   arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
			   arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X,
			   localIndex const nq,
			   std::string path);
  
  void writeInitialConditionsPOD(arrayView1d< real32 const > const stiffnessVector,
				 int const ordF,
				 int const cycle);

  void writeInitialConditionsPODFinal(arrayView1d< integer const > const nodeGhostRank);
  
protected:

  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  /**
   * @brief Locate sources and receivers position in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param mesh mesh of the computational domain
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


  /// Pressure_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_pressureNp1AtReceivers;
  localIndex m_orderFrechet;
  localIndex m_orderGS;
  real32 m_epsilonGS;
  array1d< int > m_count_q;
  localIndex m_totcount_q;
  array2d< int > m_selectionOrder;
  array2d< int > m_cycleOrder;
};


namespace fields
{

DECLARE_FIELD( Pressure_nm1,
               "pressure_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n-1." );

DECLARE_FIELD( Pressure_n,
               "pressure_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n." );

DECLARE_FIELD( Pressure_np1,
               "pressure_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar pressure at time n+1." );

DECLARE_FIELD( PressureFrechet_nm1,
               "pressureFrechet_nm1",
               array2d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar Frechet derivative of pressure at time n-1." );

DECLARE_FIELD( PressureFrechet_n,
               "pressureFrechet_n",
               array2d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar Frechet derivative of pressure at time n." );

DECLARE_FIELD( PressureFrechet_np1,
               "pressureFrechet_np1",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar Frechet derivative of pressure at time n+1." );

DECLARE_FIELD( ForcingRHS,
               "rhs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS" );

DECLARE_FIELD( ForcingRHS_fp1,
               "rhs_fp1",
               array1d< real32 >,
	       0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS_fp1" );
  
DECLARE_FIELD( AcousticMassVectorFrechet,
               "acousticMassVectorFrechet",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Mass Matrix with gradient." );

DECLARE_FIELD( AcousticMassVector,
               "acousticMassVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Mass Matrix." );

DECLARE_FIELD( DampingVector,
               "dampingVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix." );

DECLARE_FIELD( AcousticVelocity,
               "acousticVelocity",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium velocity of the cell" );

DECLARE_FIELD( AcousticDensity,
               "acousticDensity",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium density of the cell" );

DECLARE_FIELD( StiffnessVector,
               "stiffnessVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Stiffness vector contains R_h*Pressure_n." );

DECLARE_FIELD( AcousticFreeSurfaceFaceIndicator,
               "acousticFreeSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

DECLARE_FIELD( AcousticFreeSurfaceNodeIndicator,
               "acousticFreeSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a node is on free surface 0 otherwise." );

DECLARE_FIELD( PressureDoubleDerivative,
               "pressureDoubleDerivative",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Double derivative of the pressure for each node to compute the gradient" );

DECLARE_FIELD( PartialGradient,
               "partialGradient",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Partiel gradient computed during backward propagation" );

DECLARE_FIELD( AuxiliaryVar1PML,
               "auxiliaryVar1PML",
               array2d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "PML vectorial auxiliary variable 1." );

DECLARE_FIELD( AuxiliaryVar2PML,
               "auxiliaryVar2PML",
               array2d< real32 >,
               0,
               NOPLOT,
               NO_WRITE,
               "PML vectorial auxiliary variable 2." );

DECLARE_FIELD( AuxiliaryVar3PML,
               "auxiliaryVar3PML",
               array1d< real32 >,
               0,
               NOPLOT,
               NO_WRITE,
               "PML scalar auxiliary variable 3." );

DECLARE_FIELD( AuxiliaryVar4PML,
               "auxiliaryVar4PML",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "PML scalar auxiliary variable 4." );

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICROMFRECHET_HPP_ */
