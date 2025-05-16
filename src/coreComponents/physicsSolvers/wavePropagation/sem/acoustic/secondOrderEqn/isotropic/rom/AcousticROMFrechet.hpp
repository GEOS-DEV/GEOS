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

#include "physicsSolvers/wavePropagation/shared/WaveSolverBase.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticFields.hpp"

namespace geos
{

using pFieldType = real32;
  
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
  virtual void initializePML() override {};

  /**
   */
  virtual real64 computeTimeStep( real64 & dtOut ) override {return 0;};

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
    static constexpr char const * sizePOD_fString() { return "sizePOD_f"; }
    static constexpr char const * sizePODString() { return "sizePOD"; }
    static constexpr char const * selectionOrderString() { return "selectionOrder"; }
    static constexpr char const * cycleOrderString() { return "cycleOrder"; }
    static constexpr char const * solverROMString() { return "solverROM"; }
    static constexpr char const * alphaString() { return "alpha"; }
    
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
                               DomainPartition & domain );

  void computeUnknowns( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain,
                        MeshLevel & mesh,
                        arrayView1d< string const > const & regionNames );

  void synchronizeUnknowns( real64 const & time_n,
                            real64 const & dt,
                            DomainPartition & domain,
                            MeshLevel & mesh,
                            arrayView1d< string const > const & regionNames );

  void prepareNextTimestep( MeshLevel & mesh );

  
  bool gramSchmidtROMStiffness(finiteElement::FiniteElementBase const & fe,
                               arrayView1d< pFieldType const > const Ku,
                               arrayView1d< pFieldType const > const u,
                               arrayView1d< integer const > const nodeghostrank,
                               localIndex const elemRegionSize,
                               arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                               arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X,
                               localIndex const ordF);

  void modifiedGramSchmidtROMStiffness(finiteElement::FiniteElementBase const & fe,
				       arrayView1d< integer const > const nodeghostrank,
				       localIndex const elemRegionSize,
				       arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
				       arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X);
  
  void gramSchmidtROMStiffnessFinal(finiteElement::FiniteElementBase const & fe,
				    arrayView1d< integer const > const nodeghostrank,
				    localIndex const elemRegionSize,
				    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
				    arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X);

  int reorthogonalization(finiteElement::FiniteElementBase const & fe,
			  arrayView1d< integer const > const nodeghostrank,
			  localIndex const elemRegionSize,
			  arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
			  arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X,
			  localIndex const nq,
			  std::string path);


  void computeReducedMatrices( arrayView2d< real32 > const massPOD,
			       arrayView2d< real32 > const massPerturbationPOD,
			       arrayView2d< real32 > const dampingPOD,
			       arrayView2d< real32 > const dampingPerturbationPOD,
			       arrayView2d< real64 > const sourceConstantsPOD,
			       arrayView2d< real64 > const receiverConstantsPOD,
			       arrayView1d< real32 const > const mass,
			       arrayView1d< real32 const > const massPerturbation,
			       arrayView1d< real32 const > const damping,
			       arrayView1d< real32 const > const dampingPerturbation,
			       arrayView2d< real64 const > const sourceConstants,
			       arrayView2d< localIndex const > const sourceNodeIds,
			       arrayView1d< localIndex const > const sourceIsAccessible,
			       arrayView2d< real64 const > const receiverConstants,
			       arrayView2d< localIndex const > const receiverNodeIds,
			       arrayView1d< localIndex const > const receiverIsLocal,
			       arrayView1d< localIndex const > const nodesGhostRank);

  void computeSeismoTracePOD( real64 const time_n,
			      real64 const dt,
			      real64 const timeSeismo,
			      localIndex iSeismo,
			      arrayView2d< real64 const > const receiverConstants,
			      arrayView1d< localIndex const > const receiverIsLocal,
			      localIndex const nsamplesSeismoTrace,
			      arrayView1d< real32 const > const var_np1,
			      arrayView1d< real32 const > const var_n,
			      arrayView2d< real32 > varAtReceivers );

protected:

  virtual void postInputInitialization() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  /**
   * @brief Locate sources and receivers position in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param mesh mesh of the computational domain
   */
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & baseMesh, MeshLevel & mesh, arrayView1d< string const > const & regionNames ) override;

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
  virtual void applyPML( real64 const time, DomainPartition & domain ) override {};


  //Pressure_np1 at the receiver location for each time step for each receiver
  array2d< real32 > m_pressureNp1AtReceivers;

  //Order of the Frechet derivatives
  localIndex m_orderFrechet;

  //Order of Gram-Schmidt process
  localIndex m_orderGS;

  //Tolerance for Gram-Schmidt process
  real32 m_epsilonGS;

  //Table of selected selected Frechet order and size by time causality
  array2d< int > m_selectionOrder;

  //Table of selected time step for each order by time causality
  array2d< int > m_cycleOrder;

  //Flag to use the ROM solver
  int m_solverROM;

  //Mass Matrix in the POD basis forward
  array2d< real32 > m_massPOD;

  //Mass Matrix with perturbation in the pod basis
  array2d< real32 > m_massPerturbationPOD;

  //Damping Matrix in the POD basis forward
  array2d< real32 > m_dampingPOD;

  //Damping Matrix with perturbation in the pod basis
  array2d< real32 > m_dampingPerturbationPOD;

  //Distance of perturbation along the gradient
  real32 m_alpha;

  //Perturbation for Frechet derivative computation
  array1d< real32 > m_perturbation;

  //Operator to compute a_np1
  array2d< real64 > m_OpPOD;

  //Coefficient of the solution in the POD basis at time n+1
  array1d< real32 > m_a_np1;

  //Coefficient of the solution in the POD basis at time n
  array1d< real32 > m_a_n;

  //Coefficient of the solution in the POD basis at time n-1
  array1d< real32 > m_a_nm1;

  //Contribution on POD basis for the source
  array2d< real64 > m_sourceConstantsPOD;

  //Contribution on POD basis for the receivers
  array2d< real64 > m_receiverConstantsPOD;

  // Right hand side in the POD basis
  array1d< real32 > m_rhsPOD;

  //Size of POD basis for each Frechet order
  array1d< int > m_sizePOD_f;

  //Size of final POD basis
  int m_sizePOD;

};


namespace fields
{

namespace acousticfields
{

DECLARE_FIELD( Pressure64_nm1,
	       "pressure_nm1",
               array1d< pFieldType >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar of pressure at time n-1." );

DECLARE_FIELD( Pressure64_n,
               "pressure_n",
	       array1d< pFieldType >,
               0,
               NOPLOT,
               WRITE_AND_READ,
	       "Scalar of pressure at time n." );

DECLARE_FIELD( Pressure64_np1,
	       "pressure_np1",
	       array1d< pFieldType >,
	       0,
               LEVEL_0,
	       WRITE_AND_READ,
	       "Scalar of pressure at time n+1." );

DECLARE_FIELD( StiffnessVector64,
               "stiffnessVector",
	       array1d< pFieldType >,
               0,
	       LEVEL_0,
               WRITE_AND_READ,
	       "Resulting vector of stiffness matrix pressure vector product at time n" );
  
DECLARE_FIELD( PressureFrechet_nm1,
               "pressureFrechet_nm1",
               array2d< pFieldType >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar Frechet derivative of pressure at time n-1." );

DECLARE_FIELD( PressureFrechet_n,
               "pressureFrechet_n",
               array2d< pFieldType >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar Frechet derivative of pressure at time n." );

DECLARE_FIELD( PressureFrechet_np1,
               "pressureFrechet_np1",
               array2d< pFieldType >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar Frechet derivative of pressure at time n+1." );

DECLARE_FIELD( AcousticMassPerturbationVector,
               "acousticMassPerturbationVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Mass Matrix with perturbation." );

DECLARE_FIELD( DampingPerturbationVector,
               "dampingPerturbationVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix with perturbation." );

DECLARE_FIELD( Perturbation,
               "perturbation",
               array1d< real32 >,
	       0,
               NOPLOT,
               WRITE_AND_READ,
	       "Perturbation to apply to the velocity field" );
  
DECLARE_FIELD( ForcingRHS_fp1,
               "rhs_fp1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS_fp1" );


}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICROMFRECHET_HPP_ */
