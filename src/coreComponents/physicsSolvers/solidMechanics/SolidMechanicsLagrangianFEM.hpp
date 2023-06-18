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
 * @file SolidMechanicsLagrangianFEM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "common/TimingMacros.hpp"
#include "kernels/SolidMechanicsLagrangianFEMKernels.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geosx
{

/**
 * @class SolidMechanicsLagrangianFEM
 *
 * This class implements a finite element solution to the equations of motion.
 */
class SolidMechanicsLagrangianFEM : public SolverBase
{
public:

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "solid"; }

  /**
   * @enum TimeIntegrationOption
   *
   * The options for time integration
   */
  enum class TimeIntegrationOption : integer
  {
    QuasiStatic,      //!< QuasiStatic
    ImplicitDynamic,  //!< ImplicitDynamic
    ExplicitDynamic   //!< ExplicitDynamic
  };

  /**
   * Constructor
   * @param name The name of the solver instance
   * @param parent the parent group of the solver
   */
  SolidMechanicsLagrangianFEM( const string & name,
                               Group * const parent );


  SolidMechanicsLagrangianFEM( SolidMechanicsLagrangianFEM const & ) = delete;
  SolidMechanicsLagrangianFEM( SolidMechanicsLagrangianFEM && ) = default;

  SolidMechanicsLagrangianFEM & operator=( SolidMechanicsLagrangianFEM const & ) = delete;
  SolidMechanicsLagrangianFEM & operator=( SolidMechanicsLagrangianFEM && ) = delete;

  /**
   * destructor
   */
  virtual ~SolidMechanicsLagrangianFEM() override;

  /**
   * @return The string that may be used to generate a new instance from the SolverBase::CatalogInterface::CatalogType
   */
  static string catalogName() { return "SolidMechanics_LagrangianFEM"; }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override final;

  void updateIntrinsicNodalData( DomainPartition * const domain );


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

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = false ) override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override final
  {
    // There should be nothing to update
    GEOSX_UNUSED_VAR( domain );
  };

  virtual void applyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void implicitStepComplete( real64 const & time,
                                     real64 const & dt,
                                     DomainPartition & domain ) override;

  /**@}*/


  template< typename CONSTITUTIVE_BASE,
            typename KERNEL_WRAPPER,
            typename ... PARAMS >
  void assemblyLaunch( DomainPartition & domain,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs,
                       PARAMS && ... params );


  template< typename ... PARAMS >
  real64 explicitKernelDispatch( MeshLevel & mesh,
                                 arrayView1d< string const > const & targetRegions,
                                 string const & finiteElementName,
                                 real64 const dt,
                                 std::string const & elementListName );

  /**
   * Applies displacement boundary conditions to the system for implicit time integration
   * @param time The time to use for any lookups associated with this BC
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain The DomainPartition.
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   */
  void applyDisplacementBCImplicit( real64 const time,
                                    DofManager const & dofManager,
                                    DomainPartition & domain,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs );

  void applyTractionBC( real64 const time,
                        DofManager const & dofManager,
                        DomainPartition & domain,
                        arrayView1d< real64 > const & localRhs );

  void applyChomboPressure( DofManager const & dofManager,
                            DomainPartition & domain,
                            arrayView1d< real64 > const & localRhs );


  void applyContactConstraint( DofManager const & dofManager,
                               DomainPartition & domain,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs );

  virtual real64
  scalingForSystemSolution( DomainPartition const & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr char const * cflFactorString() { return "cflFactor"; }
    static constexpr char const * newmarkGammaString() { return "newmarkGamma"; }
    static constexpr char const * newmarkBetaString() { return "newmarkBeta"; }
    static constexpr char const * massDampingString() { return "massDamping"; }
    static constexpr char const * stiffnessDampingString() { return "stiffnessDamping"; }
    static constexpr char const * timeIntegrationOptionString() { return "timeIntegrationOption"; }
    static constexpr char const * maxNumResolvesString() { return "maxNumResolves"; }
    static constexpr char const * strainTheoryString() { return "strainTheory"; }
    static constexpr char const * solidMaterialNamesString() { return "solidMaterialNames"; }
    static constexpr char const * contactRelationNameString() { return "contactRelationName"; }
    static constexpr char const * noContactRelationNameString() { return "NOCONTACT"; }
    static constexpr char const * maxForceString() { return "maxForce"; }
    static constexpr char const * elemsAttachedToSendOrReceiveNodesString() { return "elemsAttachedToSendOrReceiveNodes"; }
    static constexpr char const * elemsNotAttachedToSendOrReceiveNodesString() { return "elemsNotAttachedToSendOrReceiveNodes"; }
    static constexpr char const * sendOrReceiveNodesString() { return "sendOrReceiveNodes";}
    static constexpr char const * nonSendOrReceiveNodesString() { return "nonSendOrReceiveNodes";}
    static constexpr char const * targetNodesString() { return "targetNodes";}
    static constexpr char const * forceString() { return "Force";}

    dataRepository::ViewKey newmarkGamma = { newmarkGammaString() };
    dataRepository::ViewKey newmarkBeta = { newmarkBetaString() };
    dataRepository::ViewKey massDamping = { massDampingString() };
    dataRepository::ViewKey stiffnessDamping = { stiffnessDampingString() };
    dataRepository::ViewKey timeIntegrationOption = { timeIntegrationOptionString() };
  } solidMechanicsViewKeys;

  SortedArray< localIndex > & getElemsAttachedToSendOrReceiveNodes( ElementSubRegionBase & subRegion )
  {
    return subRegion.getReference< SortedArray< localIndex > >( viewKeyStruct::elemsAttachedToSendOrReceiveNodesString() );
  }

  SortedArray< localIndex > & getElemsNotAttachedToSendOrReceiveNodes( ElementSubRegionBase & subRegion )
  {
    return subRegion.getReference< SortedArray< localIndex > >( viewKeyStruct::elemsNotAttachedToSendOrReceiveNodesString() );
  }

  real64 & getMaxForce() { return m_maxForce; }

  arrayView1d< ParallelVector > const & getRigidBodyModes() const
  {
    return m_rigidBodyModes;
  }

  array1d< ParallelVector > & getRigidBodyModes()
  {
    return m_rigidBodyModes;
  }

  /**
   * Applies displacement boundary conditions that are passed between solvers - applicable in some multi-level cases
   * @param time The time to use for any lookups associated with this BC
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain The DomainPartition.
   * @param matrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  void applyInternalDisplacementBCImplicit( real64 const time,
                                            DofManager const & dofManager,
                                            DomainPartition & domain,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs );

  /**
   * Given a set of nodes and a list of values, prepare internal boundary conditions to have u_fixedNodes = fixedValues
   * @param fixedNodes The nodes which will have the displacement prescribed as a BC
   * @param fixedValues The values of displacement to be prescribe at these nodes
   */
  void setInternalBoundaryConditions( arrayView1d< localIndex > const & fixedNodes, arrayView2d< real64 > const & fixedValues, map<globalIndex, globalIndex>* patchToBaseMap );

  void setPressureEffects()
  {
    m_pressureEffectsFlag = 1;
  }

  void setSubdomainElements( SortedArray< localIndex > const & subdomainElems )
  {
    m_subdomainElems = subdomainElems.toView();
  }  

protected:
  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override;

  real64 m_newmarkGamma;
  real64 m_newmarkBeta;
  real64 m_massDamping;
  real64 m_stiffnessDamping;
  TimeIntegrationOption m_timeIntegrationOption;
  real64 m_maxForce = 0.0;
  integer m_maxNumResolves;
  integer m_strainTheory;
  integer m_pressureEffectsFlag;
  string m_contactRelationName;
  MPI_iCommData m_iComm;
  array1d< localIndex > m_fixedDisplacementNodes;
  array2d< real64 > m_fixedDisplacementValues;
  bool m_internalBCsFlag;
  SortedArrayView< const localIndex > m_subdomainElems;
  map<globalIndex, globalIndex>* m_patchToBaseElementRelation;



  /// Rigid body modes
  array1d< ParallelVector > m_rigidBodyModes;

private:
  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

};

ENUM_STRINGS( SolidMechanicsLagrangianFEM::TimeIntegrationOption,
              "QuasiStatic",
              "ImplicitDynamic",
              "ExplicitDynamic" );

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************


template< typename CONSTITUTIVE_BASE,
          typename KERNEL_WRAPPER,
          typename ... PARAMS >
void SolidMechanicsLagrangianFEM::assemblyLaunch( DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs,
                                                  PARAMS && ... params )
{
  GEOSX_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();

    string const dofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
    arrayView1d< globalIndex const > const & dofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    KERNEL_WRAPPER kernelWrapper( dofNumber,
                                  dofManager.rankOffset(),
                                  localMatrix,
                                  localRhs,
                                  gravityVectorData,
                                  std::forward< PARAMS >( params )... );

    m_maxForce = finiteElement::
                   regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                                 CONSTITUTIVE_BASE,
                                                 CellElementSubRegion >( mesh,
                                                                         regionNames,
                                                                         this->getDiscretizationName(),
                                                                         viewKeyStruct::solidMaterialNamesString(),
                                                                         kernelWrapper );
  } );


  applyContactConstraint( dofManager, domain, localMatrix, localRhs );
}

template< typename FE_TYPE >
struct BaseDispInterpolationKernel
{
  //this is the patch subRegion
  BaseDispInterpolationKernel( CellElementSubRegion const & subRegion, map<globalIndex, globalIndex>* patchToBaseMap ):
    m_numElems( subRegion.size() ),
    m_patchToBaseMap( patchToBaseMap )
  {}

  void interpolateBaseDisp( DomainPartition & domain )
  {

    //base node manager
    NodeManager & baseNodeManager = domain.getMeshBody( 0 ).getBaseDiscretization().getNodeManager();    
    NodeManager & patchNodeManager = domain.getMeshBody( 1 ).getBaseDiscretization().getNodeManager();    
    ElementRegionBase const & baseElementRegion = domain.getMeshBody( 0 ).getBaseDiscretization().getElemManager().getRegion( 0 );
    ElementRegionBase const & patchElementRegion = domain.getMeshBody( 1 ).getBaseDiscretization().getElemManager().getRegion( 0 );
    CellElementSubRegion const & baseSubRegion = baseElementRegion.getSubRegion< CellElementSubRegion >( 0 );
    CellElementSubRegion const & patchSubRegion = patchElementRegion.getSubRegion< CellElementSubRegion >( 0 );
    arrayView1d< integer const > const & patchGhostRank = patchSubRegion.ghostRank();
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const xNodesBase = baseNodeManager.referencePosition();
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const xNodesPatch = patchNodeManager.referencePosition();
    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & baseNodalDisp = baseNodeManager.getField< fields::solidMechanics::totalDisplacement >();    
    auto const & baseElemsToNodes  = baseSubRegion.nodeList().toViewConst();
    auto const & patchElemsToNodes = patchSubRegion.nodeList().toViewConst();
    arrayView1d< real64 > nodalDispXOnMesh = patchNodeManager.getReference< array1d< real64 > >( "TranfsDispX" );
    arrayView1d< real64 > nodalDispYOnMesh = patchNodeManager.getReference< array1d< real64 > >( "TranfsDispY" );
    arrayView1d< real64 > nodalDispZOnMesh = patchNodeManager.getReference< array1d< real64 > >( "TranfsDispZ" );    
    arrayView2d< const real64 > baseElemCenters = baseSubRegion.getElementCenter();
    arrayView1d< globalIndex const > patchLocalToGlobalMap = patchSubRegion.localToGlobalMap();
    unordered_map< globalIndex, localIndex > const & baseGlobalToLocalMap = baseSubRegion.globalToLocalMap();

    //looping over all patchElems
    forAll< serialPolicy >( m_numElems, [=] ( localIndex const k )
    {
      //useful element properties
      constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
      
      //only interpolate over real elements because the mapping function doesnt work well for the others
      if(patchGhostRank[k]<0)
      {
        
        real64 xLocalPatch[ numNodesPerElement ][ 3 ];
        real64 baseNodalDispLocal[ numNodesPerElement ][3];

        //conversion from patch to base
        globalIndex K = patchLocalToGlobalMap[k];
        globalIndex baseK = (*m_patchToBaseMap)[K];
        localIndex basek = baseGlobalToLocalMap.at(baseK);
        //unordered_map< globalIndex, localIndex > const & 
        //arrayView1d< globalIndex const > 
        //these element measures are half of the element size
        real64 hx = abs(xNodesBase[baseElemsToNodes( basek, 0 )][0] - baseElemCenters[basek][0]);
        real64 hy = abs(xNodesBase[baseElemsToNodes( basek, 0 )][1] - baseElemCenters[basek][1]);
        real64 hz = abs(xNodesBase[baseElemsToNodes( basek, 0 )][2] - baseElemCenters[basek][2]);

        for( localIndex a = 0; a < numNodesPerElement; ++a )
        {
          localIndex const localNodeIndex = patchElemsToNodes( k, a );

          for( int dim=0; dim < 3; ++dim )
          {
            xLocalPatch[a][dim] = xNodesPatch[ localNodeIndex ][dim];
          }

          //element level displacement values at nodes
          baseNodalDispLocal[ a ][0] = baseNodalDisp[ localNodeIndex ][0];
          baseNodalDispLocal[ a ][1] = baseNodalDisp[ localNodeIndex ][1];
          baseNodalDispLocal[ a ][2] = baseNodalDisp[ localNodeIndex ][2];

          real64 N[ numNodesPerElement ];

          //convert x,y,z to parent coordinates
          real64 parentX[3];
          parentX[0] = (xLocalPatch[a][0] - baseElemCenters[basek][0]) / hx;
          parentX[1] = (xLocalPatch[a][1] - baseElemCenters[basek][1]) / hy;
          parentX[2] = (xLocalPatch[a][2] - baseElemCenters[basek][2]) / hz;

          //get shape function values at parent coordinates
          FE_TYPE::calcN( parentX, N );

          real64 nDisp[3];
          //interpolate disp field using shape function values
          FE_TYPE::value( N, baseNodalDispLocal, nDisp );

          //write to mesh
          nodalDispXOnMesh[localNodeIndex] = nDisp[0];
          nodalDispYOnMesh[localNodeIndex] = nDisp[1];
          nodalDispZOnMesh[localNodeIndex] = nDisp[2];
        
        }

      }
    } );
  }

  localIndex m_numElems;
  map<globalIndex, globalIndex>* m_patchToBaseMap;
};     


} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_ */
