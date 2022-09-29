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
 * @file MultiphasePoroelasticSolver.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSSOLVER_HPP_

#include "constitutive/solid/CoupledSolidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"


namespace geosx
{

class MultiphasePoromechanicsSolver : public CoupledSolver< SolidMechanicsLagrangianFEM,
                                                            CompositionalMultiphaseBase >
{
public:

  using Base = CoupledSolver< SolidMechanicsLagrangianFEM, CompositionalMultiphaseBase >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    SolidMechanics = 0,
    Flow = 1
  };

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanics"; }

  /**
   * @brief main constructor for MultiphasePoromechanicsSolver Objects
   * @param name the name of this instantiation of MultiphasePoromechanicsSolver in the repository
   * @param parent the parent group of this instantiation of MultiphasePoromechanicsSolver
   */
  MultiphasePoromechanicsSolver( const string & name,
                                 Group * const parent );

  /// Destructor for the class
  ~MultiphasePoromechanicsSolver() override {};

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new MultiphasePoromechanicsSolver object through the object catalog.
   */
  static string catalogName() { return "MultiphasePoromechanics"; }

  /**
   * @brief accessor for the pointer to the solid mechanics solver
   * @return a pointer to the solid mechanics solver
   */
  SolidMechanicsLagrangianFEM * solidMechanicsSolver() const
  {
    return std::get< toUnderlying( SolverType::SolidMechanics ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  CompositionalMultiphaseBase * flowSolver() const
  {
    return std::get< toUnderlying( SolverType::Flow ) >( m_solvers );
  }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;


  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override;

  void updateStabilizationParams(bool updateMacro, bool updateTau)
  { 

    if( m_stabilizationType == StabilizationType::Global )
    {

       DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

       // Step 1: we loop over the regions where stabilization is active and collect their name
    
       set< string > regionFilter;
       for( string const & regionName : m_stabilizationRegionNames )
       {
         regionFilter.insert( regionName );
       }

       // Step 2: loop over the target regions of the solver, and tag the elements belonging to stabilization regions
       forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & targetRegionNames )
       {
         // keep only the target regions that are in the filter
         array1d< string > filteredTargetRegionNames;
         filteredTargetRegionNames.reserve( targetRegionNames.size() );

         for( string const & targetRegionName : targetRegionNames )
         {
           if( regionFilter.count( targetRegionName ) )
           {
             filteredTargetRegionNames.emplace_back( targetRegionName );
           }
         }

         // loop over the elements and make stabilization active
         mesh.getElemManager().forElementSubRegions( filteredTargetRegionNames.toViewConst(), [&]( localIndex const,
                                                                                                ElementSubRegionBase & subRegion )

         {
           arrayView1d< integer > const macroElementIndex = subRegion.getExtrinsicData< extrinsicMeshData::flow::macroElementIndex >();
           arrayView1d< real64 > const elementStabConstant = subRegion.getExtrinsicData< extrinsicMeshData::flow::elementStabConstant >();
        
           geosx::constitutive::CoupledSolidBase const & porousSolid =
           getConstitutiveModel< geosx::constitutive::CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::porousMaterialNamesString() ) );

           arrayView1d< const real64 > const bulkModulus = porousSolid.getBulkModulus();
           arrayView1d< const real64 > const shearModulus = porousSolid.getShearModulus();
           arrayView1d< const real64 > const biotCoefficient = porousSolid.getBiotCoefficient();

           forAll< parallelHostPolicy >( subRegion.size(), [&] ( localIndex const ei )
           {
             if (updateMacro)
             {
               macroElementIndex[ei] = 1;
             }

             if (updateTau)
             {
               real64 bM = bulkModulus[ei];
               real64 sM = shearModulus[ei];
               real64 bC = biotCoefficient[ei];
               elementStabConstant[ei] = m_stabilizationMultiplier * 9.0 * (bC * bC) / (32.0 * (10.0 * sM / 3.0 + bM));
             }
           } );
         } );
 
       } );
    
    } 
  }

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override
  {
    CoupledSolver< SolidMechanicsLagrangianFEM, CompositionalMultiphaseBase > :: implicitStepComplete(time_n, dt, domain);
    updateStabilizationParams(false, true);
  }

  /**@}*/

  enum class StabilizationType : integer
  {
    None,
    Global,
    Local,
  };

protected:

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    /// Names of the porous materials
    constexpr static char const * porousMaterialNamesString() { return "porousMaterialNames"; }

    /// Type of stabilization used in the simulation
    constexpr static char const * stabilizationTypeString() { return "stabilizationType"; }

    /// Names of the regions where the stabilization is applied
    constexpr static char const * stabilizationRegionNamesString() { return "stabilizationRegionNames"; }

    /// Multiplier on stabilization
    constexpr static char const * stabilizationMultiplierString() { return "stabilizationMultiplier"; }
  };

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  /// Type of stabilization used in the simulation
  StabilizationType m_stabilizationType;

  /// Names of the regions where the stabilization is applied
  array1d< string > m_stabilizationRegionNames;

  /// Multiplier on stabilization constant
  real64 m_stabilizationMultiplier;

};

ENUM_STRINGS( MultiphasePoromechanicsSolver::StabilizationType,
              "None",
              "Global",
              "Local" );

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSSOLVER_HPP_ */
