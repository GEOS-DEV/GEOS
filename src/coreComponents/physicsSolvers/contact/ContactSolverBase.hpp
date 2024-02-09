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

/*
 *  @file ContactSolverBase.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_CONTACTSOLVERBASE_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_CONTACTSOLVERBASE_HPP_

#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"

namespace geos
{

class ContactSolverBase : public SolidMechanicsLagrangianFEM
{
public:
  ContactSolverBase( const string & name,
                     Group * const parent );

  ~ContactSolverBase() override = default;

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override;

  virtual real64
  explicitStep( real64 const & time_n,
                real64 const & dt,
                integer const cycleNumber,
                DomainPartition & domain ) override final;

  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  string const & getContactRelationName() const { return m_contactRelationName; }

  string const & getFractureRegionName() const { return m_fractureRegionNames[0]; }

  void outputConfigurationStatistics( DomainPartition const & domain ) const override final;

  //void setSolidSolverDofFlags( bool const flag ) { m_setupSolidSolverDofs = flag; }

  void synchronizeFractureState( DomainPartition & domain ) const;

protected:

  //virtual void postProcessInput() override;

  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override final;

  void computeFractureStateStatistics( MeshLevel const & mesh,
                                       globalIndex & numStick,
                                       globalIndex & numSlip,
                                       globalIndex & numOpen ) const;

  GEOS_HOST_DEVICE
  inline
  static bool compareFractureStates( integer const state0,
                                     integer const state1 )
  {
    return state0 == state1
           || ( state0 == fields::contact::FractureState::NewSlip && state1 == fields::contact::FractureState::Slip )
           || ( state0 == fields::contact::FractureState::Slip && state1 == fields::contact::FractureState::NewSlip );
  }

  void setFractureRegions(  dataRepository::Group const & domain ); 

  /// contact relation name string
  string m_contactRelationName;

  std::vector< string > m_fractureRegionNames;

  ///
  //bool m_setupSolidSolverDofs;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * contactRelationNameString() { return "contactRelationName"; }

    constexpr static char const * fractureStateString() { return "fractureState"; }

    constexpr static char const * oldFractureStateString() { return "oldFractureState"; }

    constexpr static char const * initialFractureStateString() { return "initialFractureState"; }
  };

  template< typename LAMBDA >
  void forFractureRegionOnMeshTargets( Group const & meshBodies, LAMBDA && lambda ) const
  {
    forDiscretizationOnMeshTargets( meshBodies,
                                    [&]( string const,
                                         MeshLevel const & mesh,
                                         arrayView1d< string const > const )
    {
      ElementRegionManager const & elemManager = mesh.getElemManager();

      elemManager.forElementRegions< SurfaceElementRegion >( m_fractureRegionNames,
                                                             [&] ( localIndex const,
                                                                   SurfaceElementRegion const & region )
      {
        lambda( region );
      } );
    } );
  }

  template< typename LAMBDA >
  void forFractureRegionOnMeshTargets( Group & meshBodies, LAMBDA && lambda ) const
  {
    forDiscretizationOnMeshTargets( meshBodies,
                                    [&]( string const,
                                         MeshLevel & mesh,
                                         arrayView1d< string const > const )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();

      elemManager.forElementRegions< SurfaceElementRegion >( m_fractureRegionNames,
                                                             [&] ( localIndex const,
                                                                   SurfaceElementRegion & region )
      {
        lambda( region );
      } );
    } );
  }
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_CONTACTSOLVERBASE_HPP_ */
