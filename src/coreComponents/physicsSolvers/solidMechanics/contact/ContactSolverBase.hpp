/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
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
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"

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

  /**
   * @brief Override setupSystem to ensure consistent f0/f1 face ordering for all
   *        FaceElementSubRegions before the DOF structure is built.
   *
   * Calls flipFaceMap + fixNeighboringFacesNormals on every FaceElementSubRegion,
   * then chains to PhysicsSolverBase::setupSystem.  This covers both pre-split meshes (called once
   * at startup) and SurfaceGenerator runs (called again whenever new fractures are created).
   */
  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  /**
   * @brief Compute the local rotation matrices (normal, t1, t2) for all fracture elements.
   *
   * Uses the shared GPU/CPU ComputeRotationMatricesKernel with m_tangentRefDirection to
   * build a deterministic, canonical tangent frame.  Assumes flipFaceMap has already
   * been called (e.g. via setupSystem) so that Nbar is consistently oriented.
   */
  void computeRotationMatrices( DomainPartition & domain ) const;

  string const & getUniqueFractureRegionName() const { return m_fractureRegionNames[0]; }

  void outputConfigurationStatistics( DomainPartition const & domain ) const override final;

  void synchronizeFractureState( DomainPartition & domain ) const;

  struct viewKeyStruct : SolidMechanicsLagrangianFEM::viewKeyStruct
  {
    constexpr static char const * fractureStateString() { return "fractureState"; }

    constexpr static char const * oldFractureStateString() { return "oldFractureState"; }

    constexpr static char const * frictionLawNameString() { return "frictionLawName"; }

    constexpr static char const * tangentRefDirectionString() { return "tangentReferenceDirection"; }
  };

protected:
  virtual void postInputInitialization() override;

  /// Reference direction for the canonical tangent frame (t1 = normalize(refDir x normal), t2 = normal x t1).
  R1Tensor m_tangentRefDirection;

  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override final;

  void computeFractureStateStatistics( MeshLevel const & mesh,
                                       globalIndex & numStick,
                                       globalIndex & numNewSlip,
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

  void setFractureRegions( dataRepository::Group const & domain );

  stdVector< string > m_fractureRegionNames;

  template< typename LAMBDA >
  void forFractureRegionOnMeshTargets( Group const & meshBodies, LAMBDA && lambda ) const
  {
    forDiscretizationOnMeshTargets( meshBodies,
                                    [&]( string const,
                                         MeshLevel const & mesh,
                                         string_array const & )
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
                                         string_array const & )
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
