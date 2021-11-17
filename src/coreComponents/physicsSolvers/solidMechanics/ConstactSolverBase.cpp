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
 * ContactSolverBase.cpp
 */

#include "ContactSolverBase.hpp"

#include "SolidMechanicsEFEMKernels.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactSelector.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "finiteElement/elementFormulations/FiniteElementBase.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "mesh/DomainPartition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

ContactSolverBase::ContactSolverBase( const string & name,
                                      Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_fractureRegionName(),
  m_solidSolver( nullptr )
{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver in the rock matrix" );

  registerWrapper( viewKeyStruct::fractureRegionNameString(), &m_fractureRegionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fracture region." );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );
}

ContactSolverBase::~ContactSolverBase()
{
  // TODO Auto-generated destructor stub
}

void ContactSolverBase::postProcessInput()
{
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
  SolverBase::postProcessInput();
}

void ContactSolverBase::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    ElementRegionManager & elemManager = meshLevel.getElemManager();
    {
      elemManager.forElementRegions< SurfaceElementRegion >( [&] ( SurfaceElementRegion & region )
      {
        region.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion & subRegion )
        {

          subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dispJumpString() ).
            setPlotLevel( PlotLevel::LEVEL_0 ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaDispJumpString() ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::oldDispJumpString() ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::fractureTractionString() ).
            reference().resizeDimension< 1 >( 3 );

          subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dTraction_dJumpString() ).
            reference().resizeDimension< 1, 2 >( 3, 3 );

          subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::fractureStateString() ).
            setPlotLevel( PlotLevel::LEVEL_0 ).
            setRegisteringObjects( this->getName()).
            setDescription( "An array that holds the fracture state." );
          initializeFractureState( meshLevel, viewKeyStruct::fractureStateString() );

          subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::oldFractureStateString() ).
            setPlotLevel( PlotLevel::NOPLOT ).
            setRegisteringObjects( this->getName()).
            setDescription( "An array that holds the fracture state." );
          initializeFractureState( meshLevel, viewKeyStruct::oldFractureStateString() );

        } );
      } );
    }
  } );
}

} /* namespace geosx */
