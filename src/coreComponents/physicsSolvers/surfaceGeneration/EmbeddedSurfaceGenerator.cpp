/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EmbeddedSurfaceGenerator.cpp
 */

#include "EmbeddedSurfaceGenerator.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "mpiCommunications/SpatialPartition.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEMKernels.hpp"
#include "meshUtilities/SimpleGeometricObjects/GeometricObjectManager.hpp"
#include "meshUtilities/SimpleGeometricObjects/BoundedThickPlane.hpp"

#ifdef USE_GEOSX_PTP
#include "GEOSX_PTP/ParallelTopologyChange.hpp"
#endif

#include <set>

namespace geosx
{
  using namespace dataRepository;
  using namespace constitutive;

EmbeddedSurfaceGenerator::EmbeddedSurfaceGenerator( const std::string& name,
                                    Group * const parent ):
  SolverBase( name, parent ),
  m_solidMaterialName("")
{
  registerWrapper(viewKeyStruct::solidMaterialNameString, &m_solidMaterialName, 0)->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Name of the solid material used in solid mechanic solver");

  registerWrapper( viewKeyStruct::fractureRegionNameString, &m_fractureRegionName, 0 )->
      setInputFlag(dataRepository::InputFlags::OPTIONAL)->
      setApplyDefaultValue("FractureRegion");
}

EmbeddedSurfaceGenerator::~EmbeddedSurfaceGenerator()
{
  // TODO Auto-generated destructor stub
}

void EmbeddedSurfaceGenerator::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0);

    NodeManager * const nodeManager = meshLevel->getNodeManager();
    EdgeManager * const edgeManager = meshLevel->getEdgeManager();

    nodeManager->registerWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::parentIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Parent index of node.");

    nodeManager->registerWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::childIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Child index of node.");

    edgeManager->registerWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::parentIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Parent index of the edge.");

    edgeManager->registerWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::childIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Child index of the edge.");
  }
}

void EmbeddedSurfaceGenerator::InitializePostSubGroups( Group * const problemManager )
{
  /*
   * Here we generate embedded elements for fractures (or faults) that already exist in the domain and
   * were specified in the input file.
   */

  // Get domain
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>( dataRepository::keys::domain );
  // Get geometric object manager
  GeometricObjectManager * geometricObjManager = problemManager->GetGroup<GeometricObjectManager>( "Geometry");

  // Get meshLevel
  Group     * const meshBodies = domain->getMeshBodies();
  MeshBody  * const meshBody   = meshBodies->GetGroup<MeshBody>(0);
  MeshLevel * const meshLevel  = meshBody->GetGroup<MeshLevel>(0);

  // Get managers
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  //NodeManager * const nodeManager = meshLevel->getNodeManager();
  // EdgeManager * const edgeManager = meshLevel->getEdgeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();
  array1d<R1Tensor> & faceCenter  = faceManager->faceCenter();

  // Get EmbeddedSurfaceSubRegions
  EmbeddedSurfaceRegion    * const    embeddedSurfaceRegion =
      elemManager->GetRegion<EmbeddedSurfaceRegion>(this->m_fractureRegionName);
  EmbeddedSurfaceSubRegion * const embeddedSurfaceSubRegion =
      embeddedSurfaceRegion->GetSubRegion<EmbeddedSurfaceSubRegion>(0);

  // Loop over all the fracture planes
  geometricObjManager->forSubGroups<ThickPlane>( [&]( ThickPlane * const fracture ) -> void
  {
    /* 1. Find out if an element is cut but the fracture or not.
     * Loop over all the elements and for each one of them loop over the faces and compute the
     * dot product between the distance between the plane center and the face and the normal
     * vector defining the plane. If two scalar products have different signs the plane cuts the
     * cell. To do this check multiply all dot products and then check the sign. If a face gives
     * a 0 dot product it has to be neglected or the method won't work (as a matter of fact it means
     * that the plane does cut the cell).
     *
     */
    R1Tensor planeCenter  = fracture->getCenter();
    R1Tensor normalVector = fracture->getNormal();
    // Initialize variables
    globalIndex faceIndex;
    real64 prodScalarProduct;
    R1Tensor distVec;

    elemManager->forElementRegions( [&](ElementRegionBase * const region )->void
    {
      Group * subRegions = region->GetGroup(ElementRegionBase::viewKeyStruct::elementSubRegions);
      subRegions->forSubGroups<CellElementSubRegion>( [&]( CellElementSubRegion * const subRegion ) -> void
      {
        FixedOneToManyRelation const & cellToFaces = subRegion->faceList();
        for(localIndex cellIndex =0; cellIndex<subRegion->size(); cellIndex++)
        {
          prodScalarProduct = 1;
          for(localIndex kf =0; kf<subRegion->numFacesPerElement(); kf++)
          {
            faceIndex = cellToFaces[cellIndex][kf];
            distVec  = faceCenter[faceIndex];
            distVec -= planeCenter;
            // check if the dot product is zero
            if ( std::fabs(Dot(distVec, normalVector)) > 1e-15 )
            {
              prodScalarProduct *= Dot(distVec, normalVector);
            }
          }
          if (prodScalarProduct < 0)
            // TODO in reality this condition is not sufficient because the fracture is a bounded plane. I should also check
            // that the cell is inside the actual fracture plane. For now let us assume the plane is not bounded.
          {
            /* 2. Now that you know that the element is cut by the fracture we can
             * a. add new embedded surface element
             * b. fill in the data relative to where the actual intersections are
             * c. compute the geometric quantities and the Heaviside.
            */
            // a. Add the embedded surface element
            embeddedSurfaceSubRegion->AddNewEmbeddedSurface(cellIndex, normalVector);
            // b. find actual intersections with each edge and compute area.
            // embeddedSurfaceSubRegion->ComputeElementArea();
          }
        }
      });
    });
  });

  std::cout << "Number of embedded surface elements: " << embeddedSurfaceSubRegion->size() << std::endl;
}

void EmbeddedSurfaceGenerator::InitializePostInitialConditions_PreSubGroups( Group * const  GEOSX_UNUSED_ARG ( problemManager ) )
{
  // I don't think there is  much to do here.
}


void EmbeddedSurfaceGenerator::postRestartInitialization( Group * const GEOSX_UNUSED_ARG( domain0 ) )
{
  // Not sure about this for now.
  std::cout << "postRestartInitialization \n";
}


real64 EmbeddedSurfaceGenerator::SolverStep( real64 const & GEOSX_UNUSED_ARG( time_n),
                                             real64 const & GEOSX_UNUSED_ARG( dt ),
                                             const int GEOSX_UNUSED_ARG( cycleNumber ),
                                             DomainPartition * const  GEOSX_UNUSED_ARG( domain ) )
{
  real64 rval = 0;
  /*
   * This should be method that generates new fracture elements based on the propagation criterion of choice.
   */
  return rval;
}


REGISTER_CATALOG_ENTRY( SolverBase,
                        EmbeddedSurfaceGenerator,
                        std::string const &, dataRepository::Group * const )

} /* namespace geosx */
