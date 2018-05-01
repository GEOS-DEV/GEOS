/*
 * DomainPartition.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#include "DomainPartition.hpp"

#include "../MPI_Communications/SpatialPartition.hpp"
#include "constitutive/ConstitutiveManager.hpp"

#include "fileIO/silo/SiloFile.hpp"

#include "common/SortedArray.hpp"

#include "common/Logger.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "managers/ObjectManagerBase.hpp"
namespace geosx
{
using namespace dataRepository;

DomainPartition::DomainPartition( std::string const & name,
                                  ManagedGroup * const parent ):
  ManagedGroup( name, parent ),
  m_mpiComm()
{


  this->RegisterViewWrapper< array<NeighborCommunicator> >(viewKeys.neighbors);
  MPI_Comm_dup( MPI_COMM_WORLD, &m_mpiComm );
  this->RegisterViewWrapper<SpatialPartition,PartitionBase>(keys::partitionManager)->setRestartFlags( RestartFlags::NO_WRITE );

  RegisterGroup( groupKeys.meshBodies );
  RegisterGroup<constitutive::ConstitutiveManager>( groupKeys.constitutiveManager );
  RegisterGroup<CellBlockManager>( keys::cellManager );
}

DomainPartition::~DomainPartition()
{}



void DomainPartition::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("Domain");
  docNode->setSchemaType("UniqueNode");
}


void DomainPartition::InitializationOrder( string_array & order )
{
  set<string> usedNames;
  {
    order.push_back(keys::ConstitutiveManager);
    usedNames.insert(keys::ConstitutiveManager);
  }

  {
    order.push_back(groupKeysStruct::meshBodiesString);
    usedNames.insert(groupKeysStruct::meshBodiesString);
  }


  for( auto const & subGroup : this->GetSubGroups() )
  {
    if( usedNames.count(subGroup.first) == 0 )
    {
      order.push_back(subGroup.first);
    }
  }
}


void DomainPartition::SetMaps(  )
{
  // ManagedGroup * nodeManager = this->GetGroup(keys::FEM_Nodes);
  // ElementRegionManager * elementRegionManager =
  // this->GetGroup<ElementRegionManager>(keys::FEM_Elements);

  // {
  //  integer_array & elementRegionMap =
  // nodeManager->getReference(keys::elementRegionMap);
  //  integer_array & elementSubRegionMap =
  // nodeManager->getReference(keys::elementSubRegionMap);
  //  integer_array & elementMap = nodeManager->getReference(keys::elementMap);

  //  ManagedGroup * elementRegions =
  // this->GetGroup(dataRepository::keys::elementRegions);

  //  integer elementRegionIndex = 0;
  //  elementRegionManager->forElementRegions( [&](ElementRegion&
  // elementRegion)-> void
  //  {
  //    elementRegion.forCellBlocks( [&]( CellBlockSubRegion & subRegion )->void
  //    {

  //    });
  //    ++elementRegionIndex;
  //  });
  // }
}

void DomainPartition::GenerateSets(  )
{
  MeshLevel * const mesh = this->getMeshBody(0)->getMeshLevel(0);
  ManagedGroup * nodeManager = mesh->getNodeManager();

  dataRepository::ManagedGroup const * nodeSets = nodeManager->GetGroup(dataRepository::keys::sets);

  std::map< string, integer_array > nodeInSet;
  string_array setNames;

  for( auto & viewWrapper : nodeSets->wrappers() )
  {
    string name = viewWrapper.second->getName();
    nodeInSet[name].resize( nodeManager->size() );
    nodeInSet[name] = 0;
    ViewWrapper<lSet> const * const setPtr = nodeSets->getWrapper<lSet>(name);
    if( setPtr!=nullptr )
    {
      setNames.push_back(name);
      lSet const & set = setPtr->reference();
      for( auto const a : set )
      {
        nodeInSet[name][a] = 1;
      }
    }
  }


  ElementRegionManager * elementRegionManager = mesh->getElemManager();

  for( auto & subGroup : elementRegionManager->GetGroup( dataRepository::keys::elementRegions )->GetSubGroups() )
//  elementRegionManager->forElementRegions( [&]( ElementRegion * elementRegion
// )
  {
    ElementRegion * elementRegion = subGroup.second->group_cast<ElementRegion *>();
//    elementRegion->forCellBlocks( [&]( CellBlockSubRegion * subRegion )->void
    for( auto & subRegionIter : elementRegion->GetGroup(dataRepository::keys::cellBlockSubRegions)->GetSubGroups() )
    {
      CellBlockSubRegion * subRegion = subRegionIter.second->group_cast<CellBlockSubRegion *>();
      lArray2d const & elemsToNodes = subRegion->getWrapper<FixedOneToManyRelation>(subRegion->viewKeys.nodeList)->reference();// getData<lArray2d>(keys::nodeList);
      dataRepository::ManagedGroup * elementSets = subRegion->GetGroup(dataRepository::keys::sets);
      std::map< string, integer_array > numNodesInSet;

      for( auto & setName : setNames )
      {

        lSet & set = elementSets->RegisterViewWrapper<lSet>(setName)->reference();
        for( localIndex k = 0 ; k < subRegion->size() ; ++k )
        {
          arrayView1d<localIndex const> const nodelist = elemsToNodes[k];
          integer count = 0;
          for( localIndex a = 0 ; a<elemsToNodes.size(1) ; ++a )
          {
            if( nodeInSet[setName][nodelist[a]] == 1 )
            {
              ++count;
            }
          }
          if( count == elemsToNodes.size(1) )
          {
            set.insert(k);
          }
        }
      }
    }
  }
}


void DomainPartition::SetupCommunications()
{
  PartitionBase   & partition1 = getReference<PartitionBase>(keys::partitionManager);
  SpatialPartition & partition = dynamic_cast<SpatialPartition &>(partition1);
  array<NeighborCommunicator> & allNeighbors = this->getReference< array<NeighborCommunicator> >( viewKeys.neighbors );

  //get communicator, rank, and coordinates
  MPI_Comm cartcomm;
  {
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 3, partition.m_Partitions.data(), partition.m_Periodic.data(), reorder, &cartcomm);
  }
  int rank = -1;
  int nsdof = 3;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  MPI_Comm_rank(cartcomm, &rank);
  MPI_Cart_coords(cartcomm, rank, nsdof, partition.m_coords.data());

  int ncoords[3];
  AddNeighbors(0, cartcomm, ncoords);

  MPI_Comm_free(&cartcomm);


  ManagedGroup * const meshBodies = getMeshBodies();
  MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(0);
  MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(0);

  for( auto const & neighbor : allNeighbors )
  {
    neighbor.AddNeighborGroupToMesh(meshLevel);
  }



  NodeManager * nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();

  CommunicationTools::AssignGlobalIndices( *faceManager, *nodeManager, allNeighbors );

  CommunicationTools::FindMatchedPartitionBoundaryObjects( faceManager,
                                                           allNeighbors );

  CommunicationTools::FindMatchedPartitionBoundaryObjects( nodeManager,
                                                           allNeighbors );

  CommunicationTools::FindGhosts( meshLevel, allNeighbors );

}

void DomainPartition::AddNeighbors(const unsigned int idim,
                                   MPI_Comm& cartcomm,
                                   int* ncoords)
{
  PartitionBase   & partition1 = getReference<PartitionBase>(keys::partitionManager);
  SpatialPartition & partition = dynamic_cast<SpatialPartition &>(partition1);
  array<NeighborCommunicator> & allNeighbors = this->getReference< array<NeighborCommunicator> >( viewKeys.neighbors );

  if (idim == nsdof)
  {
    bool me = true;
    for ( unsigned int i = 0 ; i < nsdof ; i++)
    {
      if (ncoords[i] != partition.m_coords(i))
      {
        me = false;
        break;
      }
    }
    if (!me)
    {
      allNeighbors.push_back(NeighborCommunicator());
      int neighborRank;
      MPI_Cart_rank(cartcomm, ncoords, &neighborRank);
      allNeighbors.back().SetNeighborRank(neighborRank);
    }
  }
  else
  {
    const int dim = partition.m_Partitions(idim);
    const bool periodic = partition.m_Periodic(idim);
    for (int i = -1 ; i < 2 ; i++)
    {
      ncoords[idim] = partition.m_coords(idim) + i;
      bool ok = true;
      if (periodic)
      {
        if (ncoords[idim] < 0)
          ncoords[idim] = dim - 1;
        else if (ncoords[idim] >= dim)
          ncoords[idim] = 0;
      }
      else
      {
        ok = ncoords[idim] >= 0 && ncoords[idim] < dim;
      }
      if (ok)
      {
        AddNeighbors(idim + 1, cartcomm, ncoords);
      }
    }
  }
}
/**
 * @brief Write to SILO
 * @author R Settgast
 * Write all objects to SILO format
 * @param[out] siloFile SILO file object
 * @param[in] cycleNum Timestep index
 * @param[in] problemTime Current simulation time
 * @param[in] writeFEMMesh Flag whether to write out the finite element mesh
 * @param[in] writeFEMFaces Flag whether to write out the finite element faces
 * @param[in] writeFEMEdges Flag whether to write out the finite element edges
 * @param[in] writeDE Flag whether to write out the discrete elements
 * @param[in] writeCP Flag whether to write out the common plane contacts
 * @param[in] writeCG Flag whether to write out the cartesian grid
 */
void DomainPartition::WriteSilo(SiloFile& siloFile,
                                const int cycleNum,
                                const realT problemTime,
                                const bool isRestart )
{

  WriteFiniteElementMesh( siloFile, cycleNum, problemTime, isRestart );


//  WriteCommonPlanes( siloFile, cycleNum, problemTime, isRestart, writeCP );

//  WriteCartesianGrid( siloFile, cycleNum, problemTime, isRestart, writeCG );
//  m_wellboreManager.WriteWellboreSilo( siloFile, cycleNum, problemTime,
// isRestart );

  if( isRestart )
  {
//    siloFile.DBWriteWrapper("m_globalDomainNumber",m_globalDomainNumber);
  }

}

void DomainPartition::ReadSilo( const SiloFile& siloFile,
                                const int cycleNum,
                                const realT problemTime,
                                const bool isRestart )
{

  ReadFiniteElementMesh( siloFile, cycleNum, problemTime, isRestart );

//  ReadCommonPlanes( siloFile, cycleNum, problemTime, isRestart );
//  ReadCartesianGrid( siloFile, cycleNum, problemTime, isRestart );
//  m_wellboreManager.ReadSilo( siloFile, "WellboreFields", "wellbore_mesh",
//                              DB_NODECENT, cycleNum, problemTime, isRestart );

}
void DomainPartition::WriteFiniteElementMesh( SiloFile& siloFile,
                                              const int cycleNum,
                                              const realT problemTime,
                                              const bool isRestart )
{
  int rank = 0;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  //--------------WRITE FE DATA-----------------
//  if (m_feElementManager->m_numElems > 0)
  {

    MeshLevel const * const mesh = this->getMeshBody(0)->getMeshLevel(0);
    NodeManager const * const nodeManager = mesh->getNodeManager();

//    NodeManager const * nodeManager =
// this->GetGroup<NodeManager>(keys::FEM_Nodes);
    localIndex numNodes = nodeManager->size();

    r1_array const & referencePosition = nodeManager->getReference<r1_array>(keys::referencePositionString);

//    r1_array const * const displacement = nodeManager->GetFieldDataPointer<r1_array>(keys::TotalDisplacement);

    bool writeArbitraryPolygon(false);
    const std::string meshName("volume_mesh");
    //set the nodal coordinate data structure
    realT* coords[3];
    dvector xcoords(numNodes);
    dvector ycoords(numNodes);
    dvector zcoords(numNodes);
    for (localIndex a = 0 ; a < numNodes ; ++a)
    {
      R1Tensor nodePosition;
      nodePosition = referencePosition[a];
//      if( displacement!=nullptr )
//      {
//        nodePosition += (*displacement)[a];
//      }

      xcoords[a] = nodePosition(0);
      ycoords[a] = nodePosition(1);
      zcoords[a] = nodePosition(2);
    }

    coords[0] = xcoords.data();
    coords[1] = ycoords.data();
    coords[2] = zcoords.data();

    ElementRegionManager const * const elementManager = mesh->getElemManager();
    const localIndex numElementRegions = elementManager->GetGroup(keys::elementRegions)->GetSubGroups().size();
    array<localIndex*> meshConnectivity(numElementRegions);
    array<int*> isGhostElement(numElementRegions);
    array<globalIndex*> globalElementNumbers(numElementRegions);
    ivector shapecnt(numElementRegions);
    ivector shapetype(numElementRegions);
    ivector shapesize(numElementRegions);

    array<FixedOneToManyRelation> elementToNodeMap;
    elementToNodeMap.resize( numElementRegions );

    int count = 0;
//    elementManager->forCellBlocks([&]( CellBlockSubRegion const * cellBlock ) -> void
    ManagedGroup const * elementRegions = elementManager->GetGroup(dataRepository::keys::elementRegions);

    for( auto const & region : elementRegions->GetSubGroups() )
    {
      ManagedGroup const * cellBlockSubRegions = region.second->GetGroup(dataRepository::keys::cellBlockSubRegions);
      for( auto const & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
      {
        CellBlockSubRegion const * cellBlock = cellBlockSubRegions->GetGroup<CellBlockSubRegion>(iterCellBlocks.first);

    {
        lArray2d const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys.nodeList)->reference();// getData<lArray2d>(keys::nodeList);

        // The following line seems to be redundant. It's actual function is to
        // size this temp array.(pfu)
        elementToNodeMap[count].resize(elemsToNodes.size(0),elemsToNodes.size(1));

        for (localIndex k = 0 ; k < cellBlock->size() ; ++k)
        {
          arrayView1d<localIndex const> const elemToNodeMap = elemsToNodes[k];

          const integer_array nodeOrdering = siloFile.SiloNodeOrdering();
          integer numNodesPerElement = elemsToNodes.size(1);
          for (localIndex a = 0 ; a < numNodesPerElement ; ++a)
          {
            elementToNodeMap[count](k, a) = elemToNodeMap[nodeOrdering[a]];
          }

        }


        //      meshConnectivity[count] =
        // elementRegion.m_ElementToNodeMap.data();
        meshConnectivity[count] = elementToNodeMap[count].data();

//      isGhostElement[count] =
// (cellBlock->GetFieldData<FieldInfo::ghostRank>()).data();

//      globalElementNumbers[count] = elementRegion.m_localToGlobalMap.data();
        shapecnt[count] = static_cast<int>(cellBlock->size());

//      if ( !elementRegion.m_elementGeometryID.compare(0, 4, "C3D8") )
//      {
        shapetype[count] = DB_ZONETYPE_HEX;
//      }
//      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "C3D6") )
//      {
//        shapetype[count] = DB_ZONETYPE_HEX;
//        writeArbitraryPolygon = true;
//      }
//      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "C3D4") )
//      {
//        shapetype[count] = DB_ZONETYPE_TET;
//      }
//      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "CPE4") ||
// !elementRegion.m_elementGeometryID.compare(0, 3, "S4R") )
//      {
//        shapetype[count] = DB_ZONETYPE_QUAD;
//      }
//      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "STRI") ||
// !elementRegion.m_elementGeometryID.compare(0, 4, "TRSH") ||
// !elementRegion.m_elementGeometryID.compare(0, 4, "CPE3"))
//      {
//        shapetype[count] = DB_ZONETYPE_TRIANGLE;
//      }
//      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "CPE2") )
//      {
//        shapetype[count] = DB_ZONETYPE_TRIANGLE;
//      }
//      else
//      {
//        GEOS_ERROR("PhysicalDomainT::WriteFiniteElementMesh: Do not recognize
// geometry type " + elementRegion.m_elementGeometryID + " \n");
//      }

        shapesize[count] = integer_conversion<int>(elemsToNodes.size(1));
        count++;
      }
      }
    }

    siloFile.WriteMeshObject(meshName, numNodes, coords,
                             nullptr, integer_conversion<int>(numElementRegions),
                             shapecnt.data(), meshConnectivity.data(), nullptr /*globalElementNumbers.data()*/,
                             isGhostElement.data(), shapetype.data(), shapesize.data(), cycleNum, problemTime);


    // write node fields in silo mesh, and all restart data as unassociated
    // variables.



    siloFile.WriteManagedGroupSilo( nodeManager, "NodalFields", meshName, DB_NODECENT, cycleNum, problemTime, isRestart, localIndex_array());



//    m_feElementManager->WriteSilo( siloFile, meshName, cycleNum, problemTime,
// isRestart );


  }//end FE write
}

void DomainPartition::ReadFiniteElementMesh( const SiloFile& siloFile,
                                             const int cycleNum,
                                             const realT problemTime,
                                             const bool isRestart )
{


//  int err = m_feNodeManager->ReadSilo( siloFile, "NodalFields", "volume_mesh",
//                                      DB_NODECENT, cycleNum, problemTime,
// isRestart );
////  err = m_feNodeManager->ReadSilo( siloFile, "NodalFieldsB", "face_mesh",
////                                      DB_NODECENT, cycleNum, problemTime,
// isRestart );
//  if(err)
//    return;
//
//  m_feElementManager->ReadSilo( siloFile, "volume_mesh",
//                               cycleNum, problemTime, isRestart );
//
//  m_feNodeManager->ConstructNodeToElementMap( m_feElementManager );

}





} /* namespace geosx */
