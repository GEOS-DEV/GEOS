/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file DofManager.cpp
 */

#include "DofManager.hpp"

namespace geosx
{

// .... COMM TOOLS :: ALL GATHER
//      TODO: move to CommunicationTools

namespace CommTools
{
  void allGather(localIndex const myValue, localIndex_array & allValues)
  {
    #ifdef GEOSX_USE_MPI
      int mpiRank; MPI_Comm_rank(MPI_COMM_GEOSX,&mpiRank);
      int mpiSize; MPI_Comm_size(MPI_COMM_GEOSX,&mpiSize);

      allValues.resize(mpiSize);
      array1d<int> tmpArray(mpiSize);

      int tmpValue = integer_conversion<int>(myValue);

      MPI_Allgather(&tmpValue,1,MPI_INT,tmpArray.data(),1,MPI_INT,MPI_COMM_GEOSX);

      for(localIndex i=0;i<tmpArray.size();++i)
        allValues[i] = integer_conversion<localIndex>(tmpArray[i]);
    #else
      allValues.resize(1);
      allValues[0] = myValue;
    #endif
  }
}


// .... DOF MANAGER :: CONSTRUCTOR

DofManager::DofManager()
{
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );

  // we pre-allocate an oversized array to store connectivity type
  // instead of resizing it dynamically as fields are added.

  m_connectivity.resize(MAX_NUM_FIELDS,MAX_NUM_FIELDS);
  
  for(localIndex i=0; i<MAX_NUM_FIELDS; ++i)
  for(localIndex j=0; j<MAX_NUM_FIELDS; ++j)
    m_connectivity[i][j] = Connectivity::None;
}


// .... DOF MANAGER :: SET MESH

void DofManager::setMesh(DomainPartition * const domain, 
                         localIndex const meshLevelIndex, 
                         localIndex const meshBodyIndex)
{
  GEOS_ERROR_IF(m_meshLevel != nullptr,"A mesh is already assigned to this DofManager.");
  m_domain = domain;
  m_meshLevel = m_domain->getMeshBodies()->GetGroup<MeshBody>(meshBodyIndex)->getMeshLevel(meshLevelIndex);;
}


// .... DOF MANAGER :: FIELD INDEX

localIndex DofManager::fieldIndex(string const & key)
{
  for(localIndex i=0; i<m_fields.size(); ++i)
    if(m_fields[i].name == key)
      return i;
  GEOS_ERROR("Field's string key not found in list of active fields.");
  return -1;
}


// .... DOF MANAGER :: KEY IN USE

bool DofManager::keyInUse(string const & key)
{
  bool check = false;
  for(localIndex i=0; i<m_fields.size(); ++i)
     check = (m_fields[i].name == key);
  return check;
}



// .... DOF MANAGER :: ADD FIELD

void DofManager::addField(string const & field,
                          Location const location, 
                          Connectivity const connectivity, 
                          localIndex const components,
                          string_array const & regions)
{
  // check if the field name is already being used

  GEOS_ERROR_IF(keyInUse(field),"Requested field name matches an existing field in the DofManager.");

  // save field description to list of active fields

  FieldDescription description;

  description.name = field;
  description.location = location;
  description.numComponents = components;
  description.regions = regions;                        // TODO: if regions is empty, add list of all regions
  description.key = field + "_dof_indices";
  description.docstring = field + " dof indices";

  if(components > 1)
    description.docstring += " (with " + std::to_string(components) + "-component blocks)";

  m_fields.push_back(description);

  GEOS_ERROR_IF(m_fields.size() > MAX_NUM_FIELDS, "Limit on DofManager's MAX_NUM_FIELDS exceeded.");

  // based on location, allocate an index array for this field

  localIndex last = m_fields.size()-1;

  switch(location)
  {
    case Location::Elem :
      createElemIndexArray(m_fields[last]);
      break;
    case Location::Face :
      createFaceIndexArray(m_fields[last]);
      break;
    case Location::Node :
      createNodeIndexArray(m_fields[last]);
      break;
    default:
      GEOS_ERROR("DoF support location is not yet supported");
  }
  
  // determine field's global offset

  if(last > 0)
    m_fields[last].fieldOffset = m_fields[last-1].fieldOffset + m_fields[last-1].numGlobalRows;
  else
    m_fields[last].fieldOffset = 0;

  // save field's connectivity type (self-to-self) 

  m_connectivity[last][last] = connectivity;

  // log some basic info
  
  GEOS_LOG_RANK_0("DofManager :: Added field .... " << m_fields[last].docstring);
  GEOS_LOG_RANK_0("DofManager :: Global dofs .... " << m_fields[last].numGlobalRows); 
  GEOS_LOG_RANK_0("DofManager :: Field offset ... " << m_fields[last].fieldOffset); 
}


// .... DOF MANAGER :: CREATE NODE INDEX ARRAY

void DofManager::createNodeIndexArray(FieldDescription & field)
{
  // step 0. register an index array with default = -1

  NodeManager * const nodeManager = m_meshLevel->getNodeManager();

  nodeManager->RegisterViewWrapper<globalIndex_array>( field.key )->
    setApplyDefaultValue(-1)->
    setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
    setDescription(field.docstring);
      
  globalIndex_array & indexArray = nodeManager->getReference<globalIndex_array>(field.key);

  // step 1. loop over all active regions 
  //         determine number of local rows
  //         and sequentially number objects

  ElementRegionManager * const elemManager = m_meshLevel->getElemManager();

  field.numLocalRows = 0;
  for(localIndex er=0; er<field.regions.size(); ++er)
  {
    ElementRegion * const elemRegion = elemManager->GetRegion(field.regions[er]);
    GEOS_ERROR_IF(elemRegion == nullptr,"Specified element region not found: " + field.regions[er]);

    for(localIndex esr=0; esr<elemRegion->numSubRegions(); esr++)
    {
      CellBlockSubRegion * subRegion = elemRegion->GetSubRegion(esr);

      integer_array const & ghostRank = subRegion->m_ghostRank;
      localIndex_array2d const & nodeMap = subRegion->getWrapper<FixedOneToManyRelation>(subRegion->viewKeys().nodeList)->reference();

      for(localIndex e=0; e<nodeMap.size(0); ++e)
      if(ghostRank[e] < 0)
      {
        for(localIndex n=0; n<nodeMap.size(1); ++n)
        {
          localIndex node = nodeMap[e][n];
          if(indexArray[node] == -1)
          {
            indexArray[node] = field.numLocalRows;
            field.numLocalRows++;
          }
        }
      }
    }
  }

  // step 2. gather row counts across ranks

  localIndex_array gather;

  CommTools::allGather(field.numLocalRows,gather); 

  field.numGlobalRows = 0;
  for(localIndex p=0; p<mpiSize; ++p)
    field.numGlobalRows += gather[p];

  field.firstLocalRow = 0;
  for(localIndex p=0; p<mpiRank; ++p)
    field.firstLocalRow += gather[p];

  // step 3. adjust local values to reflect processor offset

  for(localIndex node=0; node < indexArray.size(); ++node)
  if(indexArray[node] != -1)
    indexArray[node] += field.firstLocalRow;

  // step 4. synchronize across ranks

  std::map<string,string_array> fieldNames;
  fieldNames["node"].push_back(field.key);

  CommunicationTools::SynchronizeFields(fieldNames,m_meshLevel,
    m_domain->getReference< array1d<NeighborCommunicator> >( m_domain->viewKeys.neighbors ) );

  // step 5. scale row counts by number of vector components

  if(field.numComponents > 1)
  {
    field.numGlobalRows *= field.numComponents;
    field.numLocalRows  *= field.numComponents;
    field.firstLocalRow *= field.numComponents;
  }
}


// .... DOF MANAGER :: CREATE ELEM INDEX ARRAY
//      TODO: revise to look more like node version.
//            may even be able to condense to one function.

void DofManager::createElemIndexArray(FieldDescription & field)
{
  // step 1. loop over all active regions 
  //         determine number of local rows

  ElementRegionManager * const elemManager = m_meshLevel->getElemManager();

  field.numLocalRows = 0;
  for(localIndex er=0; er<field.regions.size(); ++er)
  {
    ElementRegion * const elemRegion = elemManager->GetRegion( field.regions[er]);
    GEOS_ERROR_IF(elemRegion == nullptr,"Specified element region not found: " + field.regions[er]);

    elemRegion->forCellBlocks([&]( CellBlockSubRegion * const subRegion )
    {
      localIndex numGhost = subRegion->GetNumberOfGhosts();
      field.numLocalRows  += subRegion->size() - numGhost;
    });
  }

  // step 2. gather row counts across ranks

  localIndex_array gather;

  CommTools::allGather(field.numLocalRows,gather); 

  field.numGlobalRows = 0;
  for(localIndex p=0; p<mpiSize; ++p)
    field.numGlobalRows += gather[p];

  field.firstLocalRow = 0;
  for(localIndex p=0; p<mpiRank; ++p)
    field.firstLocalRow += gather[p];

  // step 3. loop again (sequential policy)
  //         allocate the index array
  //         set unique global indices

  globalIndex const isUnset = -1;
  globalIndex count = 0;

  for(localIndex er=0; er<field.regions.size(); ++er)
  {
    ElementRegion * const elemRegion = elemManager->GetRegion( field.regions[er]);

    for(localIndex esr=0; esr < elemRegion->numSubRegions(); esr++)
    {
      CellBlockSubRegion * subRegion = elemRegion->GetSubRegion(esr);

      subRegion->RegisterViewWrapper<globalIndex_array>( field.key )->
        setApplyDefaultValue(isUnset)->
        setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->   // TODO: level 1 or 2?
        setDescription(field.docstring);
   
      globalIndex_array & indexArray = subRegion->getReference<globalIndex_array>(field.key);
      integer_array const & ghostRank = subRegion->m_ghostRank;

      GEOS_ERROR_IF(indexArray.size() != ghostRank.size(),"Mismatch in ghost rank and index array sizes.");

      for(localIndex elem=0; elem<ghostRank.size(); ++elem)
      if(ghostRank[elem] < 0)
      {
        indexArray[elem] = field.firstLocalRow+count; 
        count++;
      }
    };
  }

  GEOS_ERROR_IF(count != field.numLocalRows,"Mismatch during assignment of local row indices");

  // step 4. synchronize across ranks

  std::map<string,string_array> fieldNames;
  fieldNames["elems"].push_back(field.key);

  CommunicationTools::SynchronizeFields(fieldNames,m_meshLevel,
    m_domain->getReference< array1d<NeighborCommunicator> >( m_domain->viewKeys.neighbors ) );

  // step 5. scale row counts by number of vector components

  field.numGlobalRows *= field.numComponents;
  field.numLocalRows  *= field.numComponents;
  field.firstLocalRow *= field.numComponents;
}




// .... DOF MANAGER :: CREATE FACE INDEX ARRAY

void DofManager::createFaceIndexArray(FieldDescription & field)
{

  FaceManager * const faceManager = m_meshLevel->getFaceManager();

  faceManager->RegisterViewWrapper<globalIndex_array>( field.key )->
    setApplyDefaultValue(-1)->
    setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
    setDescription(field.docstring);

  // .... to be written ....
}


void DofManager::addCoupling(string const & rowField,
                             string const & colField,
                             Connectivity const connectivity,
                             bool const symmetric)
{
  // ... to be written ...
} 


}

