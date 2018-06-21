// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * DiscreteElementManagerBaseT.cpp
 *
 *  Created on: July 12, 2011
 *      Author: Scott Johnson
 */
#include "DiscreteElementManagerBaseT.h"

/**
 * @brief Constructor to set internal pointers to the external node and face
 * managers
 * @author Scott Johnson
 */
DiscreteElementManagerBaseT::DiscreteElementManagerBaseT():
  ObjectDataStructureBaseT(ObjectDataStructureBaseT::DiscreteElementManager),
  writeVTK(true),
#if USECPP11==1
  m_mat()
#else
  m_mat(NULL)
#endif
{
#if USECPP11!=1
  if(m_mat)
    delete m_mat;
#endif
  m_mat = MaterialFactory::NewMaterial("LinearElasticDEMMaterial");
  m_mat->resize(0,1);
}

DiscreteElementManagerBaseT::DiscreteElementManagerBaseT( const ObjectType objectType ):
  ObjectDataStructureBaseT(objectType),
  writeVTK(true),
  #if USECPP11==1
  m_mat()
#else
  m_mat(NULL)
#endif

{
  m_mat = MaterialFactory::NewMaterial("LinearElasticDEMMaterial");
  m_mat->resize(0,1);
}

DiscreteElementManagerBaseT::~DiscreteElementManagerBaseT()
{
#if USECPP11!=1
  if(m_mat)
    delete m_mat;
#endif
}

globalIndex DiscreteElementManagerBaseT::insert(const localIndex i, const bool assignGlobals )
{
  globalIndex gi = ObjectDataStructureBaseT::insert(i, assignGlobals);
  m_mat->insert(i);
  return gi;
}
void DiscreteElementManagerBaseT::erase( const localIndex i )
{
  ObjectDataStructureBaseT::erase(i);
  m_mat->erase(i);
}
globalIndex DiscreteElementManagerBaseT::resize( const localIndex size, const bool assignGlobals )
{
  globalIndex gi = ObjectDataStructureBaseT::resize( size, assignGlobals );
  m_mat->resize(size);
  return gi;
}


/**
 * @brief Adds fields associated with all derived classes
 * @author Scott Johnson
 */
void DiscreteElementManagerBaseT::AddBaseFields()
{
  this->AddKeyedDataField<FieldInfo::incrementalDisplacement>();//-> specified
                                                                // in
                                                                // NodeManagerT
                                                                // constructor
  this->AddKeyedDataField<FieldInfo::displacement>();//-> specified in
                                                     // NodeManagerT constructor
  this->AddKeyedDataField<FieldInfo::referencePosition>();
  this->AddKeyedDataField<FieldInfo::currentPosition>();
  this->AddKeyedDataField<FieldInfo::rotationAxis>();
  this->AddKeyedDataField<FieldInfo::rotationMagnitude>();
  this->AddKeyedDataField<FieldInfo::rotationalInertia>();
  this->AddKeyedDataField<FieldInfo::velocity>();
  this->AddKeyedDataField<FieldInfo::acceleration>();
  this->AddKeyedDataField<FieldInfo::force>();
  this->AddKeyedDataField<FieldInfo::mass>();
  this->AddKeylessDataField<realT>("kinetic_energy", true, true);
}

/**
 * @author Scott Johnson
 * @param[in,out] buffer the buffer to unpack elements from
 * @param[in,out] elementReceiveLocalIndices the local indices of the elements
 * @return
 */
unsigned int DiscreteElementManagerBaseT::Unpack( const char*& buffer, lArray1d& elementReceiveLocalIndices )
{
  unsigned int sizeOfUnpacked = 0;

  int numReceivedElements;
  sizeOfUnpacked += bufvector::Unpack( buffer, numReceivedElements );
  const localIndex oldSize = this->m_DataLengths;
  elementReceiveLocalIndices.resize( numReceivedElements );

  // local variable to store global indices of new elements so that they can be
  // used to fill the localToGlobalMap
  gArray1d newGlobalIndices;

  int numNewElements = 0;
  for( int a=0 ; a<numReceivedElements ; ++a )
  {
    // read the global index from the buffer
    globalIndex gElementIndex;
    sizeOfUnpacked += bufvector::Unpack( buffer, gElementIndex );

    // check to see if the object already exists by checking for the global
    // index in m_globalToLocalMap. If it doesn't,
    // then add the element to the maps, and increment numNewElements.
    std::map<globalIndex,localIndex>::iterator iterG2L = m_globalToLocalMap.find(gElementIndex);
    if( iterG2L == m_globalToLocalMap.end() )
    {
      // the global index is not contained in this object, so we should add the
      // element.

      // add to the element information to the globalToLocalMap
      m_globalToLocalMap[gElementIndex] = this->DataLengths() + numNewElements;

      // add the local index of the element do the elementReceiveLocalIndices
      // array
      elementReceiveLocalIndices(a) = oldSize + numNewElements;

      // add entry to newGlobalIndices
      newGlobalIndices.push_back( gElementIndex );

      // increment the number of new elements
      ++numNewElements;
    }
    else
    {
      // the global index is contained in this object, so we should point
      // everything at the existing element

      // get the local index of the element
      localIndex b = iterG2L->second;
      elementReceiveLocalIndices(a) = b;

    }
  }

  // figure out new sizes, and resize object
  const localIndex newSize = oldSize + numNewElements;
  this->resize( newSize );

  // add entires to the localToGlobalMap
  for( int a=0 ; a<numNewElements ; ++a )
  {
    localIndex b = oldSize + a;
    m_localToGlobalMap[b] = newGlobalIndices(a);
  }

  // Unpack fields
  sizeOfUnpacked += UnpackAllFieldsFromBuffer( buffer, elementReceiveLocalIndices );

  sizeOfUnpacked += UnpackSets( buffer );

  return sizeOfUnpacked;
}

/**
 * @author Scott Johnson
 * @param sendElements local indices of elements to pack and send
 * @param buffer the buffer to pack the elements into
 * @return size of characters packed into the buffer.
 *
 * This function packs complete elements into a buffer. this should include all
 * information needed to reconstruct the element
 * on a remote process domain. This does not include maps to other objects, as
 * those are locally indexed relations and
 * must be constructed on the receiving domain.
 */
unsigned int DiscreteElementManagerBaseT::Pack( const lArray1d& sendElements, bufvector& buffer ) const
{
  unsigned int sizeOfPacked = 0;

  // the number of elements that will be packed is the size of the sendElements
  // array
  int numElements = sendElements.size();

  // pack the number of elements
  sizeOfPacked += buffer.Pack( numElements );

  // iterate over the entries in sendElements
  for(lArray1d::const_iterator elementIndex=sendElements.begin() ; elementIndex!=sendElements.end() ; ++elementIndex )
  {
    // pack the global index
    globalIndex gNodeIndex = this->m_localToGlobalMap(*elementIndex);
    sizeOfPacked += buffer.Pack(gNodeIndex);
  }

  // pack all fields
  sizeOfPacked += ObjectDataStructureBaseT::PackAllFieldsIntoBuffer( buffer, sendElements );

  PackSets( sendElements, buffer);

  return sizeOfPacked;
}

void DiscreteElementManagerBaseT::ReadXML(TICPP::HierarchicalDataNode* hdn)
{
  this->writeVTK = hdn->GetAttributeOrDefault<int> ("writeVTK", 1) > 0;
}

void DiscreteElementManagerBaseT::WriteVTK(const int cycleNum, ContactManagerBaseT& contacts)
{
  //open file
  char buffer[30];
#if GPAC_MPI
  {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > 1)
    {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0)
      {
        sprintf(buffer, "demall_%06d.visit", cycleNum);
        std::ofstream outAll(buffer);
        outAll << "!NBLOCKS " << size << std::endl;
        for (int i = 0 ; i < size ; i++)
        {
          sprintf(buffer, "dem_%06d_%06d.vtk", i, cycleNum);
          outAll << buffer << std::endl;
        }
        outAll.close();
      }
      sprintf(buffer, "dem_%06d_%06d.vtk", rank, cycleNum);
    }
    else
    {
      sprintf(buffer, "dem_%06d.vtk", cycleNum);
    }
  }
#else
  sprintf(buffer, "dem_%06d.vtk",cycleNum);
#endif
  std::ofstream out(buffer);

  //write header
  out << "# vtk DataFile Version 3.0" << std::endl << "vtk output" << std::endl << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl <<std::endl;

  //write points
  out << "POINTS " << m_DataLengths << " double" << std::endl;
  {
    const array<R1Tensor>& deCurrentPosition         = GetFieldData<FieldInfo::currentPosition> ();
    for(array<R1Tensor>::const_iterator it = deCurrentPosition.begin() ; it != deCurrentPosition.end() ; ++it)
    {
      out << (*it)(0) << " " << (*it)(1) << " " << (*it)(2) << std::endl;
    }
  }

//  //write vertices
//  out << std::endl << "VERTICES " << m_DataLengths << " " << (2*m_DataLengths)
// << std::endl;
//  for(localIndex i = 0; i < this->m_DataLengths; i++)
//  {
//    out << "1 " << i << std::endl;
//  }

  //write lines
  out << std::endl << "LINES " << contacts.DataLengths() << " " << (3*contacts.DataLengths()) << std::endl;
  {
    const lArray1d& v0 = contacts.GetFieldData<localIndex>( "face1");
    const lArray1d& v1 = contacts.GetFieldData<localIndex>( "face2");
    for(localIndex i = 0 ; i < contacts.DataLengths() ; i++)
    {
      out << "2 " << v0(i) << " " << v1(i) << std::endl;
    }
  }

  //write point data
  out << std::endl << "POINT_DATA " << m_DataLengths << std::endl;
  WriteVTKPointData(out);

  //write line data
  out << std::endl << "CELL_DATA " << contacts.DataLengths() << std::endl;
  contacts.WriteVTKCellData(out);

  //close file
  out.close();
}

void DiscreteElementManagerBaseT::WriteVTKPointData(std::ofstream& out)
{
  for( std::map<std::string, array<real64> >::const_iterator itn=m_realData.begin() ; itn!=m_realData.end() ; ++itn )
  {
    out << "SCALARS pt" << itn->first << " double" << std::endl << "LOOKUP_TABLE default" << std::endl;
    const array<real64>& scalar = itn->second;
    for(array<real64>::const_iterator it = scalar.begin() ; it != scalar.end() ; ++it)
    {
      out << (*it) << " ";
    }
    out << std::endl;
  }

  for( std::map<std::string, array<R1Tensor> >::const_iterator itn=m_R1TensorData.begin() ; itn!=m_R1TensorData.end() ; ++itn )
  {
    out << "SCALARS pt" << itn->first << "_Magnitude double" << std::endl << "LOOKUP_TABLE default" << std::endl;
    const array<R1Tensor>& scalar = itn->second;
    for(array<R1Tensor>::const_iterator it = scalar.begin() ; it != scalar.end() ; ++it)
    {
      out << it->L2_Norm() << " ";
    }
    out << std::endl;
  }

  for( std::map<std::string, array<integer> >::const_iterator itn=m_IntegerData.begin() ; itn!=m_IntegerData.end() ; ++itn )
  {
    out << "SCALARS pt" << itn->first << " integer" << std::endl << "LOOKUP_TABLE default" << std::endl;
    const array<integer>& scalar = itn->second;
    for(array<integer>::const_iterator it = scalar.begin() ; it != scalar.end() ; ++it)
    {
      out << (*it) << " ";
    }
    out << std::endl;
  }
}


void DiscreteElementManagerBaseT::ForceToCouple(const R1Tensor& position,
                                                const R1Tensor& cforce,
                                                const R1Tensor& currentPosition,
                                                const R2Tensor& rotation,
                                                R1Tensor& force,
                                                R1Tensor& moment)
{
  //get the local frame force
  R1Tensor fl;
  fl.AijBj(rotation, cforce);

  R1Tensor tmp(position);
  tmp -= currentPosition;
  R1Tensor xl;
  xl.AijBj(rotation, tmp);

  //get the incremental moment
  tmp.Cross(xl, fl);

  moment += tmp;
  force += cforce;
}
