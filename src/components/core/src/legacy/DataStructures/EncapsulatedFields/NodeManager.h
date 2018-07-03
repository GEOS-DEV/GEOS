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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file NodeManager.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */


#ifndef NODEMANAGER_H_
#define NODEMANAGER_H_

#include "Common/Common.h"
#include "ObjectManager.h"


/**
 * @author Randolph Settgast
 *
 * The NodeManager class manages the node data as a collection of NodeT objects,
 * each containing
 * their own field data.
 */
class NodeManager : public ObjectManager
{
public:


  /// default constructor
  NodeManager();

  /// copy constructor
  NodeManager( const NodeManager& init );

  /// default destructor
  ~NodeManager();

  virtual const std::type_info& get_typeid() const override final
  {
    return typeid(*this);
  }

  /// reads nodal input from the preliminary geometry file format
  void ReadAsciiNodeInput( std::ifstream& geometryStream );

  /// copy fields from one node to another
  void CopyNode( const int destination, const int source );

  /// translate a node
  void TranslateNode( const int nodenum, const R1Tensor& offset );

  /// copy some field data from the global array to a local array
  template < FieldKey FIELDNAME, typename T>
  void CopyGlobalFieldToLocalField( const int* const nodelist,
                                    array< T >& localField ) const;

  /// copy a pair of field data from their global arrays to local arrays
  template < FieldKey FIELDNAME1, FieldKey FIELDNAME2, typename T >
  void CopyGlobalFieldsToLocalFields( const int* __restrict__ const nodelist,
                                      array< T >& localField1,
                                      array< T >& localField2 ) const;

  /// add the contents of a local field to a global field
  template < FieldKey FIELDNAME, typename T>
  void AddLocalFieldToGlobalField( const int* const nodelist,
                                   const array< T >& localField );

  /// write the node data to ensight6
  template < FieldKey FIELDKEY, typename T >
  void WriteFieldToEnsight( const int cycle,
                            const realT time,
                            const std::string& fileroot ) const;

  /// number of nodes
  const size_t& m_numNodes;

  /// vector of NodeT objects
//  std::vector< NodeT* > m_nodes;


protected:

private:
  NodeManager& operator=( const NodeManager&);

};



/**
 * @author R. Settgast
 * @tparam FIELDNAME name of the field to copy
 * @tparam T type of the field
 * @param nodelist the list of local node numbers to copy
 * @param localField array to copy field data to
 */
template < FieldKey FIELDNAME, typename T>
inline void NodeManager::CopyGlobalFieldToLocalField( const int* __restrict__ const nodelist,
                                                      array< T >& localField ) const
{

  const int N = localField.size();

  for( int a=0 ; a<N ; ++a )
  {
    localField[a] = m_objects[ nodelist[a] ]->GetField<FIELDNAME>();
  }

}

/**
 * @author R. Settgast
 * @tparam FIELDNAME1 name of the first field to copy
 * @tparam FIELDNAME2 name of the second field to copy
 * @tparam T type of the fields
 * @param nodelist the list of local node numbers to copy
 * @param localField1 array to copy first field data to
 * @param localField2 array to copy second field data to
 */
template < FieldKey FIELDNAME1, FieldKey FIELDNAME2, typename T>
inline void NodeManager::CopyGlobalFieldsToLocalFields( const int* __restrict__ const nodelist,
                                                        array< T >& localField1,
                                                        array< T >& localField2 ) const
{

  const int N = localField1.size();

  for( int a=0 ; a<N ; ++a )
  {
    const Object& node = *(m_objects[ nodelist[a] ]);
    localField1[a] = node.GetField<FIELDNAME1>();
    localField2[a] = node.GetField<FIELDNAME2>();


  }

}

/**
 * @author R. Settgast
 * @tparam FIELDNAME name of the field to add data to
 * @tparam T type of the field
 * @param nodelist the list of local node numbers to copy
 * @param localField array of data to add
 */
template < FieldKey FIELDNAME, typename T>
inline void NodeManager::AddLocalFieldToGlobalField( const int* const nodelist,
                                                     const array< T >& localField )
{
  const int N = localField.size();

  for( int a=0 ; a<N ; ++a )
  {
    m_objects[ nodelist[a] ]->GetField<FIELDNAME>() += localField[a];
  }

}



#endif /* NodeManager_H_ */
