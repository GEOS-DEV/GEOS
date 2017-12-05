//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
