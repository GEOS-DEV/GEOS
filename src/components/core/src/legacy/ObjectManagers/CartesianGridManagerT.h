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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file CartesianGridManagerT.h
 * @author walsh24
 * @date Apr 6, 2011
 */

#ifndef CartesianGridManager_H_
#define CartesianGridManager_H_

//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "../../dataRepository/Group.hpp"
#include "IO/ticpp/HierarchicalDataNode.h"

class FaceManagerT;
class NodeManagerT;
class EdgeManagerT;

class CartesianGridManagerT: public ObjectDataStructureBaseT
{
public:
  CartesianGridManagerT();
  ~CartesianGridManagerT();

  void Initialize(){}

  void SetDomainBoundaryObjects(const ObjectDataStructureBaseT* const referenceObject ){
    // get the "isDomainBoundary" field from for *this, and set it to zero
    iArray1d& isDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();
    isDomainBoundary = 0;
    /*
     * Loop over array boundaries and set to true if borders a partition
     * */

    int kL = m_numGhostCells[2][0]-1;
    int kU = m_nZ-m_numGhostCells[2][1];
    for(localIndex j = 0;j < m_nY; ++j){
    for(localIndex i = 0;i < m_nX; ++i){
      if(kL >= 0){
        localIndex indxL = lxyz2ind(i,j,kL);
        isDomainBoundary[indxL] = 1;
      }
      if(kU < int(m_nZ)){
        localIndex indxU = lxyz2ind(i,j,kU);
        isDomainBoundary[indxU] = 1;
      }
    }
    }

    int jL = m_numGhostCells[1][0]-1;
    int jU = m_nY-m_numGhostCells[1][1];
    for(localIndex k = 0;k < m_nZ; ++k){
    for(localIndex i = 0;i < m_nX; ++i){
      if(jL >= 0){
        localIndex indxL = lxyz2ind(i,jL,k);
        isDomainBoundary[indxL] = 1;
      }
      if(jU < int(m_nY)){
        localIndex indxU = lxyz2ind(i,jU,k);
        isDomainBoundary[indxU] = 1;
      }
    }
    }

    int iL = m_numGhostCells[0][0]-1;
    int iU = m_nX-m_numGhostCells[0][1];
    for(localIndex k = 0;k < m_nZ; ++k){
    for(localIndex j = 0;j < m_nY; ++j){
      if(iL >= 0){
        localIndex indxL = lxyz2ind(iL,j,k);
        isDomainBoundary[indxL] = 1;
      }
      if(iU < int(m_nX)){
        localIndex indxU = lxyz2ind(iU,j,k);
        isDomainBoundary[indxU] = 1;
      }
    }
    }
  };
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject  = NULL){

    geosx::dataRepository::ViewWrapper<int32_array>::rtype m_isExternal = getData<int32_array>("IsExternal");

    // get the "isExternal" field from for *this, and set it to zero
    m_isExternal = 0;

    // Loop over array boundaries and set to true if external

    for(localIndex j = 0;j < m_nY; ++j){
    for(localIndex i = 0;i < m_nX; ++i){
      if(m_numGhostCells[2][0] == 0){
        localIndex indxL = lxyz2ind(i,j,0);
        m_isExternal[indxL] = 1;
      }
      if(m_numGhostCells[2][1] == 0){
        localIndex indxU = lxyz2ind(i,j,m_nZ-1);
        m_isExternal[indxU] = 1;
      }
    }
    }

    for(localIndex k = 0;k < m_nZ; ++k){
    for(localIndex i = 0;i < m_nX; ++i){
      if(m_numGhostCells[1][0] == 0){
        localIndex indxL = lxyz2ind(i,0,k);
        m_isExternal[indxL] = 1;
      }
      if(m_numGhostCells[1][1] == 0){
        localIndex indxU = lxyz2ind(i,m_nY-1,k);
        m_isExternal[indxU] = 1;
      }
    }
    }

    for(localIndex k = 0;k < m_nZ; ++k){
    for(localIndex j = 0;j < m_nY; ++j){
      if(m_numGhostCells[0][0] == 0){
        localIndex indxL = lxyz2ind(0,j,k);
        m_isExternal[indxL] = 1;
      }
      if(m_numGhostCells[0][1] == 0){
        localIndex indxU = lxyz2ind(m_nX-1,j,k);
        m_isExternal[indxU] = 1;
      }
    }
    }

   };
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager ,
                                                         Array1dT<gArray1d>& objectToCompositionObject  )
  {
    (void)compositionObjectManager;
    (void)objectToCompositionObject;
    throw GPException("CartesianGridManagerT::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n");
   };

  void ReadXML(TICPP::HierarchicalDataNode* hdn);

  // ind2sub
  // convert local index to local (current processor) x,y,z coordinates
  void ind2lxyz(int ind, int x[3]){
    x[2] = ind/m_nXY; ind -= x[2]*m_nXY;
    x[1] = ind/m_nX; ind -= x[1]*m_nX;
    x[0] = ind;
  }

  // convert local x,y,z coordinates to local index
  int lxyz2ind(int x, int y, int z){
    return m_nXY*z + m_nX*y +x;
  }

  // get value at local x,y,z coordinate
  template<class T>
  T GetValue(Array1dT<T>& values, int x, int y, int z){
     return values[ lxyz2ind(x, y, z) ];
  }

  // data
  realT m_dx;                  // regular grid spacing

  localIndex m_global_dims[3]; // dimensions of grid    
  localIndex m_local_dims[3];  // dimensions of grid on this processor
  localIndex m_offset[3];      // xyz offset of origin grid on this processor

  const localIndex& m_nX; // references to local dimensions
  const localIndex& m_nY; 
  const localIndex& m_nZ; 
  localIndex m_nXY;
  const localIndex& m_nXYZ;


  /// Reference position
  Array1dT< R1Tensor >* const m_refposition;


private:

  void SetReferencePosition();

  int m_numGhostCells[3][2]; // number of ghost cells along each axis
  R1Tensor   m_global_origin;  // global origin of grid

};

#endif /* EDGEMANAGERT_H_ */
