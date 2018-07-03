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
 * @file CartesianGridManagerT.h
 * @author walsh24
 * @date Apr 6, 2011
 */

#ifndef CartesianGridManager_H_
#define CartesianGridManager_H_

//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "../../dataRepository/Group.hpp"
#include "../IO/ticpp/HierarchicalDataNode.h.old"

class FaceManagerT;
class NodeManager;
class EdgeManagerT;

class CartesianGridManagerT : public ObjectDataStructureBaseT
{
public:
  CartesianGridManagerT();
  ~CartesianGridManagerT();

  void Initialize(){}

  void SetDomainBoundaryObjects(const ObjectDataStructureBaseT* const referenceObject ){
    // get the "isDomainBoundary" field from for *this, and set it to zero
    array<integer>& isDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();
    isDomainBoundary = 0;
    /*
     * Loop over array boundaries and set to true if borders a partition
     * */

    int kL = m_numGhostCells[2][0]-1;
    int kU = m_nZ-m_numGhostCells[2][1];
    for(localIndex j = 0 ; j < m_nY ; ++j)
    {
      for(localIndex i = 0 ; i < m_nX ; ++i)
      {
        if(kL >= 0)
        {
          localIndex indxL = lxyz2ind(i,j,kL);
          isDomainBoundary[indxL] = 1;
        }
        if(kU < int(m_nZ))
        {
          localIndex indxU = lxyz2ind(i,j,kU);
          isDomainBoundary[indxU] = 1;
        }
      }
    }

    int jL = m_numGhostCells[1][0]-1;
    int jU = m_nY-m_numGhostCells[1][1];
    for(localIndex k = 0 ; k < m_nZ ; ++k)
    {
      for(localIndex i = 0 ; i < m_nX ; ++i)
      {
        if(jL >= 0)
        {
          localIndex indxL = lxyz2ind(i,jL,k);
          isDomainBoundary[indxL] = 1;
        }
        if(jU < int(m_nY))
        {
          localIndex indxU = lxyz2ind(i,jU,k);
          isDomainBoundary[indxU] = 1;
        }
      }
    }

    int iL = m_numGhostCells[0][0]-1;
    int iU = m_nX-m_numGhostCells[0][1];
    for(localIndex k = 0 ; k < m_nZ ; ++k)
    {
      for(localIndex j = 0 ; j < m_nY ; ++j)
      {
        if(iL >= 0)
        {
          localIndex indxL = lxyz2ind(iL,j,k);
          isDomainBoundary[indxL] = 1;
        }
        if(iU < int(m_nX))
        {
          localIndex indxU = lxyz2ind(iU,j,k);
          isDomainBoundary[indxU] = 1;
        }
      }
    }
  };
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject  = NULL){

    geosx::dataRepository::ViewWrapper<integer_array>::rtype m_isExternal = getData<integer_array>("IsExternal");

    // get the "isExternal" field from for *this, and set it to zero
    m_isExternal = 0;

    // Loop over array boundaries and set to true if external

    for(localIndex j = 0 ; j < m_nY ; ++j)
    {
      for(localIndex i = 0 ; i < m_nX ; ++i)
      {
        if(m_numGhostCells[2][0] == 0)
        {
          localIndex indxL = lxyz2ind(i,j,0);
          m_isExternal[indxL] = 1;
        }
        if(m_numGhostCells[2][1] == 0)
        {
          localIndex indxU = lxyz2ind(i,j,m_nZ-1);
          m_isExternal[indxU] = 1;
        }
      }
    }

    for(localIndex k = 0 ; k < m_nZ ; ++k)
    {
      for(localIndex i = 0 ; i < m_nX ; ++i)
      {
        if(m_numGhostCells[1][0] == 0)
        {
          localIndex indxL = lxyz2ind(i,0,k);
          m_isExternal[indxL] = 1;
        }
        if(m_numGhostCells[1][1] == 0)
        {
          localIndex indxU = lxyz2ind(i,m_nY-1,k);
          m_isExternal[indxU] = 1;
        }
      }
    }

    for(localIndex k = 0 ; k < m_nZ ; ++k)
    {
      for(localIndex j = 0 ; j < m_nY ; ++j)
      {
        if(m_numGhostCells[0][0] == 0)
        {
          localIndex indxL = lxyz2ind(0,j,k);
          m_isExternal[indxL] = 1;
        }
        if(m_numGhostCells[0][1] == 0)
        {
          localIndex indxU = lxyz2ind(m_nX-1,j,k);
          m_isExternal[indxU] = 1;
        }
      }
    }

  };
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                         array<gArray1d>& objectToCompositionObject  )
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
  T GetValue(array<T>& values, int x, int y, int z){
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
  array< R1Tensor >* const m_refposition;


private:

  void SetReferencePosition();

  int m_numGhostCells[3][2]; // number of ghost cells along each axis
  R1Tensor   m_global_origin;  // global origin of grid

};

#endif /* EDGEMANAGERT_H_ */
