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
 * @file CartesianGridManagerT.cpp
 * @author walsh24
 * @date April 6, 2012
 */

#include "CartesianGridManagerT.h"

#include "DataStructures/VectorFields/NodeManagerT.h"
#include "FaceManagerT.h"

CartesianGridManagerT::CartesianGridManagerT():
  ObjectDataStructureBaseT(ObjectDataStructureBaseT::CartesianGridManager),
  m_dx(1.0),
  m_nX(m_local_dims[0]),
  m_nY(m_local_dims[1]),
  m_nZ(m_local_dims[2]),
  m_nXY(0),
  m_nXYZ(this->m_DataLengths),
  m_refposition(NULL)
{
  m_global_dims[0] = 1;
  m_global_dims[1] = 1;
  m_global_dims[2] = 1;
  m_local_dims[0] = 1;
  m_local_dims[1] = 1;
  m_local_dims[2] = 1;
  m_offset[0] = 0;
  m_offset[1] = 0;
  m_offset[2] = 0;

  for(int d =0 ; d < 3 ; ++d )
  {
    for(int lu = 0 ; lu < 2 ; ++lu)
    {
      m_numGhostCells[d][lu] = 0;
    }
  }


  this->AddKeyedDataField<FieldInfo::referencePosition>();
  SetConstPointer<FieldInfo::referencePosition>( m_refposition );
}

CartesianGridManagerT::~CartesianGridManagerT()
{
  // TODO Auto-generated destructor stub
}


void CartesianGridManagerT::ReadXML(TICPP::HierarchicalDataNode* hdn){
  m_global_dims[0] = hdn->GetAttributeOrDefault<localIndex>("nx",1);
  m_global_dims[1] = hdn->GetAttributeOrDefault<localIndex>("ny",1);
  m_global_dims[2] = hdn->GetAttributeOrDefault<localIndex>("nz",1);
  m_dx = hdn->GetAttributeOrDefault<realT>("dx",1.0);

  {
    R1Tensor zero(0.0);
    m_global_origin = hdn->GetAttributeTensorOrDefault("origin",zero);
  }

  // fixme assumes 1 processor
  m_local_dims[0] = m_global_dims[0];
  m_local_dims[1] = m_global_dims[1];
  m_local_dims[2] = m_global_dims[2];
  m_offset[0] = 0;
  m_offset[1] = 0;
  m_offset[2] = 0;



  // fixme need to set number of ghost cells and check for periodic
  // bc's/external bcs
  for(int d =0 ; d < 3 ; ++d )
  {
    for(int lu = 0 ; lu < 2 ; ++lu)
    {
      m_numGhostCells[d][lu] = 0; // fixme
    }
  }

  m_nXY = m_local_dims[0]*m_local_dims[1];

  this->resize(m_local_dims[0]*m_local_dims[1]*m_local_dims[2]);

  SetReferencePosition();


  std::cout << "Cartesian Grid:"<< std::endl;
  std::cout << "    Spacing:    " << m_dx   << std::endl;
  std::cout << "    Dimensions: "<< std::endl;
  std::cout << "        x: " << m_global_dims[0] << std::endl;
  std::cout << "        y: " << m_global_dims[1] << std::endl;
  std::cout << "        z: " << m_global_dims[2] << std::endl;
}

void CartesianGridManagerT::SetReferencePosition(){

  localIndex indx = 0;

  for(localIndex k = 0 ; k < m_nZ ; ++k)
  {
    realT z = (k+ m_offset[2])*m_dx + m_global_origin[2];
    for(localIndex j = 0 ; j < m_nY ; ++j)
    {
      realT y = (j+ m_offset[1])*m_dx + m_global_origin[1];
      for(localIndex i = 0 ; i < m_nX ; ++i)
      {
        realT x = (i+ m_offset[0])*m_dx + m_global_origin[0];

        (*m_refposition)[indx] = R1Tensor(x,y,z);

        ++indx;
      }
    }
  }

}
