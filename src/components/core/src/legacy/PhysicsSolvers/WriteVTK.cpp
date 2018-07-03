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

/*
 * WriteVTK.cpp
 *
 *  Created on: May 27, 2015
 *      Author: stuartwalsh
 */



#include "WriteVTK.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"

#include "SolverFactory.h"

#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"


#include "ObjectManagers/FaceManagerT.h"
#include "ObjectManagers/EdgeManagerT.h"
#include "DataStructures/VectorFields/NodeManagerT.h"


#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"

std::string GetVTKFieldTypePrimitive(FieldType type){

  std::string rvalue;
  switch(type)
  {
  case FieldInfo::integerField:
  case FieldInfo::localIndexField:
  case FieldInfo::globalIndexField:
  case FieldInfo::integerParameter:
    rvalue = "long";  break;
  case FieldInfo::realField:
  case FieldInfo::R1TensorField:
  case FieldInfo::R2TensorField:
  case FieldInfo::R2SymTensorField:
  case FieldInfo::realParameter:
    rvalue = "double"; break;
  case FieldInfo::numFieldTypes:
  default:
    throw GPException("AddKeylessDataField: Unrecognized field type "+ type);
  }
  return rvalue;

}

WriteVTK::WriteVTK(  const std::string& name, ProblemManagerT* const pm ):
  WriteFieldToFile(name,pm)
{

  m_appendToFile = false;
  m_isFirstTime = true;
}

void WriteVTK::ReadXML(TICPP::HierarchicalDataNode* const hdn){
  WriteFieldToFile::ReadXML(hdn);
  m_appendToFile = false;
  m_isFirstTime = true;
}


/**
 *
 *
 **/
double WriteVTK::TimeStep( const realT& time,
                           const realT& dt,
                           const int cycleNumber,
                           PhysicalDomainT& domain,
                           const array<string>& namesOfSolverRegions,
                           SpatialPartition& partition,
                           FractunatorBase* const fractunator)
{

  realT dt_return = dt;

  m_stabledt.m_maxdt = std::numeric_limits<double>::max()*0.9;

  if(time >= m_outputTimes.front()*(1.0-1e-15) )
  {
    if(partition.m_rank == 0)
      std::cout << "Writing data field(s) to vtk file." << std::endl;

    std::string filename = m_filePrefix;
    if(partition.m_size > 1)
      filename += "_" + toString(partition.m_rank);
    filename += "_" + toString(time) + ".vtk";


    // output header
    std::ofstream fStream(filename.c_str(),std::ios::out);
    fStream << "# vtk DataFile Version 2.0\n";
    fStream << m_headerString;
    //fStream << "# Time: " << time << "\n"; // slightly redundant as time is in
    // filename
    fStream << "ASCII\n"
            << "DATASET UNSTRUCTURED_GRID\n";


    NodeManager& nodeManager = domain.m_feNodeManager;
    FaceManagerT& faceManager = domain.m_feFaceManager;
    // output node locations

    localIndex numNodes = nodeManager.DataLengths();
    fStream << "POINTS " << numNodes << " double\n";
    for( localIndex kn=0 ; kn<numNodes ; ++kn )
    {
      R1Tensor nodePosition = (*nodeManager.m_refposition)[kn];
      nodePosition += (*nodeManager.m_displacement)[kn];
      fStream << nodePosition[0] << " " << nodePosition[1] << " " << nodePosition[2] << " "  << "\n";
    }
    fStream << "\n";

    // output faces

    localIndex numFaces = faceManager.DataLengths();
    localIndex cumEntryCount = 0;
    for( localIndex kf=0 ; kf<numFaces ; ++kf )
    {
      const lArray1d& faceNodeMap = faceManager.m_toNodesRelation[kf];
      cumEntryCount += faceNodeMap.size()+1;
    }

    fStream << "CELLS " << numFaces << " " << cumEntryCount << "\n";

    // connectivity
    for( localIndex kf=0 ; kf<numFaces ; ++kf )
    {

      const lArray1d& faceNodeMap = faceManager.m_toNodesRelation[kf];
      fStream << faceNodeMap.size();
      for( localIndex a=0 ; a<faceNodeMap.size() ; ++a )
      {
        //const localIndex nd = faceNodeMap[a];
        localIndex nd = faceNodeMap[a];

        fStream << " " << nd;

      }
      fStream << "\n";

    }

    fStream << "CELL_TYPES " << numFaces << "\n";

    for( localIndex kf=0 ; kf<numFaces ; ++kf )
    {

      const lArray1d& faceNodeMap = faceManager.m_toNodesRelation[kf];
      localIndex numFaceNodes = faceNodeMap.size();
      int ctype = 7;   // polygon
      if(numFaceNodes == 3)
      {
        ctype = 5;   // tri
      }
      else if(numFaceNodes == 4)
      {
        ctype = 9;   // quad
      }
      fStream << ctype <<"\n";

    }

    ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectType,m_regionName);

    if(m_objectType == PhysicalDomainT::FiniteElementFaceManager)
    {
      fStream << "CELL_DATA " << numFaces << "\n";
    }
    else if(m_objectType == PhysicalDomainT::FiniteElementNodeManager)
    {
      fStream << "POINT_DATA " << numNodes << "\n";
    }

    localIndex numFaceFields = m_fieldNames.size();
    fStream << "FIELD fieldData " << numFaceFields << "\n";

    for(localIndex i =0 ; i < numFaceFields ; ++i)
    {
      fStream << m_fieldNames[i] << " " << FieldInfo::FieldSize( m_fieldTypes[i] ) << " " << numFaces << " " << GetVTKFieldTypePrimitive(m_fieldTypes[i]) <<
        "\n";

      fStream.close();

      objectManager.WriteAsciiFieldData( m_fieldTypes[i], m_fieldNames[i], filename, true);

      fStream.open(filename.c_str(), std::ios::out | std::ios::app);
    }



    fStream.close();

    /*

       ObjectDataStructureBaseT& objectManager =
          domain.GetObjectDataStructure(m_objectType,m_regionName);

       if(m_objectType == PhysicalDomainT::FiniteElementElementRegion){
       ElementRegionT *elementRegion = dynamic_cast<ElementRegionT
     *>(&objectManager);
       elementRegion->UpdateElementFieldsWithGaussPointData();
       }

       if( m_setNames.empty() )
       {

       objectManager.WriteAsciiFieldData(m_fieldTypes,
           m_fieldNames,filename,true);

       }else{
       for(unsigned i = 0; i < m_setNames.size(); ++i)
       {
        lSet& set = objectManager.GetSet(m_setNames[i]);
        objectManager.WriteAsciiFieldData(m_fieldTypes,
            m_fieldNames,filename,set,true);
       }
       }
     */


    if(m_outputTimes.back() < std::numeric_limits<double>::max() )
    {
      m_outputTimes.push_back(m_outputTimes.back() + m_dt);
    }

    m_outputTimes.pop_front();
  }

  return dt_return;
}

REGISTER_SOLVER(WriteVTK)
