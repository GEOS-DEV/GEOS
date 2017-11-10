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
 * File: FileManagerT.cpp
 * Class provides file IO
 * created : RRS (10/11/2001)
 */

#include "FileManagerT.h"

#if GPAC_MPI
#include <mpi.h>
#endif
#include <stdlib.h>
#include <map>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include "metis.h"

#include "Common/GlobalIndexManager.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ElementManagerT.h"
#include "ObjectManagers/DiscreteElementManagerT.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "MPI_Communications/SpatialPartition.h"

#include "ElementLibrary/FiniteElement.h"

#include <limits.h>

namespace GPAC_IO
{
  static array<integer> AbaqusNodeOrdering( const std::string& elementGeometry ){

    array<integer> nodeOrdering;



    if( !elementGeometry.compare(0, 4, "CPE2") )
    {
      nodeOrdering.resize(2);
      nodeOrdering[0] = 0;
      nodeOrdering[1] = 1;
    }
    else if( !elementGeometry.compare(0, 4, "CPE3") )
    {
      nodeOrdering.resize(3);
      nodeOrdering[0] = 0;
      nodeOrdering[1] = 1;
      nodeOrdering[2] = 2;
      //      throw GPException("ElementRegionT::AllocateElementLibrary(): CPE3 unimplemented");
    }
    else if( !elementGeometry.compare(0, 4, "CPE4") )
    {
      nodeOrdering.resize(4);
      nodeOrdering[0] = 0;
      nodeOrdering[1] = 1;
      nodeOrdering[2] = 3;
      nodeOrdering[3] = 2;
    }
    else if( !elementGeometry.compare(0, 4, "C3D4") )
    {
      nodeOrdering.resize(4);
      nodeOrdering[0] = 0;
      nodeOrdering[1] = 1;
      nodeOrdering[2] = 2;
      nodeOrdering[3] = 3;
    }
    else if( !elementGeometry.compare(0, 4, "C3D8") )
    {
      nodeOrdering.resize(8);
      nodeOrdering[0] = 0;
      nodeOrdering[1] = 1;
      nodeOrdering[2] = 3;
      nodeOrdering[3] = 2;
      nodeOrdering[4] = 4;
      nodeOrdering[5] = 5;
      nodeOrdering[6] = 7;
      nodeOrdering[7] = 6;
    }
    else if( !elementGeometry.compare(0, 4, "C3D6") )
    {
      nodeOrdering.resize(8);
      nodeOrdering[0] = 0;
      nodeOrdering[1] = 1;
      nodeOrdering[2] = 3;
      nodeOrdering[3] = 2;
      nodeOrdering[4] = 4;
      nodeOrdering[5] = 5;
      nodeOrdering[6] = 4;
      nodeOrdering[7] = 5;
    }
    else if( !elementGeometry.compare(0, 4, "STRI") )
    {
      nodeOrdering.resize(3);
      nodeOrdering[0] = 0;
      nodeOrdering[1] = 1;
      nodeOrdering[2] = 2;
    }
    else if( !elementGeometry.compare(0, 3, "S4R") )
    {
      nodeOrdering.resize(4);
      nodeOrdering[0] = 0;
      nodeOrdering[1] = 1;
      nodeOrdering[2] = 2;
      nodeOrdering[3] = 3;
    }
    else if( !elementGeometry.compare(0, 4, "TRSH") )
    {
      nodeOrdering.resize(4);
      nodeOrdering[0] = 0;
      nodeOrdering[1] = 1;
      nodeOrdering[2] = 2;
    }
    return nodeOrdering;
  }






  /**
   * Constructor sets refences
   * @param pman reference to problem manager
   */
  FileManagerT::FileManagerT() :
    ensight_sequence_num(0), output_precision(14),geometryUnits(1.0),geometryMessageSize(100000)
  {
    //ts_output.resize(10);
  }

  /** Destructor */
  FileManagerT::~FileManagerT()
  {
  }

  /**
   * Set fileroot std::string
   * @param root file root
   */
  void FileManagerT::SetRoot(const char* const root)
  {
    //fileroot = root ;
    inputfilename = std::string(root) + ".xml";
    geometryfilename = std::string(root) + ".geom";
    degeometryfilename = std::string(root) + "_de.geom";
    fpgeometryfilename = std::string(root) + "_fp.geom";
  }

  /**
   * Set input filename std::string
   * @param filename input filename
   */
  void FileManagerT::SetInputFilename(const char* const filename)
  {
    inputfilename = filename;
  }

  /**
   * Set discrete element geometry filename std::string
   * @param filename discrete element geometry filename
   */
  void FileManagerT::SetDiscreteElementGeometryFilename(const char* const filename)
  {
    degeometryfilename = filename;
  }

  /**
   * Set discrete element geometry filename std::string
   * @param filename discrete element geometry filename
   */
  void FileManagerT::SetEllipsoidalDiscreteElementGeometryFilename(const char* const filename)
  {
    edegeometryfilename = filename;
  }

#ifdef SRC_EXTERNAL
  /**
   * Set fault patch geometry filename std::string
   * @param filename fault patch element geometry filename
   */
  void FileManagerT::SetFaultPatchElementGeometryFilename(const char* const filename)
  {
    fpgeometryfilename = filename;
  }
#endif
  /**
   * Set geometry filename std::string
   * @param filename geometry filename
   */
  void FileManagerT::SetGeometryFilename(const char* const filename)
  {
    geometryfilename = filename;
  }

  /** Read Input data from an Ascii file */
  void FileManagerT::ReadAsciiInput()
  {

  }

  /** Read Input data from an XML file */
  bool FileManagerT::ReadXMLInput(TICPP::HierarchicalDataNode& hdn) const
  {
    return TICPP::TinyXMLParser::Load(inputfilename.c_str(), &hdn);
  }

  bool FileManagerT::ReadXMLInput(const char* const filename, TICPP::HierarchicalDataNode& hdn) const
  {
    return TICPP::TinyXMLParser::Load(filename, &hdn);
  }

  /** Read Mesh data from an XML file */
  bool FileManagerT::ReadMeshXML(TICPP::HierarchicalDataNode* meshNode)
  {
    bool fileNameSpecified = false;

    realT units = meshNode->GetAttributeOrDefault<realT> ("units", 1.0);
    globalIndex messageSize = meshNode->GetAttributeOrDefault<realT> ("messagesize", 10000000);

    SetGeometryUnits(units);
    SetMessageSize(messageSize);

    std::string meshFileString = meshNode->GetAttributeString("file");
    if(!meshFileString.empty()) fileNameSpecified = true;
    SetGeometryFilename(meshFileString.c_str());

    meshFileString = meshNode->GetAttributeString("de_file");
    if(!meshFileString.empty()) fileNameSpecified = true;
    SetDiscreteElementGeometryFilename(meshFileString.c_str());

    meshFileString = meshNode->GetAttributeString("ellipsoid_file");
    if(!meshFileString.empty()) fileNameSpecified = true;
    SetEllipsoidalDiscreteElementGeometryFilename(meshFileString.c_str());

#ifdef SRC_EXTERNAL
    meshFileString = meshNode->GetAttributeString("faultpatch_file");
    if(!meshFileString.empty()) fileNameSpecified = true;
    SetFaultPatchElementGeometryFilename(meshFileString.c_str());
#endif
    return fileNameSpecified;
  }

  #if 1
  /**
   * @brief Converts an element type string to the number of nodes
   * @author Scott Johnson
   * @param[in] elementType type string
   * @return Number of nodes associated with the element
   */
  unsigned int FileManagerT::ElementTypeToNumberOfNodes(const std::string& elementType)
  {
    if (!elementType.compare(0, 4, "CPE2") )//line
      return 2;
    else if (!elementType.compare(0, 4, "CPE3") )//tri
      return 3;
    else if (!elementType.compare(0, 4, "CPE4") )//quad
      return 4;
    else if (!elementType.compare(0, 4, "C3D8") || !elementType.compare(0, 4, "C3D6") )//hex and prism
      return 8;
    else if (!elementType.compare(0, 4, "C3D4") )//tet
      return 4;
    else if (!elementType.compare(0, 4, "STRI") )//tet
      return 3;
    else if (!elementType.compare(0, 3, "S4R") )//quadrilateral shell
      return 4;
    else if (!elementType.compare(0, 4, "TRSH") )//triangular shell
      return 3;
    else
      throw GPException("FileManagerT::ElementTypeToNumberOfNodes: cannot parse element type " + elementType);
  }

  #if 0
  /**
   * @brief Converts an element type string to the number of nodes per face
   * @author walsh24
   * @param[in] elementType type string
   * @return Number of nodes associated with the element
   */
  unsigned int FileManagerT::ElementTypeToNumberOfNodesPerFace(const std::string& elementType)
  {
    if (!elementType.compare(0, 4, "CPE3") )//tri
      return 2;
    else if (!elementType.compare(0, 4, "CPE4") )//quad
      return 2;
    else if (!elementType.compare(0, 4, "C3D8")  || !elementType.compare(0, 4, "C3D6") )//hex
      return 4;
    else if (!elementType.compare(0, 4, "C3D4") )//tet
      return 3;
    else
      throw GPException("FileManagerT::ElementTypeToNumberOfNodesPerFace: cannot parse element type " + elementType);
  }

  /**
   * @brief Converts an element type string to the number of faces
   * @author Scott Johnson
   * @param[in] elementType type string
   * @return Number of faces associated with the element
   */
  unsigned int FileManagerT::ElementTypeToNumberOfFaces(const std::string& elementType)
  {
    if (!elementType.compare(0, 4, "CPE3") )//tri
      return 1;
    else if (!elementType.compare(0, 4, "CPE4") )//quad
      return 1;
    else if (!elementType.compare(0, 4, "C3D8")  || !elementType.compare(0, 4, "C3D6") )//hex
      return 6;
    else if (!elementType.compare(0, 4, "C3D4") )//tet
      return 4;
    else
      throw GPException("FileManagerT::ElementTypeToNumberOfFaces: cannot parse element type " + elementType);
  }


  /**
   * @brief Converts an element type string to the number of dimensions
   * @author walsh24
   * @param[in] elementType type string
   * @return Number of nodes associated with the element
   */
  unsigned int FileManagerT::ElementTypeToNumberOfDimensions(const std::string& elementType)
  {
    if (!elementType.compare(0, 4, "CPE3") )//tri
      return 2;
    else if (!elementType.compare(0, 4, "CPE4") )//quad
      return 2;
    else if (!elementType.compare(0, 4, "C3D8")  || !elementType.compare(0, 4, "C3D6") )//hex
      return 3;
    else if (!elementType.compare(0, 4, "C3D4") )//tet
      return 3;
    else
      throw GPException("FileManagerT::ElementTypeToNumberOfDimensions: cannot parse element type " + elementType);
  }
  #endif
  #endif

  /**
   * @brief Read geometry definitions
   * @author Scott Johnson
   * Reads Abaqus files for finite element (and/or) discrete element definitions
   * For the discrete element mesh file, discrete elements are defined in nodesets within the file
   * Centers of the discrete elements are calculated afterwards using an arbitrarily selected
   * node on the surface and taking the volume-weighted average of the constituent tetrathedra
   * formed by the surface faces and the arbitrary node
   * Also, reads geometry definition for ellipsoidal DE's
   * @param[in,out] domain Physical domain object where element collections is stored
   * @param[in,out] partition Spatial partitioning object
   *
   *
   * @modified author: Xiao Chen for external mesh reading using MPI
   */
  void FileManagerT::ReadMesh(PhysicalDomainT& domain, SpatialPartition& partition)
  {
    //setup variables to store FE, DE, and FPE properties, respectively
    array<FileManagerDataT*> fd;
    AbaqusFileManagerDataT fe, de, fp;
    EllipsoidFileManagerDataT ee;

    {
      fd.push_back(&fe);
      fd.back()->filename = this->geometryfilename;
      fd.push_back(&de);
      fd.back()->filename = this->degeometryfilename;
      fd.push_back(&fp);
      fd.back()->filename = this->fpgeometryfilename;
      fd.push_back(&ee);
      fd.back()->filename = this->edegeometryfilename;
    }
    const unsigned int nManagers = fd.size();
    const unsigned int iee = fd.size() - 1;

    //Read FE, DE, FPE, EDE
    //The very first time is to set up the following basic information
    //without broadcasting the node and element information to other processors

    {
      localIndex ifd = 0;
      for (array<FileManagerDataT*>::iterator it = fd.begin(); it != fd.end(); ++it, ++ifd)
      {
        if((*it)->OpenFile())
        {
          if(ifd < iee)
          {
            if (partition.m_rank == 0)
              std::cout << "Reading Mesh from " << (*it)->filename << "\n";
            AbaqusFileManagerDataT& afd = static_cast<AbaqusFileManagerDataT&>(**it);

            if (ifd == 1) ReadDiscreteAbaqusMeshA(afd);//To be fixed

            else
            {
              //Only master processor reads the data
              if (partition.m_rank == 0)
              {
                ReadAbaqusMeshA(afd);
              }
            }


            //Pack up
            realT bufSpatialMin[nsdof];
            realT bufSpatialMax[nsdof];

            //spatialMin and spatialMax
            if (partition.m_rank == 0)
            {
              for(unsigned int i = 0; i < nsdof; i++)
                {
                  bufSpatialMin[i] = afd.spatialMin(i);
                  bufSpatialMax[i] = afd.spatialMax(i);
                }
            }

            //spatialMin and spatialMax
            {
              MPI_Barrier(MPI_COMM_WORLD);
              MPI_Bcast(bufSpatialMin, nsdof, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(bufSpatialMax, nsdof, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }

            //Copy buffer content into memory
            if (partition.m_rank != 0)
            {
              //spatialMin and spatialMax
              for(unsigned int i = 0; i < nsdof; i++)
                {
                  afd.spatialMin(i) = bufSpatialMin[i];
                  afd.spatialMax(i) = bufSpatialMax[i];
                }
            }

             //Broadcast
            {
              MPI_Barrier(MPI_COMM_WORLD);
              MPI_Bcast(&(afd.maxGlobalNodeID), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
            }

            //Broadcast
            {
              //numElementRegions
              MPI_Barrier(MPI_COMM_WORLD);
              MPI_Bcast(&(afd.numElementRegions), 1, MPI_INT, 0, MPI_COMM_WORLD);
            }


            {
              bufvector bufferElemRegion;
              bufferElemRegion.Pack(afd.elementRegionTypes);
              bufferElemRegion.Pack(afd.elementRegionNames);

              int bufferLength = bufferElemRegion.size();
              MPI_Bcast(&bufferLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
              bufferElemRegion.resize( bufferLength);
              MPI_Bcast(bufferElemRegion.data(), bufferElemRegion.size(), MPI_CHAR, 0, MPI_COMM_WORLD);

              const char* pbuffer = bufferElemRegion.data();
              bufvector::Unpack( pbuffer, afd.elementRegionTypes );
              bufvector::Unpack( pbuffer, afd.elementRegionNames );
            }



            //Copy buffer content into memory
            if (partition.m_rank != 0)
            {
              for (localIndex i = 0; i < afd.numElementRegions; i++)
              {
                //do some update work
                afd.numElementsInRegion[afd.elementRegionNames[i]];
                afd.elemsInRegion[afd.elementRegionNames[i]];
                afd.maxNumElementsInRegion[afd.elementRegionNames[i]];
              }
            }


          }
          else
          {
            if (partition.m_rank == 0)
              std::cout << "Reading ellipsoids from " << (*it)->filename << std::endl;
            EllipsoidFileManagerDataT& efd = static_cast<EllipsoidFileManagerDataT&>(**it);
            ReadEllipsoidFileA(efd);
          }
        }
      }
    }

    //set the spatial partitioning at this point to facilitate the MPI broadcast
    {
      for(localIndex i = 1; i < nManagers; i++)
        fd[0]->spatialMin.SetMin(fd[i]->spatialMin);
      for(localIndex i = 1; i < nManagers; i++)
        fd[i]->spatialMin = fd[0]->spatialMin;
      for(localIndex i = 1; i < nManagers; i++)
        fd[0]->spatialMax.SetMax(fd[i]->spatialMax);
      for(localIndex i = 1; i < nManagers; i++)
        fd[i]->spatialMax = fd[0]->spatialMax;
      //All of the spatial extents have been synchronized at this point
      partition.setSizes(fd[0]->spatialMin, fd[0]->spatialMax);
    }

   //Read FE, DE, FPE, EDE
    {
      localIndex ifd = 0;
      for (array<FileManagerDataT*>::iterator itt = fd.begin(); itt != fd.end(); ++itt, ++ifd)
      {
        if((*itt)->OpenFile())
        {
          if(ifd < iee)
          {
            AbaqusFileManagerDataT& afd = static_cast<AbaqusFileManagerDataT&>(**itt);


            /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            /*-----------------------MPI-NODE-DATA-INTENSIVE-STARTED---------------------------*/
            /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

            afd.GoToSection(NODE);

            const globalIndex mpiNodeLimit = geometryMessageSize; //Define mpiNodeLimit to be 1M later on

            globalIndex sizeofBuf;

            array<real64> bufNodalPositions( (mpiNodeLimit+1)*nsdof );
            gArray1d bufNodes(mpiNodeLimit+1);

            std::map<globalIndex,R1Tensor> tempNodalPositionsMap;

            while (afd.maxGlobalNodeID > 0)
            {

              tempNodalPositionsMap.clear();
              sizeofBuf = 0;

              //Only master processor reads the data
              if (partition.m_rank == 0)
              {
                //Only master processor need to have nodalPositions to read an array of nodal coordinates
                ReadAbaqusMeshA(afd, mpiNodeLimit, sizeofBuf, tempNodalPositionsMap);
              }

              //Broadcast
              {
                //bufNodalPositionsLength
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Bcast(&sizeofBuf, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
              }

              R1Tensor nodePosition; //nodalPositions

              //Pack up
              if (partition.m_rank == 0)
              {
                globalIndex globalNodeNumber = 0;
                std::map<globalIndex,R1Tensor> :: iterator it = tempNodalPositionsMap.begin();

                for(; (globalNodeNumber < sizeofBuf)||(it != tempNodalPositionsMap.end()); globalNodeNumber++, ++it)
                {
                  bufNodes[globalNodeNumber] = it->first;
                  nodePosition = it->second;
                  for(unsigned int i = 0; i < nsdof; i++)
                    bufNodalPositions[globalNodeNumber*nsdof+i] = nodePosition(i);
                }
              }

              //Broadcast
              {
                //nodalPositions
                //to be changed to resize with MPI_NODE_LIMIT * nsdof (it depends on the size of bufNodalPositions)
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Bcast(bufNodes.data(), sizeofBuf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Bcast(bufNodalPositions.data(), sizeofBuf*nsdof, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              }

              //Update
              if (partition.m_rank != 0)
              {
                for(globalIndex globalNodeNumber = 0; globalNodeNumber < sizeofBuf; globalNodeNumber++)
                {
                  for(unsigned int i = 0; i < nsdof; i++)
                    nodePosition(i) = bufNodalPositions[globalNodeNumber*nsdof+i];
                  tempNodalPositionsMap[bufNodes[globalNodeNumber]] = nodePosition;
                }
              }


              globalIndex globalNodeNumber = 0;

              //Update potentialNodalPositionsMap based on bufNodalPositions built by tempNodalPositionsMap
              {
                //nodalPositions
                for(std::map<globalIndex,R1Tensor> :: iterator it = tempNodalPositionsMap.begin(); it != tempNodalPositionsMap.end(); ++it)
                {
                  globalNodeNumber = it->first;
                  nodePosition = it->second;

                  //one neighborhood partition beside
                  if( partition.IsCoordInPartition(nodePosition, 1) ) afd.potentialNodalPositionsMap[globalNodeNumber] = nodePosition;
                }
              }
              if(globalNodeNumber == afd.maxGlobalNodeID) break;
            }

            /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            /*----------------------MPI-NODE-DATA-INTENSIVE-FINISISHED-------------------------*/
            /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
          }
        }
      }
    }



    //read the FE, DE, and FPE files more
    {
      localIndex ifd = 0;
      for (array<FileManagerDataT*>::iterator it = fd.begin(); it != fd.end(); ++it, ++ifd)
      {
        if((*it)->exist)
        {
          //FE, DE, FPE
          if(ifd < iee)
          {
            AbaqusFileManagerDataT& afd = static_cast<AbaqusFileManagerDataT&>(**it);
            if(ifd == 0)
              ReadFiniteElementAbaqusMesh(domain, partition, afd);
            else if(ifd == 1)
              ReadDiscreteElementAbaqusMesh(domain, partition, afd);
            else
#ifdef SRC_EXTERNAL
              ReadFaultPatchElementAbaqusMesh(domain, partition, afd);
#else
              throw GPException("ReadFaultPatchElementAbaqusMesh() is not included in this build");
#endif
          }
          else
          {
            EllipsoidFileManagerDataT& efd = static_cast<EllipsoidFileManagerDataT&>(**it);
            ReadEllipsoidFile(domain, partition, efd);
          }
        }
      }
    }
  }

  /**
   * @brief Read through the finite element Abaqus mesh
   * @author Rui Wang
   * @param[in,out] domain, spatial partition
   */
  void FileManagerT::ReadMeshforMetis(PhysicalDomainT& domain, SpatialPartition& partition)
  {

    std::map<globalIndex,R1Tensor> tempNodalPositionsMap;
    std::map<globalIndex,globalIndex> resortNodeMap;
    AbaqusFileManagerDataT finiteElementFileData;

    finiteElementFileData.filename = this->geometryfilename;

    //Read the finite element file for number of nodes, nodal positions and element regions
    if(finiteElementFileData.OpenFile())
    {
      //Only master processor reads the data
      if (partition.m_rank == 0)
      {
        std::cout << "Reading Mesh from " << finiteElementFileData.filename << "\n";
        ReadAbaqusNodeRegion(finiteElementFileData,tempNodalPositionsMap,resortNodeMap);
      }
      //Broadcast
      {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&(finiteElementFileData.maxGlobalNodeID), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      }

      //Broadcast
      {
        //numElementRegions
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&(finiteElementFileData.numElementRegions), 1, MPI_INT, 0, MPI_COMM_WORLD);
      }


      {
        bufvector bufferElemRegion;
        bufferElemRegion.Pack(finiteElementFileData.elementRegionTypes);
        bufferElemRegion.Pack(finiteElementFileData.elementRegionNames);

        int bufferLength = bufferElemRegion.size();
        MPI_Bcast(&bufferLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
        bufferElemRegion.resize( bufferLength);
        MPI_Bcast(bufferElemRegion.data(), bufferElemRegion.size(), MPI_CHAR, 0, MPI_COMM_WORLD);

        const char* pbuffer = bufferElemRegion.data();
        bufvector::Unpack( pbuffer, finiteElementFileData.elementRegionTypes );
        bufvector::Unpack( pbuffer, finiteElementFileData.elementRegionNames );
      }



      //Copy buffer content into memory
      if (partition.m_rank != 0)
      {
        for (localIndex i = 0; i < finiteElementFileData.numElementRegions; i++)
        {
          //do some update work
          finiteElementFileData.numElementsInRegion[finiteElementFileData.elementRegionNames[i]];
          finiteElementFileData.elemsInRegion[finiteElementFileData.elementRegionNames[i]];
          finiteElementFileData.maxNumElementsInRegion[finiteElementFileData.elementRegionNames[i]];
        }
      }
    }

    /****Read and partition element connectivity, region and nodal position*******/
    if (finiteElementFileData.exist)
    {
      ReadAbaqusElement(domain, partition, finiteElementFileData, resortNodeMap, tempNodalPositionsMap);
    }

    /******set local to global node ID********/

    finiteElementFileData.geometry.clear();

    // allocate node objects
    domain.m_feNodeManager.resize(finiteElementFileData.numNodes*1.5);
    domain.m_feNodeManager.resize(finiteElementFileData.numNodes);

    globalIndex globalNodeNumber = 0;
    localIndex localNodeNumber = 0;
    {
      std::map<globalIndex,R1Tensor>::iterator npos = finiteElementFileData.nodalPositionsMap.begin();

      gArray1d& ltog = domain.m_feNodeManager.m_localToGlobalMap;
      array<R1Tensor>& rpos = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition>();
      for (; npos != finiteElementFileData.nodalPositionsMap.end(); ++npos)
      {

        globalNodeNumber = npos->first;
        ltog[localNodeNumber] = globalNodeNumber;
        rpos[localNodeNumber] = npos->second;

        finiteElementFileData.GlobalToLocalNodeMap[globalNodeNumber] = localNodeNumber;
        ++localNodeNumber;

      }
    }

    std::cout << "rank " << partition.m_rank << " has " << finiteElementFileData.numNodes << " of nodes." << std::endl;
    std::cout << "rank " << partition.m_rank << " has " << finiteElementFileData.numElementRegions
        << " of element regions." << std::endl;

    /********set local to global element ID************/
    lvector numElements(finiteElementFileData.numElementRegions);
    lvector numElements2(finiteElementFileData.numElementRegions);
    for (localIndex i = 0; i < finiteElementFileData.numElementRegions; ++i)
    {

      numElements[i] = finiteElementFileData.numElementsInRegion[finiteElementFileData.elementRegionNames[i]];
      std::cout << "  element region " << finiteElementFileData.elementRegionNames[i] << " is of type \""
          << finiteElementFileData.elementRegionTypes[i] << "\", and has " << numElements[i] << " number of elements."
          << std::endl;
      numElements2[i] = numElements[i] * 1.25;
    }

    // allocate element objects
    domain.m_feElementManager.resize(numElements2, finiteElementFileData.elementRegionNames, finiteElementFileData.elementRegionTypes );
    domain.m_feElementManager.resize(numElements, finiteElementFileData.elementRegionNames, finiteElementFileData.elementRegionTypes );
    // now fill the elementfToNode array
    for( localIndex count=0 ; count<finiteElementFileData.elementRegionNames.size() ; ++count )
    {
      const std::string& regionName = finiteElementFileData.elementRegionNames[count];
      ElementRegionT& elemRegion = domain.m_feElementManager.m_ElementRegions[regionName];
      elemRegion.m_maxGlobalNumber = finiteElementFileData.maxNumElementsInRegion[regionName];

      const localIndex numNodes = ElementTypeToNumberOfNodes(finiteElementFileData.elementRegionTypes[count]);
      elemRegion.m_numNodesPerElem = numNodes;
      if (elemRegion.m_numElems > 0)
      {
        const array<integer> nodeOrdering = AbaqusNodeOrdering(elemRegion.m_elementGeometryID);
        gArray1d& elemLocalToGlobal = elemRegion.m_localToGlobalMap;

        // fill element node array
        localIndex localElemIndexInRegion = 0;


        for (gArray1d::iterator elemNumber = finiteElementFileData.elemsInRegion[regionName].begin();
            elemNumber != finiteElementFileData.elemsInRegion[regionName].end(); ++elemNumber)
        {
          elemLocalToGlobal(localElemIndexInRegion) = *elemNumber;
          localIndex* elemToNodeMap = elemRegion.ElementToNodeMap(localElemIndexInRegion++);
          for(localIndex a=0 ; a<elemRegion.m_numNodesPerElem ; ++a )
          {
            elemToNodeMap[a] = finiteElementFileData.GlobalToLocalNodeMap[ finiteElementFileData.elemToNodesMap[*elemNumber][nodeOrdering[a]] ];
          }

        }
      }
    }
    /*********Generate the neighbors for each partition*****************/
    partition.AddNeighborsMetis(finiteElementFileData.neighborList);
    partition.SetRankOfNeighborNeighbors();
    partition.SetDomain(domain);

    /*********Read Node Set*********/
    ReadAbaqusNodeSet(domain, partition, finiteElementFileData, resortNodeMap);


    // close the geometry file
    finiteElementFileData.CloseFile();

    MPI_Barrier(MPI_COMM_WORLD);
    if (partition.m_rank == 0)
      std::cout << "Done Reading Finite Element Mesh.\n";


  }

  /**
   * @brief Read through the finite element Abaqus mesh for nodes and element regions
   * @author Rui Wang
   * @param[in,out] femData File data manager object, temporary nodal position map
   */
  void FileManagerT::ReadAbaqusNodeRegion(AbaqusFileManagerDataT& femData, std::map<globalIndex,R1Tensor> &tempNodalPositionsMap, std::map<globalIndex,globalIndex> &resortNodeMap)
  {
    //Reset the data structure
    femData.Reset();

    // first pass through file get node positions and set spatial extent of the mesh
    std::string inputline;
    while (femData.OK())
    {
      const bool modeChange = femData.AdvanceLine(inputline);
      if (modeChange)
      {
        if(femData.mode == ELEM)
        {
          femData.AddElementRegionLine(inputline);
        }
      }
      else if (femData.mode == NODE)
      {
        femData.AddNodeLine(inputline, geometryUnits, tempNodalPositionsMap, resortNodeMap);
      }
      else if (femData.mode == NSET)
      {
        break;
      }
    }
  }

  /**
   * @brief Read the Abaqus mesh file for finite element connectivity, partitions, send/recv element, region, node data
   * @author Rui Wang
   * @param[in,out] domain Physical domain object
   * @param[in,out] partition Spatial partition object
   * @param[in,out] femData File data manager object
   * @param[in,out] resorting map for global node ID
   */
  void FileManagerT::ReadAbaqusElement( PhysicalDomainT& domain , SpatialPartition& partition, AbaqusFileManagerDataT& femData,
                                        std::map<globalIndex,globalIndex> &resortNodeMap, std::map<globalIndex,R1Tensor> &tempNodalPositionsMap)
  {
    globalIndex numNodesPerElem;
    std::map<std::string,globalIndex> numNodesPerElemMap, numActualNodesPerElemMap;

    std::string regionName, inputline;

    std::vector<std::string> tempGEleToRegionMap;
    std::map<globalIndex,gArray1d> tempElemToNodesMap;
    std::map<globalIndex,gArray1d> tempNodeToElemsMap;

    std::map<globalIndex,gSet> tempNodeToPartsMap;
    array<gArray1d> bufRegionGElem;
    std::vector <idx_t> elementConnectVector, elementStartVector;

    std::vector<std::map<std::string,gArray1d>> tempPartRegionToGElemMap;
    std::vector<std::map<globalIndex,gArray1d>> tempPartElemToNodesMap;
    std::vector<std::map<globalIndex,R1Tensor>> tempPartNodalPositionMap;
    std::vector<gSet> tempTempPartNeighborList, tempPartNeighborList;
    tempPartRegionToGElemMap.resize(partition.m_size);
    tempPartElemToNodesMap.resize(partition.m_size);
    tempPartNodalPositionMap.resize(partition.m_size);
    tempTempPartNeighborList.resize(partition.m_size);
    tempPartNeighborList.resize(partition.m_size);

    /************* allow mixed type elements **************/

    globalIndex maxNumNodesPerElem = 0;

    if (partition.m_rank == 0)
    {
      for(std::map<std::string, std::string>::const_iterator iElementRegionNameType = femData.elementRegionNameTypes.begin(); iElementRegionNameType != femData.elementRegionNameTypes.end(); ++iElementRegionNameType)
      {
        numNodesPerElem = ElementTypeToNumberOfNodes(iElementRegionNameType->second);
        if(numNodesPerElem > maxNumNodesPerElem)
          maxNumNodesPerElem = numNodesPerElem;
      }
    }

    /************* allow mixed type elements **************/
    //Master reads element part

    if (partition.m_rank == 0)
    {
      unsigned int totalElementNodeCount = 0;
      femData.GoToSection(ELEM);
      while (femData.OK())
      {
        const bool modeChange = femData.AdvanceLine(inputline);
        if (modeChange)
        {
          if(femData.mode == ELEM)
          {
            regionName = FileManagerDataT::ExtractValue("ELSET", inputline);
            femData.maxNumElementsInRegion[regionName] = 0;
            //retrieve number of nodes for the given element type
            std::map<std::string, std::string>::const_iterator iElementRegionNameType = femData.elementRegionNameTypes.find(regionName);
            if (iElementRegionNameType == femData.elementRegionNameTypes.end())
            {
              std::cout<< regionName<<std::endl;
              throw GPException("Cannot find the requested region");
            }
            numNodesPerElem = ElementTypeToNumberOfNodes(iElementRegionNameType->second);
            numNodesPerElemMap[regionName] = numNodesPerElem;
            if (iElementRegionNameType->second == "C3D6")
            {
              numActualNodesPerElemMap[regionName] = 6;
            }
            else
            {
              numActualNodesPerElemMap[regionName] = numNodesPerElem;
            }
          }
          else
            femData.mode = UNDEF;
        }
        else if (femData.mode == ELEM)
        {
          //retrieve global element index
          globalIndex globalElemNum;
          std::istringstream linestream(inputline);
          linestream >> globalElemNum;
          //resort the global element numbering
          globalElemNum = femData.totalNumElem;
          if ((globalElemNum+1) > femData.maxNumElementsInRegion[regionName])
            femData.maxNumElementsInRegion[regionName] = globalElemNum+1;


          femData.totalNumElem = femData.totalNumElem + 1;


          tempGEleToRegionMap.push_back(regionName);

          //Read the nodal connectivity of each element
          for (globalIndex a = 0; a < numActualNodesPerElemMap[regionName]; ++a)
          {
            globalIndex globalNodeNum;
            linestream >> globalNodeNum;
            globalNodeNum=resortNodeMap[globalNodeNum];
            tempElemToNodesMap[globalElemNum].push_back(globalNodeNum);
            tempNodeToElemsMap[globalNodeNum].push_back(globalElemNum);
            //Metis style vector
            elementConnectVector.push_back(globalNodeNum);
          }

          //Special treatment for prisms.
          if (numActualNodesPerElemMap[regionName] == 6 && numNodesPerElemMap[regionName] == 8)
          {
            gArray1d hexNodeMapForPrism;
            hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][4]);
            hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][1]);
            hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][0]);
            hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][3]);
            hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][5]);
            hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][2]);
            hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][5]);
            hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][2]);
            tempElemToNodesMap[globalElemNum].clear();
            tempElemToNodesMap[globalElemNum] = hexNodeMapForPrism;
          }

          //Generate the vectors that Metis needs
          elementStartVector.push_back(totalElementNodeCount);
          totalElementNodeCount = totalElementNodeCount + numActualNodesPerElemMap[regionName];

        }
      }
      elementStartVector.push_back(totalElementNodeCount);

      /***********METIS PARTITION**************/
      idx_t nn = femData.maxGlobalNodeID + 1;
      idx_t ne = femData.totalNumElem;
      idx_t *vwgt=NULL, *vsize=NULL;
      idx_t ncommon=1;
      real_t *tpwgts = NULL;
      idx_t *options = NULL;
      idx_t objval;

      array<idx_t> elePart(ne);
      array<idx_t> nodePart(nn);
//      elePart.resize(ne);
//      nodePart.resize(nn);


      if(partition.m_size != 1)
      {
        METIS_PartMeshDual(&ne, &nn, elementStartVector.data(), elementConnectVector.data(),
                           vwgt, vsize, &ncommon, &partition.m_size, tpwgts, options, &objval, elePart.data(), nodePart.data());
      }
      else
      {
        for(int i = 0; i < ne; i++)
        {
          elePart[i]=0;
        }
      }



      /***********DISTRIBUTE PARTITIONS**************/
      std::map<globalIndex,gArray1d>::iterator itElem = tempElemToNodesMap.begin();
      for (; itElem != tempElemToNodesMap.end(); itElem++)
      {
        tempPartElemToNodesMap[elePart[itElem->first]][itElem->first] = itElem->second;
        tempPartRegionToGElemMap[elePart[itElem->first]][tempGEleToRegionMap[itElem->first]].push_back(itElem->first);
        unsigned int tempNodeNum = (itElem->second).size();
        for (globalIndex a = 0; a < tempNodeNum; ++a)
        {
          tempPartNodalPositionMap[elePart[itElem->first]][(itElem->second)[a]] = tempNodalPositionsMap[(itElem->second)[a]];
          tempNodeToPartsMap[(itElem->second)[a]].insert(elePart[itElem->first]);
        }

      }
      /*******FIND NEIGHBORS************/
      //R.W. note: the way that the neighbors are set right now will only guarantee to be correct for 1 layer of ghost
      if(partition.m_size > 1)
      {
        std::map<globalIndex,gSet>::iterator itNode = tempNodeToPartsMap.begin();
        for (; itNode != tempNodeToPartsMap.end(); itNode++)
        {
          gSet::iterator itPart1 = tempNodeToPartsMap[itNode->first].begin();
          for (; itPart1 != tempNodeToPartsMap[itNode->first].end(); itPart1++)
          {
            gSet::iterator itPart2 = tempNodeToPartsMap[itNode->first].begin();
            for (; itPart2 != tempNodeToPartsMap[itNode->first].end(); itPart2++)
            {
              if (*itPart1 != *itPart2)
              {
                tempTempPartNeighborList[*itPart1].insert(*itPart2);
              }
            }
          }
        }

        for (int itNeighborRank = 0; itNeighborRank < partition.m_size; itNeighborRank++)
        {

          gSet::iterator itNeighbor = tempTempPartNeighborList[itNeighborRank].begin();
          for (; itNeighbor != tempTempPartNeighborList[itNeighborRank].end(); itNeighbor++)
          {
            tempPartNeighborList[itNeighborRank].insert(*itNeighbor);
            tempPartNeighborList[itNeighborRank].insert(tempTempPartNeighborList[*itNeighbor].begin(), tempTempPartNeighborList[*itNeighbor].end());
          }
        }
        for (int itNeighborRank = 0; itNeighborRank < partition.m_size; itNeighborRank++)
        {
          tempPartNeighborList[itNeighborRank].erase(tempPartNeighborList[itNeighborRank].find(itNeighborRank));
        }
      }
    }



    gArray1d bufMaxNumElementsInRegion(femData.numElementRegions);
    if (partition.m_rank == 0)
    {
      std::map<std::string,globalIndex>::iterator it = femData.maxNumElementsInRegion.begin();
      localIndex i = 0;

      for (; (i < femData.numElementRegions)||(it!=femData.maxNumElementsInRegion.end()) ; i++,++it)
      {
        bufMaxNumElementsInRegion[i] = it->second;
      }
    }

    //Broadcast
    {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(bufMaxNumElementsInRegion.data(), femData.numElementRegions, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    }

    if (partition.m_rank != 0)
    {
      std::map<std::string,globalIndex>::iterator it = femData.maxNumElementsInRegion.begin();
      localIndex i = 0;

      for (; i < (femData.numElementRegions)||(it!=femData.maxNumElementsInRegion.end()) ; i++,++it)
      {
        it->second = bufMaxNumElementsInRegion[i];
      }
    }

    /*****MPI Send and Recv********/
    {
      MPI_Status stat;
      bufvector bufferElemNodeRegionMap;
      int bufferLength;

      if(partition.m_rank == 0)
      {
        femData.elemToNodesMap = tempPartElemToNodesMap[partition.m_rank];
        femData.elemsInRegion = tempPartRegionToGElemMap[partition.m_rank];
        femData.nodalPositionsMap = tempPartNodalPositionMap[partition.m_rank];
        femData.neighborList = tempPartNeighborList[partition.m_rank];
        array<string>::const_iterator iRegionName = femData.elementRegionNames.begin();
        for (; iRegionName != femData.elementRegionNames.end(); ++iRegionName)
        {
          femData.numElementsInRegion[*iRegionName] = femData.elemsInRegion[*iRegionName].size();
        }
      }
      // Send data from rank 0
      if (partition.m_rank == 0 && partition.m_size > 1)
      {

        for (int iRank = 1; iRank < partition.m_size; iRank++)
        {
          bufferElemNodeRegionMap.clear();
          bufferElemNodeRegionMap.Pack(tempPartElemToNodesMap[iRank]);
          bufferElemNodeRegionMap.Pack(tempPartRegionToGElemMap[iRank]);
          bufferElemNodeRegionMap.Pack(tempPartNodalPositionMap[iRank]);
          bufferElemNodeRegionMap.Pack(tempPartNeighborList[iRank]);
          int bufferLengthSend = bufferElemNodeRegionMap.size();
          MPI_Send(&bufferLengthSend, 1, MPI_UNSIGNED_LONG_LONG, iRank, iRank, MPI_COMM_WORLD);
          MPI_Send(bufferElemNodeRegionMap.data(), bufferLengthSend, MPI_CHAR, iRank, iRank, MPI_COMM_WORLD);

        }

      }
      else if (partition.m_rank != 0)
      {
        MPI_Recv(&bufferLength, 1, MPI_UNSIGNED_LONG_LONG, 0, partition.m_rank, MPI_COMM_WORLD, &stat);
        bufferElemNodeRegionMap.resize(bufferLength);
        MPI_Recv(bufferElemNodeRegionMap.data(), bufferElemNodeRegionMap.size(), MPI_CHAR, 0, partition.m_rank, MPI_COMM_WORLD, &stat);
        const char* pbuffer = bufferElemNodeRegionMap.data();
        bufvector::Unpack( pbuffer, femData.elemToNodesMap);
        bufvector::Unpack( pbuffer, femData.elemsInRegion);
        bufvector::Unpack( pbuffer, femData.nodalPositionsMap);
        bufvector::Unpack( pbuffer, femData.neighborList);
        array<string>::const_iterator iRegionName = femData.elementRegionNames.begin();
        for (; iRegionName != femData.elementRegionNames.end(); ++iRegionName)
        {
          femData.numElementsInRegion[*iRegionName] = femData.elemsInRegion[*iRegionName].size();
        }
      }
    }

    SetNumNodes(femData);

  }

  /**
   * @brief Read the Abaqus mesh file for nodesets
   * @author Rui Wang
   * @param[in,out] domain Physical domain object
   * @param[in,out] partition Spatial partition object
   * @param[in,out] femData File data manager object
   */
  void FileManagerT::ReadAbaqusNodeSet( PhysicalDomainT& domain , SpatialPartition& partition, AbaqusFileManagerDataT& femData, std::map<globalIndex,globalIndex> &resortNodeMap)
  {
    std::string inputline;
    std::string setName;
    std::map<std::string,gSet> mymap;


    int numSet;

    // now read the nodeset data
    if (partition.m_rank == 0)
    {
      femData.GoToSection(NSET);
      while (femData.OK())
      {
        //Pack up set name
        const bool modeChange = femData.AdvanceLine(inputline);
        if (modeChange && femData.mode == NSET)
        {
          setName = FileManagerDataT::ExtractValue("NSET", inputline);
          domain.m_feNodeManager.m_Sets[setName];
          mymap[setName];
        }
        else
        {
          if (modeChange)
          {
            femData.mode = UNDEF;
          }
          else if (femData.mode == NSET)
          {
            std::istringstream linestream(inputline);
            localIndex lnodeNum;
            globalIndex gnodeNum = GLOBALINDEX_MAX;
            while (linestream >> gnodeNum)
            {
              gnodeNum=resortNodeMap[gnodeNum];
              if (gnodeNum != GLOBALINDEX_MAX)
              {
                mymap[setName].insert(gnodeNum);//pack up data sets
                lnodeNum = femData.GlobalToLocalNodeMap[gnodeNum];
                if (lnodeNum != GLOBALINDEX_MAX)
                  domain.m_feNodeManager.m_Sets[setName].insert(lnodeNum);
              }
            }
          }
        }
      }
    }

    if (partition.m_rank == 0)
      numSet = mymap.size();

    //Broadcast
    {
      //numSet
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&numSet, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    {
      bufvector bufferSets;
      bufferSets.Pack(mymap);

      int bufferLength = bufferSets.size();
      MPI_Bcast(&bufferLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
      bufferSets.resize( bufferLength);
      MPI_Bcast(bufferSets.data(), bufferSets.size(), MPI_CHAR, 0, MPI_COMM_WORLD);

      const char* pbuffer = bufferSets.data();
      bufvector::Unpack( pbuffer, mymap );
    }

    //Update
    if (partition.m_rank != 0)
    {
      for ( std::map<std::string,gSet>::const_iterator i = mymap.begin() ; i != mymap.end() ; ++i )
      {
        setName = i->first;
        const gSet& globalSet = i->second;
        domain.m_feNodeManager.m_Sets[setName];
        for ( gSet::const_iterator gNode=globalSet.begin() ; gNode!=globalSet.end() ; ++gNode)
        {
          if ( *gNode != GLOBALINDEX_MAX)
          {
            localIndex lnodeNum = femData.GlobalToLocalNodeMap[*gNode];
            if (lnodeNum != GLOBALINDEX_MAX)
              domain.m_feNodeManager.m_Sets[setName].insert(lnodeNum);
          }
        }
      }
    }
  }
  /**
   * @brief Read the ellipsoidal discrete element file before spatial partitioning
   * @author Scott Johnson
   * @param[in,out] fd File data manager object
   */
  void FileManagerT::ReadEllipsoidFileA(EllipsoidFileManagerDataT& fd)
  {
    //Reset the data structure
    fd.Reset();

    // first pass through file get node positions and set spatial extent of the mesh
    std::string inputline;
    while (fd.OK())
    {
      fd.AdvanceLine(inputline);
      if (inputline.compare(0, 1, "*") != 0 && inputline.size() != 0)
      {
        fd.AddNodalPositionLine(inputline, geometryUnits);
      }
    }
  }

  /**
   * @brief Read the ellipsoid definition file
   * @author Scott Johnson
   * @param[in,out] partition Spatial partition object
   * @param[in,out] fd File data manager object
   */
  void FileManagerT::ReadEllipsoidFileB(SpatialPartition& partition, EllipsoidFileManagerDataT& fd)
  {

    std::string inputline;

    // second pass through file get node positions and set spatial extent of the mesh
    fd.numNodes = 0; // Number of nodes
    fd.geometry.clear();
    fd.geometry.seekg (0, std::ios::beg);
    while (fd.OK())
    {
      fd.AdvanceLine(inputline);
      if (inputline.compare(0, 1, "*") != 0 && inputline.size() != 0)
      {
        globalIndex globalNodeNumber;
        R1Tensor position;
        {
          EllipsoidFileManagerDataT::ReadLine(inputline, globalNodeNumber, position);
          position = fd.nodalPositions[globalNodeNumber];
        }
        if(partition.IsCoordInPartition(position))
        {
          fd.isNodeInDomain[globalNodeNumber] = 1;
          ++fd.numNodes;
        }
        else
        {
          fd.isNodeInDomain[globalNodeNumber] = 0;
        }
      }
    }
    fd.GlobalToLocalNodeMap.resize(fd.maxGlobalNodeID + 1, GLOBALINDEX_MAX);
  }

  /**
   * @brief Read the ellipsoid discrete element file
   * @author Scott Johnson
   * @param[in,out] domain Physical domain object
   * @param[in,out] partition Spatial partition object
   * @param[in,out] fd File data manager object
   */
  void FileManagerT::ReadEllipsoidFile(PhysicalDomainT& domain, SpatialPartition& partition, EllipsoidFileManagerDataT& fd)
  {
    ReadEllipsoidFileB(partition, fd);

    std::string inputline;

    // allocate ellipsoid objects
    domain.m_ellipsoidalDiscreteElementManager.resize(fd.numNodes);

    array<R1Tensor>& pradii = domain.m_ellipsoidalDiscreteElementManager.GetFieldData<R1Tensor>("principalRadii");
    array<R1Tensor>& rotationAxis = domain.m_ellipsoidalDiscreteElementManager.GetFieldData<FieldInfo::rotationAxis> ();//qx,qy,qz
    array<real64>& rotationMagnitude = domain.m_ellipsoidalDiscreteElementManager.GetFieldData<FieldInfo::rotationMagnitude> ();//qw
    array<R1Tensor>& refpos = domain.m_ellipsoidalDiscreteElementManager.GetFieldData<FieldInfo::referencePosition> ();//x,y,z
    array<R1Tensor>& curpos = domain.m_ellipsoidalDiscreteElementManager.GetFieldData<FieldInfo::currentPosition> ();//x,y,z

    // third pass through file set node positions, rotations, and radii
    fd.geometry.clear();
    fd.geometry.seekg (0, std::ios::beg);
    localIndex localNodeNumber = 0;
    while (fd.OK())
    {
      const bool modeChange = fd.AdvanceLine(inputline);
      if (!modeChange)
      {
        globalIndex globalNodeNumber;
        R1Tensor position, principalRadii;
        R1TensorT<4> rotation;
        EllipsoidFileManagerDataT::ReadLine(inputline, globalNodeNumber, position, principalRadii, rotation);
        if(fd.isNodeInDomain[globalNodeNumber])
        {
          position = fd.nodalPositions[globalNodeNumber]; // rescale to problem units
          principalRadii *= geometryUnits;//rescale to problem units

          realT rmag = rotation(0);
          R1Tensor raxis;
          {
            raxis(0) = rotation(1);
            raxis(1) = rotation(2);
            raxis(2) = rotation(3);
          }

          const realT mag = sqrt(rmag*rmag + Dot(raxis,raxis) );
          if( mag > 0.0 )
          {
            rmag /= mag;
            raxis /= mag;
          }

          domain.m_ellipsoidalDiscreteElementManager.m_globalToLocalMap[globalNodeNumber] = localNodeNumber;
          domain.m_ellipsoidalDiscreteElementManager.m_localToGlobalMap[localNodeNumber] = globalNodeNumber;

          pradii[localNodeNumber] = principalRadii;
          refpos[localNodeNumber] = position;
          curpos[localNodeNumber] = position;
          rotationMagnitude[localNodeNumber] = rmag;
          rotationAxis[localNodeNumber] = raxis;

          {
            std::string setName;
            std::stringstream out;
            out << globalNodeNumber;
            setName = out.str();
            domain.m_ellipsoidalDiscreteElementManager.m_Sets[setName].insert(localNodeNumber);
          }

          fd.GlobalToLocalNodeMap[globalNodeNumber] = localNodeNumber++;
        }
        else
        {
          {
            std::string setName;
            std::stringstream out;
            out << globalNodeNumber;
            setName = out.str();
            domain.m_ellipsoidalDiscreteElementManager.m_Sets[setName];
          }
        }
      }
    }

    // close the geometry file
    fd.CloseFile();

    if (partition.m_rank == 0)
      std::cout << "Done Reading Ellipsoid Discrete Element Mesh.\n";
  }


  bool FileManagerT::ExtractElementRegionFromAbaqus( PhysicalDomainT& domain, const PartitionBase& partition )
  {
    bool rval = false;


    bufvector buffer;
    array<string> inputLines;

    AbaqusFileManagerDataT fd;
    fd.filename = this->geometryfilename;

    if( partition.m_rank == 0 )
    {
      if(fd.OpenFile())
      {
        // first pass through file get node positions and set spatial extent of the mesh
        std::string inputline;
        while (fd.OK())
        {
          const bool modeChange = fd.AdvanceLine(inputline);
          if(modeChange)
          {
            if (fd.mode == ELEM)
            {
              inputLines.push_back( inputline );
            }
            else if (fd.mode == NSET)
            {
              break;
            }
          }
        }
      }
      if( !(inputLines.empty()) )
      {
        buffer.Pack( inputLines );
      }
    }


    {

      unsigned long long bufferLength = buffer.size();
      MPI_Bcast(&bufferLength, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      if( bufferLength>0 )
      {
        rval = true;

        buffer.resize( bufferLength);

        MPI_Bcast(buffer.data(), buffer.size(), MPI_CHAR, 0, MPI_COMM_WORLD);

        const char* pbuffer = buffer.data();
        bufvector::Unpack( pbuffer, inputLines );

        for( array<string>::iterator inputLine=inputLines.begin() ; inputLine!=inputLines.end() ; ++inputLine )
        {
          fd.AddElementRegionLine(*inputLine);
        }
        // allocate element objects
        lvector temp2;
        temp2.resize(fd.elementRegionNames.size());
        domain.m_feElementManager.resize(temp2, fd.elementRegionNames, fd.elementRegionTypes );
      }
    }


    return rval;
  }

  /**
   * @brief Read the finite element Abaqus mesh before spatial partitioning
   * @author Scott Johnson
   * @param[in,out] fd File data manager object
   */
  void FileManagerT::ReadAbaqusMeshA(AbaqusFileManagerDataT& fd)
  {
    //Reset the data structure
    fd.Reset();

    // first pass through file get node positions and set spatial extent of the mesh
    std::string inputline;
    while (fd.OK())
    {
      const bool modeChange = fd.AdvanceLine(inputline);
      if (modeChange)
      {
        if(fd.mode == ELEM)
        {
          fd.AddElementRegionLine(inputline);
        }
      }
      else if (fd.mode == NODE)
      {
        fd.AddNodalPositionLine(inputline, geometryUnits);
      }
      else if (fd.mode == NSET)
      {
        break;
      }
    }
  }


  void FileManagerT::ReadAbaqusMeshA(AbaqusFileManagerDataT& fd, const globalIndex mpiNodeLimit, globalIndex &countGlobalNodeNum, std::map<globalIndex,R1Tensor> &tempNodalPositionsMap)
  {

    // first pass through file get node positions and set spatial extent of the mesh
    std::string inputline;

    while (fd.OK())
    {
      const bool modeChange = fd.AdvanceLine(inputline);
      if (modeChange)
      {
        ;
      }
      else if (fd.mode == NODE)
      {
        fd.AddNodalPositionLine(inputline, geometryUnits, tempNodalPositionsMap);
        countGlobalNodeNum++;
        if(countGlobalNodeNum == mpiNodeLimit) break;

      }
      else if (fd.mode == NSET)
      {
        break;
      }
    }
  }

  void FileManagerT::ReadDiscreteAbaqusMeshA(AbaqusFileManagerDataT& fd)
  {
    //Reset the data structure
    fd.Reset();

    // first pass through file get node positions and set spatial extent of the mesh
    std::string inputline;
    while (fd.OK())
    {
      const bool modeChange = fd.AdvanceLine(inputline);
      if (modeChange)
      {
        if(fd.mode == ELEM)
        {
          fd.AddElementRegionLine(inputline);
        }
      }
      else if (fd.mode == NODE)
      {
        fd.AddNodalPositionLineDiscrete(inputline, geometryUnits);
      }
      else if (fd.mode == NSET)
      {
        break;
      }
    }
  }

  void FileManagerT::SetNumNodes(FileManagerDataT& fd)
  {
    // find number of nodes in the computational domain
    fd.numNodes = fd.nodalPositionsMap.size();

    fd.GlobalToLocalNodeMap.resize(fd.maxGlobalNodeID + 1, GLOBALINDEX_MAX);
  }

  void FileManagerT::SetNumNodesDiscrete(FileManagerDataT& fd)
  {
    // find number of nodes in the computational domain
    fd.numNodes = 0;
    for (array<integer>::const_iterator i = fd.isNodeInDomain.begin(); i != fd.isNodeInDomain.end(); ++i)
    {
      if (*i == 1)
        ++fd.numNodes;
    }

    fd.GlobalToLocalNodeMap.resize(fd.maxGlobalNodeID + 1, GLOBALINDEX_MAX);
  }

  void FileManagerT::ReadAbaqusMeshDEM_PartitionDEM(SpatialPartition& partition, AbaqusFileManagerDataT& fd, std::map<std::string, gSet>& dems)
  {
    std::string inputline, setName;
    globalIndex nodeIndex;

    //track center of DEM
    bool isDEM = false;
    R1Tensor center(0.0);
    localIndex npts = 0;
    gArray1d indices;

    fd.GoToSection(NSET);
    while (fd.OK())
    {
      const bool modeChange = fd.AdvanceLine(inputline);
      if(modeChange)
      {
        //first, do we need to determine partition ...
        if(npts > 0)
        {
          center /= npts;
          if(partition.IsCoordInPartition(center))
          {
            dems[setName].insert(indices.begin(), indices.end());
            for(gArray1d::const_iterator it = indices.begin(); it != indices.end(); ++it)
              fd.isNodeInDomain[*it] = 1;
          }
          npts = 0;
          center = 0.0;
          indices.clear();
          isDEM = false;
        }
        //second, for the mode change, do we need to set the name
        if (fd.mode == NSET)
        {
          setName = FileManagerDataT::ExtractValue("NSET", inputline);
          //---DE's are delimited by the reserved "de" keyword---
          isDEM = setName.length() > 1 && !setName.compare(0, 2, "de");

          npts = 0;
          center = 0.0;
        }
      }
      else if(fd.mode == NSET && isDEM)
      {
        std::istringstream iss(inputline);
        do
        {
          iss >> nodeIndex;
          indices.push_back(nodeIndex);
          center += fd.nodalPositions[nodeIndex];
          ++npts;
        } while (iss);
      }
    }
    if(npts > 0)
    {
      center /= npts;
      if(partition.IsCoordInPartition(center))
      {
          dems[setName].insert(indices.begin(), indices.end());
          for(gArray1d::const_iterator it = indices.begin(); it != indices.end(); ++it)
            fd.isNodeInDomain[*it] = 1;
      }
    }

    //set the number of nodes now ...
    SetNumNodesDiscrete(fd);
  }


  /**
   * @brief Read the Abaqus mesh file
   * @author Scott Johnson
   * @param[in,out] partition Spatial partition object
   * @param[in,out] fd File data manager object
   */
  void FileManagerT::ReadDiscreteElementAbaqusMesh(PhysicalDomainT& domain,
                                                   SpatialPartition& partition,
                                                   AbaqusFileManagerDataT& fd)
  {
    std::map<std::string, gSet> dems;

    //set the DE's and nodes that we own (incl. numNodes)
    ReadAbaqusMeshDEM_PartitionDEM(partition, fd, dems);

    // assign the nodes to the DEM Manager
    localIndex localNodeIndex = 0;
    {
      domain.m_discreteElementManager.m_nodeManager->resize(fd.numNodes);

      array<R1Tensor>& rpos = domain.m_discreteElementManager.m_nodeManager->GetFieldData<FieldInfo::referencePosition>();
      array<R1Tensor>& cpos = domain.m_discreteElementManager.m_nodeManager->GetFieldData<FieldInfo::currentPosition>();

      globalIndex globalNodeNumber = 0;
      array<R1Tensor>::iterator npos = fd.nodalPositions.begin();
      for (array<integer>::const_iterator it = fd.isNodeInDomain.begin();
          it != fd.isNodeInDomain.end(); ++it, ++npos, ++globalNodeNumber)
      {
        if (*it == 1)
        {
          //node: localToGlobal
          domain.m_discreteElementManager.m_nodeManager->m_localToGlobalMap[localNodeIndex] = globalNodeNumber;
          //node: globalToLocal
          domain.m_discreteElementManager.m_nodeManager->m_globalToLocalMap[globalNodeNumber] = localNodeIndex;
          //node: position
          rpos[localNodeIndex] = (*npos);
          cpos[localNodeIndex] = (*npos);

          fd.GlobalToLocalNodeMap[globalNodeNumber] = localNodeIndex;
          ++localNodeIndex;
        }
      }
    }

    // set the node maps
    {
      domain.m_discreteElementManager.resize(dems.size());

      localIndex demIndex = 0;
      lArray1d& nodeToDE = domain.m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::demIndex>();
      for(std::map<std::string, gSet>::const_iterator it = dems.begin(); it != dems.end(); ++it, ++demIndex)
      {
        domain.m_discreteElementManager.m_Sets[it->first].insert(demIndex);
        for(gSet::const_iterator gi = it->second.begin(); gi != it->second.end(); ++gi)
        {
          localNodeIndex = fd.GlobalToLocalNodeMap[*gi];
          domain.m_discreteElementManager.m_discreteElementToExternalNodesMap[demIndex].push_back(localNodeIndex);
          nodeToDE[localNodeIndex] = demIndex;
        }
      }
    }

    //---------------------------------------------------------------
    //1) FILL FACES

    // second pass through file determine what elements are in the domain (based on nodes)
    localIndex numElements = 0; //number of total faces which is also the number of external faces
    {
      std::string regionName, inputline;
      fd.GoToSection(ELEM);
      while (fd.OK())
      {
        const bool modeChange = fd.AdvanceLine(inputline);
        if (modeChange)
        {
          if(fd.mode == ELEM)
            regionName = FileManagerDataT::ExtractValue("ELSET", inputline);
          else
            fd.mode = UNDEF;
        }
        else if (fd.mode == ELEM)
        {
          //retrive global element index
          globalIndex globalElemNum, globalNodeNum;
          std::istringstream linestream(inputline);
          linestream >> globalElemNum;

          //add element to the region map
          fd.elemsInRegion[regionName].push_back(globalElemNum);

          //retrieve number of nodes for the given element type
          std::map<std::string, std::string>::const_iterator iernt = fd.elementRegionNameTypes.find(regionName);
          if (iernt == fd.elementRegionNameTypes.end())
            throw GPException("Cannot find the requested region");
          const localIndex numNodesPerElem = ElementTypeToNumberOfNodes(iernt->second);

          //update maximum and array sizes if necessary
          if (globalElemNum > fd.maxGlobalElemID)
          {
            fd.maxGlobalElemID = globalElemNum;
            fd.isElemInDomain.resize(fd.maxGlobalElemID + 1);
            fd.elemToNodes.resize2(fd.maxGlobalElemID + 1, numNodesPerElem);
          }

          linestream >> globalNodeNum;
          fd.isElemInDomain[globalElemNum] = fd.isNodeInDomain[globalNodeNum];
          if (fd.isElemInDomain[globalElemNum] == 1)
          {
            ++fd.numElementsInRegion[regionName];
            ++numElements;
          }
        }
      }
    }

    // allocate external face objects (note: nodeManager has already been allocated)
    domain.m_discreteElementManager.m_faceManager->resize(numElements);

    //---------------------------------------------------------------
    //2) Now assign those faces to the domain.m_discreteElementManager
    {
      std::string regionName, inputline;
      localIndex localElemNum = 0;
      globalIndex globalElemNum, globalNodeNum;
      fd.GoToSection(ELEM);
      while (fd.OK())
      {
        const bool modeChange = fd.AdvanceLine(inputline);
        if (modeChange)
        {
          if(fd.mode == ELEM)
            regionName = FileManagerDataT::ExtractValue("ELSET", inputline);
          else
            fd.mode = UNDEF;
        }
        else if (fd.mode == ELEM)
        {
          //retrive global element index
          std::istringstream linestream(inputline);
          linestream >> globalElemNum;
          if (fd.isElemInDomain[globalElemNum] == 1)
          {
            //retrieve number of nodes for the given element type
            std::map<std::string, std::string>::const_iterator iernt = fd.elementRegionNameTypes.find(regionName);
            if (iernt == fd.elementRegionNameTypes.end())
              throw GPException("Cannot find the requested region");
            const localIndex numNodesPerElem = ElementTypeToNumberOfNodes(iernt->second);

            domain.m_discreteElementManager.m_faceManager->m_toNodesRelation[localElemNum].clear();
            for(localIndex i = 0; i < numNodesPerElem; i++)
            {
              linestream >> globalNodeNum;
              localNodeIndex = fd.GlobalToLocalNodeMap[globalNodeNum];
              domain.m_discreteElementManager.m_faceManager->m_toNodesRelation[localElemNum].push_back(
                localNodeIndex);
            }
            ++localElemNum;
          }
        }
      }
    }

    // close the geometry file
    fd.CloseFile();

    //---------------------------------------------------------------
    //3) Last map to fill: DE to face map and vice versa
    {
      const lArray1d& nodeToDE = domain.m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::demIndex>();
      lArray1d& faceToDE = domain.m_discreteElementSurfaceFaces.GetFieldData<FieldInfo::demIndex>();
      for (localIndex localFaceIndex = 0; localFaceIndex < domain.m_discreteElementSurfaceFaces.DataLengths();
          ++localFaceIndex)
      {
        localNodeIndex = domain.m_discreteElementSurfaceFaces.m_toNodesRelation[localFaceIndex][0];
        const localIndex idem = nodeToDE[localNodeIndex];
        domain.m_discreteElementManager.m_discreteElementToExternalFacesMap[idem].push_back(localFaceIndex);
        faceToDE(localFaceIndex) = idem;
      }
    }

    //---------------------------------------------------------------
    //4) CALCULATE DE CENTROIDS, ROTATIONS, AND ALSO SET RELATIVE POSITION
    domain.m_discreteElementManager.RecalculatePhysicalProperties();

    if (partition.m_rank == 0)
      std::cout << "Done Reading Discrete Element Mesh.\n";

    //domain.Initialize();
  }



  /**
   * @brief Read the Abaqus mesh file and determine which are in the partition based on element center
   * @author Scott Johnson, Eric Herbold
   * @param[in,out] partition Spatial partition object
   * @param[in,out] fd File data manager object
   */
  void FileManagerT::ReadAbaqusMeshB(SpatialPartition& partition, AbaqusFileManagerDataT& fd)
  {
    // second pass through file determine what elements and nodes are in the domain
    globalIndex numNodesPerElem;
    std::map<std::string,globalIndex> numNodesPerElemMap, numActualNodesPerElemMap;

    std::string regionName, inputline;

    std::map<std::string,gSet> tempRegionGElemMap;
    std::map<globalIndex,gArray1d> tempElemToNodesMap;
    array<gArray1d> bufRegionGElem;

    /************* allow mixed type elements **************/

    globalIndex maxNumNodesPerElem = 0;

    if (partition.m_rank == 0)
    {
      for(std::map<std::string, std::string>::const_iterator iernt = fd.elementRegionNameTypes.begin(); iernt != fd.elementRegionNameTypes.end(); ++iernt)
      {
        numNodesPerElem = ElementTypeToNumberOfNodes(iernt->second);
        if(numNodesPerElem > maxNumNodesPerElem)
          maxNumNodesPerElem = numNodesPerElem;
      }
    }

    /************* allow mixed type elements **************/

    //Read element part the very first time to know how many elements and their regions and their names then broadcast
    unsigned int ii = 0;
    if (partition.m_rank == 0)
    {
      fd.GoToSection(ELEM);
      while (fd.OK())
      {
        const bool modeChange = fd.AdvanceLine(inputline);
        if (modeChange)
        {
          if(fd.mode == ELEM)
          {
            regionName = FileManagerDataT::ExtractValue("ELSET", inputline);
            ii = 0;
            fd.maxNumElementsInRegion[regionName] = 0;
          }
          else
            fd.mode = UNDEF;
        }
        else if (fd.mode == ELEM)
        {
          //retrieve global element index
          globalIndex globalElemNum;
          std::istringstream linestream(inputline);
          linestream >> globalElemNum;

          //add element to the region map
          if (fd.maxNumElementsInRegion[regionName] < globalElemNum)
          {
            fd.maxNumElementsInRegion[regionName] = globalElemNum;
            ii++;
            fd.numElementsInRegion[regionName] = ii;
          }

          //retrieve number of nodes for the given element type
          std::map<std::string, std::string>::const_iterator iernt = fd.elementRegionNameTypes.find(regionName);
          if (iernt == fd.elementRegionNameTypes.end())
          {
            throw GPException("Cannot find the requested region");
          }
          numNodesPerElem = ElementTypeToNumberOfNodes(iernt->second);
          numNodesPerElemMap[regionName] = numNodesPerElem;
          if (iernt->second == "C3D6")
          {
            numActualNodesPerElemMap[regionName] = 6;
          }
          else
          {
            numActualNodesPerElemMap[regionName] = numNodesPerElem;
          }
          //update maximum and array sizes if necessary
          if (globalElemNum > fd.maxGlobalElemID)
            fd.maxGlobalElemID = globalElemNum;
        }
      }
    }

    //Broadcast
    {
      MPI_Barrier(MPI_COMM_WORLD);

      //numNodesPerElem
      {
        bufvector bufferNumNodesPerElemMap;
        bufferNumNodesPerElemMap.Pack(numNodesPerElemMap);

        int bufferLength = bufferNumNodesPerElemMap.size();
        MPI_Bcast(&bufferLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
        bufferNumNodesPerElemMap.resize(bufferLength);
        MPI_Bcast(bufferNumNodesPerElemMap.data(), bufferNumNodesPerElemMap.size(), MPI_CHAR, 0, MPI_COMM_WORLD);

        const char* pbuffer = bufferNumNodesPerElemMap.data();
        bufvector::Unpack( pbuffer, numNodesPerElemMap );
      }

      MPI_Bcast(&maxNumNodesPerElem, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      //maxGlobalElemID
      MPI_Bcast(&(fd.maxGlobalElemID), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    }

    const globalIndex mpiElemLimit = geometryMessageSize;

    bufRegionGElem.resize(fd.numElementRegions);

    for (localIndex i = 0; i < fd.numElementRegions; i++)
    {
      bufRegionGElem[i].resize(mpiElemLimit + 1);
    }

    gArray1d bufMaxNumElementsInRegion(fd.numElementRegions);
    if (partition.m_rank == 0)
    {
      std::map<std::string,globalIndex>::iterator it = fd.maxNumElementsInRegion.begin();
      localIndex i = 0;

      for (; (i < fd.numElementRegions)||(it!=fd.maxNumElementsInRegion.end()) ; i++,++it)
      {
        bufMaxNumElementsInRegion[i] = it->second;
      }
    }

    //Broadcast
    {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(bufMaxNumElementsInRegion.data(), fd.numElementRegions, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    }

    if (partition.m_rank != 0)
    {
      std::map<std::string,globalIndex>::iterator it = fd.maxNumElementsInRegion.begin();
      localIndex i = 0;

      for (; i < (fd.numElementRegions)||(it!=fd.maxNumElementsInRegion.end()) ; i++,++it)
      {
        it->second = bufMaxNumElementsInRegion[i];
      }
    }

    gArray1d bufElemToNodes( (mpiElemLimit+1)*maxNumNodesPerElem );//prepare bufElemToNodes here and prepare actual fd.ElemToNodes later on


    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*----------------------MPI-ELEMENT-DATA-INTENSIVE-STARTED-------------------------*/
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    globalIndex currentRegion = 0;//record which element region we are working on currently

    if (partition.m_rank == 0) fd.GoToSection(ELEM);

    globalIndex sizeofBuf;
    globalIndex numberOfRegions = fd.numElementRegions;
    gArray1d numElementsInRegion;
    numElementsInRegion.resize(numberOfRegions);
    numElementsInRegion = 0;

    int completedFlag = 0;
    while ( fd.maxGlobalElemID > 0 && completedFlag==0 )
    {
      tempRegionGElemMap.clear();
      tempElemToNodesMap.clear();
      int modeChange;

      sizeofBuf = 0;

      //Only master processor reads the data
      if (partition.m_rank == 0)
      {
        while (fd.OK())
        {

          if( currentRegion > fd.numElementRegions && fd.mode == UNDEF) break;
          modeChange = fd.AdvanceLine(inputline);
          if (modeChange)
          {
            if(fd.mode == ELEM)
            {
              currentRegion++;
              regionName = FileManagerDataT::ExtractValue("ELSET", inputline);
              ii = 0;
            }
            else
              fd.mode = UNDEF;
          }
          else if (fd.mode == ELEM)
          {
            //retrieve global element index
            globalIndex globalElemNum;
            std::istringstream linestream(inputline);
            linestream >> globalElemNum;

            //add element to the region map with size mpiElemLimit
            tempRegionGElemMap[regionName].insert(globalElemNum);

            for (globalIndex a = 0; a < numActualNodesPerElemMap[regionName]; ++a)
            {
              globalIndex globalNodeNum;
              linestream >> globalNodeNum;
              tempElemToNodesMap[globalElemNum].push_back(globalNodeNum);
            }

            //Special treatment for prisms.
            if (numActualNodesPerElemMap[regionName] == 6 && numNodesPerElemMap[regionName] == 8)
            {
              gArray1d hexNodeMapForPrism;
              hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][4]);
              hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][1]);
              hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][0]);
              hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][3]);
              hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][5]);
              hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][2]);
              hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][5]);
              hexNodeMapForPrism.push_back(tempElemToNodesMap[globalElemNum][2]);
              tempElemToNodesMap[globalElemNum].clear();
              tempElemToNodesMap[globalElemNum] = hexNodeMapForPrism;
            }



            sizeofBuf++;
            ii++;
            numElementsInRegion[currentRegion-1] = ii;
            if(sizeofBuf >= mpiElemLimit)
            {
              break;
            }
          }
        }
        if( !(fd.OK()) )
        {
          completedFlag = 1;
        }
      }

      {

        MPI_Bcast(&completedFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Bcast(&sizeofBuf, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&currentRegion, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numberOfRegions, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        // EBH: the size of numElementsInRegion will be the same as the number incremented in currentRegion
        // because it is incremented each time a new element region is found in the file.  Processors other
        // than the root will need to resize here to accept the proper amount.
        if (partition.m_rank != 0)
        {
          numElementsInRegion.resize(numberOfRegions);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(numElementsInRegion.data(), numberOfRegions, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      }





      if (partition.m_rank == 0)
      {
        unsigned int i=0;
        for( std::map<std::string,gSet>::iterator iterGElemMap=tempRegionGElemMap.begin() ; iterGElemMap!=tempRegionGElemMap.end() ; ++iterGElemMap, ++i )
        {
          const gSet& gElemMap = iterGElemMap->second;
          unsigned int j=0;
          for( gSet::iterator iterElemMap=gElemMap.begin() ; iterElemMap!=gElemMap.end() ; ++iterElemMap, ++j )
          {
            bufRegionGElem[i][j] = *iterElemMap;
          }
        }
      }

      // EBH: if the cpu is not the master, then it must allocate space in tempRegionGElemMap. This is done in a stupid way by just allocating
      // zeros in the spaces
      if (partition.m_rank != 0)
      {
        for (localIndex i = 0; i < fd.numElementRegions; ++i)
        {
          regionName = fd.elementRegionNames(i);
          fd.numElementsInRegion[regionName] = numElementsInRegion[i];
          for (globalIndex j = 0; j < numElementsInRegion[i]; ++j)
          {
            tempRegionGElemMap[regionName].insert(0);
          }
        }
      }

      //Broadcast to all other processors with the same iterator(s)
      {
        unsigned int i=0;
        for( std::map<std::string,gSet>::iterator iterGElemMap=tempRegionGElemMap.begin() ; iterGElemMap!=tempRegionGElemMap.end() ; ++iterGElemMap, ++i )
        {
          MPI_Barrier(MPI_COMM_WORLD);
          MPI_Bcast(bufRegionGElem[i].data(), numElementsInRegion[i], MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        }
      }

      //Pack up
      if (partition.m_rank == 0)
      {
        std::map<globalIndex,gArray1d>::iterator it = tempElemToNodesMap.begin();
        //elemToNodes
        for (globalIndex globalElemNum = 0; ( globalElemNum < sizeofBuf)||(it != tempElemToNodesMap.end() ); globalElemNum++, ++it)
        {
          gArray1d::iterator it2 = it->second.begin();
          for (localIndex a = 0; (a < maxNumNodesPerElem)||(it2 != it->second.end()); ++a, ++it2)
          {
            bufElemToNodes[globalElemNum*maxNumNodesPerElem+a] = *it2;
          }
        }
      }

      //Broadcast
      {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(bufElemToNodes.data(), (sizeofBuf)*maxNumNodesPerElem, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      }

      //Copy buffer content into memory
      if (partition.m_rank != 0)
      {
        unsigned int i=0, elemIdx=0;
        for( std::map<std::string,gSet>::iterator iterGElemMap=tempRegionGElemMap.begin() ; iterGElemMap!=tempRegionGElemMap.end() ; ++iterGElemMap, ++i )
        {
          //tempElemToNodesMap
          for (globalIndex globalElemNum = 0; globalElemNum < numElementsInRegion[i]; globalElemNum++)
          {
            for (localIndex a = 0; a < maxNumNodesPerElem; ++a)
            {
              if(bufElemToNodes[elemIdx*maxNumNodesPerElem+a] != 0)
                tempElemToNodesMap[bufRegionGElem[i][globalElemNum]].push_back(bufElemToNodes[elemIdx*maxNumNodesPerElem+a]);
            }
            ++elemIdx;
          }
        }
      }


      globalIndex globalElemNum;
      unsigned int i=0;
      //for( std::map<std::string,gSet>::iterator iterGElemMap=tempRegionGElemMap.begin() ; iterGElemMap!=tempRegionGElemMap.end() ; ++iterGElemMap, ++i )
      for( i = 0 ; i<fd.elementRegionNames.size() ; ++i )
      {
        // reset numElementsInRegion for fd since each processor needs to figure out how many elements it contains
        regionName = fd.elementRegionNames[i];
        fd.numElementsInRegion[regionName] = 0;
        for (globalIndex j = 0; j < numElementsInRegion[i]; j++)
        {
          globalElemNum = bufRegionGElem[i][j];

          //get element center
          R1Tensor elemCenter;//Warning: has to be defined inside so it can be initialized to be zero whenever it is declared here.

          bool findGlobalElem = 1;

          //Find if each node of a global element lies in the potential NodalPositionsMap
          for(gArray1d::iterator it = tempElemToNodesMap[globalElemNum].begin(); it != tempElemToNodesMap[globalElemNum].end(); ++it)
          {
            if (fd.potentialNodalPositionsMap.find(*it) != fd.potentialNodalPositionsMap.end())
              continue;
            else
            {
              findGlobalElem = 0;
              break;
            }
          }

          //The element will only be tested when all of its nodes lies in the potential NodalPositionsMap
          if(findGlobalElem == 1)
          {
            for(gArray1d::iterator it = tempElemToNodesMap[globalElemNum].begin(); it != tempElemToNodesMap[globalElemNum].end(); ++it)
              elemCenter += fd.potentialNodalPositionsMap[*it];
            elemCenter /= numNodesPerElemMap[regionName];
            //determine whether the element lies within this partition
            if (partition.IsCoordInPartition(elemCenter))
            {
              ++fd.numElementsInRegion[regionName];
              fd.elemsInRegion[regionName].push_back(globalElemNum);

              for(gArray1d::iterator it = tempElemToNodesMap[globalElemNum].begin(); it != tempElemToNodesMap[globalElemNum].end(); ++it)
              {
                fd.nodalPositionsMap[*it] = fd.potentialNodalPositionsMap[*it];
                fd.elemToNodesMap[globalElemNum].push_back(*it);
              }
            }
          }
        }
      }
//      if(globalElemNum == fd.maxGlobalElemID) break;
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*----------------------MPI-ELEM-DATA-INTENSIVE-FINISISHED-------------------------*/
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


    //set the number of nodes now ...
    SetNumNodes(fd);

  }
  /**
   * @brief Read the finite element Abaqus mesh
   * @author Scott Johnson
   * @param[in,out] domain Physical domain object
   * @param[in,out] partition Spatial partition object
   * @param[in,out] fd File data manager object
   */
  void FileManagerT::ReadFiniteElementAbaqusMesh(PhysicalDomainT& domain, SpatialPartition& partition, AbaqusFileManagerDataT& fd)
  {
    //do the first file sweep
    ReadAbaqusMeshB(partition, fd);

    std::string inputline;
    fd.geometry.clear();

    // allocate node objects
    domain.m_feNodeManager.resize(fd.numNodes*1.5);
    domain.m_feNodeManager.resize(fd.numNodes);

    globalIndex globalNodeNumber = 0;
    localIndex localNodeNumber = 0;
    {
      std::map<globalIndex,R1Tensor>::iterator npos = fd.nodalPositionsMap.begin();

      gArray1d& ltog = domain.m_feNodeManager.m_localToGlobalMap;
      array<R1Tensor>& rpos = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition>();
      for (; npos != fd.nodalPositionsMap.end(); ++npos)
      {

        globalNodeNumber = npos->first;
        ltog[localNodeNumber] = globalNodeNumber;
        rpos[localNodeNumber] = npos->second;

        fd.GlobalToLocalNodeMap[globalNodeNumber] = localNodeNumber;
        ++localNodeNumber;

      }
    }

    std::cout << "rank " << partition.m_rank << " has " << fd.numNodes << " of nodes." << std::endl;
    std::cout << "rank " << partition.m_rank << " has " << fd.numElementRegions
        << " of element regions." << std::endl;
    lvector numElements(fd.numElementRegions);
    lvector numElements2(fd.numElementRegions);
    for (localIndex i = 0; i < fd.numElementRegions; ++i)
    {

      numElements[i] = fd.numElementsInRegion[fd.elementRegionNames[i]];
      std::cout << "  element region " << fd.elementRegionNames[i] << " is of type \""
          << fd.elementRegionTypes[i] << "\", and has " << numElements[i] << " number of elements."
          << std::endl;
      numElements2[i] = numElements[i] * 1.25;
    }

    // allocate element objects
    domain.m_feElementManager.resize(numElements2, fd.elementRegionNames, fd.elementRegionTypes );
    domain.m_feElementManager.resize(numElements, fd.elementRegionNames, fd.elementRegionTypes );

    // now fill the elementfToNode array
//    for (std::map<ElementManagerT::RegKeyType, ElementRegionT>::iterator i =
//        domain.m_feElementManager.m_ElementRegions.begin(); i
//        != domain.m_feElementManager.m_ElementRegions.end(); ++i, ++count)
    for( localIndex count=0 ; count<fd.elementRegionNames.size() ; ++count )
    {
      const std::string& regionName = fd.elementRegionNames[count];
      ElementRegionT& elemRegion = domain.m_feElementManager.m_ElementRegions[regionName];
      elemRegion.m_maxGlobalNumber = fd.maxNumElementsInRegion[regionName];

      const localIndex numNodes = ElementTypeToNumberOfNodes(fd.elementRegionTypes[count]);
      elemRegion.m_numNodesPerElem = numNodes;
      if (elemRegion.m_numElems > 0)
      {
        const array<integer> nodeOrdering = AbaqusNodeOrdering(elemRegion.m_elementGeometryID);
        gArray1d& elemLocalToGlobal = elemRegion.m_localToGlobalMap;

        // fill element node array
        localIndex localElemIndexInRegion = 0;


        for (gArray1d::iterator elemNumber = fd.elemsInRegion[regionName].begin();
                                elemNumber != fd.elemsInRegion[regionName].end(); ++elemNumber)
        {
          elemLocalToGlobal(localElemIndexInRegion) = *elemNumber;
          localIndex* elemToNodeMap = elemRegion.ElementToNodeMap(localElemIndexInRegion++);
          for(localIndex a=0 ; a<elemRegion.m_numNodesPerElem ; ++a )
          {
            elemToNodeMap[a] = fd.GlobalToLocalNodeMap[ fd.elemToNodesMap[*elemNumber][nodeOrdering[a]] ];
          }

        }
      }
    }

    std::string setName;
    std::map<std::string,gSet> mymap;

//    globalIndex *bufSetGnodeNumLength = NULL;
//    globalIndex **bufSetGnode = NULL;

    int numSet;

    // now read the nodeset data
    if (partition.m_rank == 0)
    {
      fd.GoToSection(NSET);
      while (fd.OK())
      {
        //Pack up set name
        const bool modeChange = fd.AdvanceLine(inputline);
        if (modeChange && fd.mode == NSET)
        {
          setName = FileManagerDataT::ExtractValue("NSET", inputline);
          domain.m_feNodeManager.m_Sets[setName];
          mymap[setName];
        }
        else
        {
          if (modeChange)
          {
            fd.mode = UNDEF;
          }
          else if (fd.mode == NSET)
          {
            std::istringstream linestream(inputline);
            localIndex lnodeNum;
            globalIndex gnodeNum = GLOBALINDEX_MAX;
            while (!linestream.eof())
            {
              linestream >> gnodeNum;
              if (gnodeNum != GLOBALINDEX_MAX)
              {
                mymap[setName].insert(gnodeNum);//pack up data sets
                lnodeNum = fd.GlobalToLocalNodeMap[gnodeNum];
                if (lnodeNum != GLOBALINDEX_MAX)
                  domain.m_feNodeManager.m_Sets[setName].insert(lnodeNum);
              }
            }
          }
        }
      }
    }

    if (partition.m_rank == 0)
      numSet = mymap.size();

    //Broadcast
    {
      //numSet
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&numSet, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    {
      bufvector bufferSets;
      bufferSets.Pack(mymap);

      int bufferLength = bufferSets.size();
      MPI_Bcast(&bufferLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
      bufferSets.resize( bufferLength);
      MPI_Bcast(bufferSets.data(), bufferSets.size(), MPI_CHAR, 0, MPI_COMM_WORLD);

      const char* pbuffer = bufferSets.data();
      bufvector::Unpack( pbuffer, mymap );
    }

    //Update
    if (partition.m_rank != 0)
    {
      for ( std::map<std::string,gSet>::const_iterator i = mymap.begin() ; i != mymap.end() ; ++i )
      {
        setName = i->first;
        const gSet& globalSet = i->second;
        domain.m_feNodeManager.m_Sets[setName];
        for ( gSet::const_iterator gNode=globalSet.begin() ; gNode!=globalSet.end() ; ++gNode)
        {
          if ( *gNode != GLOBALINDEX_MAX)
          {
            localIndex lnodeNum = fd.GlobalToLocalNodeMap[*gNode];
            if (lnodeNum != GLOBALINDEX_MAX)
              domain.m_feNodeManager.m_Sets[setName].insert(lnodeNum);
          }
        }
      }
    }

    // close the geometry file
    fd.CloseFile();

    MPI_Barrier(MPI_COMM_WORLD);
    if (partition.m_rank == 0)
      std::cout << "Done Reading Finite Element Mesh.\n";


    //  std::cout<<"Initializing Domain...\n";
    //domain.Initialize();
    //  std::cout<<"Done Initializing Domain.\n";
  }

#ifdef SRC_EXTERNAL

  /**
   * @brief Read the fault patch element Abaqus mesh
   * @author Scott Johnson
   * @param[in,out] domain Physical domain object
   * @param[in,out] partition Spatial partition object
   * @param[in,out] fd File data manager object
   */
  void FileManagerT::ReadFaultPatchElementAbaqusMesh(PhysicalDomainT& domain, SpatialPartition& partition, AbaqusFileManagerDataT& fd)
  {
    //do the first sweep through the file
    ReadAbaqusMeshB(partition, fd);

    //ALLOCATE AND ASSIGN NODAL POSITIONS
    std::map<globalIndex, localIndex> globalToLocal;
    {
      domain.m_faultPatchNodes.resize(fd.numNodes);


      std::map<globalIndex,R1Tensor>::iterator npos = fd.nodalPositionsMap.begin();

      array<R1Tensor>& nref = domain.m_faultPatchNodes.GetFieldData<FieldInfo::referencePosition>();

      globalIndex globalNodeNumber = 0;
      localIndex localNodeNumber = 0;
      for(; npos != fd.nodalPositionsMap.end(); ++npos)
      {
        globalNodeNumber = npos->first;
        nref[localNodeNumber] = npos->second;
        fd.GlobalToLocalNodeMap[globalNodeNumber] = localNodeNumber;
        ++localNodeNumber;
      }
    }

    // allocate element objects
    {
      localIndex numLocalFaces = 0;
      for (std::map<std::string, localIndex>::const_iterator it = fd.numElementsInRegion.begin();
          it != fd.numElementsInRegion.end(); ++it)
        numLocalFaces += it->second;
      domain.m_faultPatchFaces.resize(numLocalFaces);
    }

    // now fill the element array
    {
      domain.m_faultPatchFaces.m_maxGlobalNumber = 0;
      localIndex localElementIndex = 0;
      for (std::map<std::string, gArray1d>::const_iterator it = fd.elemsInRegion.begin();
          it != fd.elemsInRegion.end(); ++it)
      {
        for(gArray1d::const_iterator ei = it->second.begin(); ei != it->second.end(); ++ei)
        {
          const globalIndex globalElementNumber = *ei;
          domain.m_faultPatchFaces.m_maxGlobalNumber = globalElementNumber > domain.m_faultPatchFaces.m_maxGlobalNumber ?
              globalElementNumber : domain.m_faultPatchFaces.m_maxGlobalNumber;

          domain.m_faultPatchFaces.m_toNodesRelation[localElementIndex].resize(4);
          for(localIndex nn = 0; nn < 4; nn++)
          {
            domain.m_faultPatchFaces.m_toNodesRelation[localElementIndex][nn] = fd.GlobalToLocalNodeMap[ fd.elemToNodesMap[globalElementNumber][nn] ];
            //std::cout << " ... read int toNodesRelation[" << localElementIndex << "][" << nn << "] -> " << domain.m_faultPatchFaces.m_toNodesRelation[localElementIndex][nn] << std::endl;
          }
          domain.m_faultPatchFaces.m_localToGlobalMap[localElementIndex] = globalElementNumber;
          domain.m_faultPatchFaces.m_globalToLocalMap[globalElementNumber] = localElementIndex;
          ++localElementIndex;

        }
      }
    }

    // close the geometry file
    fd.CloseFile();
    if (partition.m_rank == 0)
      std::cout << "Done Reading Fault Patch Element Mesh.\n";
  }
#endif

  void FileManagerT::SetTemp(const int filenum)
  {
    char junk[10];
    sprintf(junk, "_%03d", filenum);
    temp = junk;
  }
}

