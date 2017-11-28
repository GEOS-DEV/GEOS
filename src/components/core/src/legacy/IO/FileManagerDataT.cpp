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
 * File: FileManagerDataT.cpp
 * Class provides file IO data structure
 * created : SJ (11/16/2012)
 * ported from work by RRS (10/2001)
 */

#include "FileManagerDataT.h"

namespace GPAC_IO
{
/**
 * @brief Constructor
 * @author Scott Johnson
 */
FileManagerDataT::FileManagerDataT():
  exist(false),
  numNodes(0),
  maxGlobalNodeID(0),
  maxGlobalElemID(0),
  totalNumElem(0),
  numElementRegions(0),
  numNodeSets(0),
  numElementsInRegion(),
  elementRegionNameTypes(),
  elementRegionNames(),
  elementRegionTypes(),
  elemToNodes(),
  nodalPositions(),
  spatialMin(),
  spatialMax(),
  isElemInDomain(),
  isNodeInDomain(),
  elemsInRegion(),
  beginElems(0),
  beginNodes(0),
  beginNodesets(0),
  current(0),
  GlobalToLocalNodeMap(),
  geometry(),
  filename()
{
  spatialMin = 1e100;
  spatialMax = -1e100;
}

/**
 * @brief Brief
 * @author Scott Johnson
 * @param[in] geometryfilename File name for the Abaqus geometry file
 * @return Flag for successfully opening the file
 */
bool FileManagerDataT::OpenFile()
{
  if (filename.empty())
  {
    this->exist = false;
  }
  else
  {
    geometry.open(filename.c_str());
    this->exist = geometry.is_open();
    if(!this->exist)
      throw GPException("Could not open mesh file!");
  }
  return this->exist;
}

/**
 * @brief Write the current state to std out
 * @author Scott Johnson
 */
void FileManagerDataT::Write()
{
  std::cout << "elementRegionNameTypes " << elementRegionNameTypes.size() << "\n";
  for (std::map<std::string, std::string>::const_iterator it = elementRegionNameTypes.begin() ;
       it != elementRegionNameTypes.end() ; ++it)
    std::cout << " " << it->first << " " << it->second << "\n";
  std::cout << "nodalPositions " << nodalPositions.size() << "\n";
  for (array<R1Tensor>::size_type i = 0 ; i < nodalPositions.size() ; ++i)
    std::cout << " " << nodalPositions[i](0) << " " << nodalPositions[i](1) << " "
              << nodalPositions[i](2) << "\n";
  std::cout << "maxGlobalNodeID " << maxGlobalNodeID << "\n";
  std::cout << "numElementsInRegion " << numElementsInRegion.size() << "\n";
  for (std::map<std::string, localIndex>::const_iterator it = numElementsInRegion.begin() ;
       it != numElementsInRegion.end() ; ++it)
    std::cout << " " << it->first << " " << it->second << "\n";
  std::cout << "elemsInRegion " << elemsInRegion.size() << "\n";
  for (std::map<std::string, gArray1d>::const_iterator it = elemsInRegion.begin() ;
       it != elemsInRegion.end() ; ++it)
    for (gArray1d::size_type j = 0 ; j < it->second.size() ; ++j)
      std::cout << " " << it->first << " " << it->second[j] << "\n";
  std::cout << "maxGlobalElemID " << maxGlobalElemID << "\n";
  std::cout << "isNodeInDomain " << isNodeInDomain.size() << "\n";
  for (array<integer>::size_type i = 0 ; i < isNodeInDomain.size() ; ++i)
    std::cout << " " << isNodeInDomain[i] << "\n";
  std::cout << "isElemInDomain " << isElemInDomain.size() << "\n";
  for (array<integer>::size_type i = 0 ; i < isElemInDomain.size() ; ++i)
    std::cout << " " << isElemInDomain[i] << "\n";
  std::cout << "numNodes " << numNodes << "\n";
  std::cout << "elemToNodes " << elemToNodes.size() << "\n";
  //    for(int i = 0; i < elemToNodes.size(); ++i)
  //      std::cout << " " << elemToNodes[i,j] << "\n";
  std::cout << "numElementRegions " << numElementRegions << "\n";
  std::cout << "elementRegionNames " << elementRegionNames.size() << "\n";
  for (array<string>::size_type i = 0 ; i < elementRegionNames.size() ; ++i)
    std::cout << " " << elementRegionNames[i] << "\n";
  std::cout << "elementRegionTypes " << elementRegionTypes.size() << "\n";
  for (array<string>::size_type i = 0 ; i < elementRegionTypes.size() ; ++i)
    std::cout << " " << elementRegionTypes[i] << "\n";
  std::cout << "elementRegionNameTypes " << elementRegionNameTypes.size() << "\n";
  for (std::map<std::string, std::string>::const_iterator it = elementRegionNameTypes.begin() ;
       it != elementRegionNameTypes.end() ; ++it)
    std::cout << " " << it->first << " " << it->second << "\n";
  std::cout << "elemsInRegion " << elemsInRegion.size() << "\n";
  for (std::map<std::string, gArray1d>::const_iterator it = elemsInRegion.begin() ;
       it != elemsInRegion.end() ; ++it)
    for (gArray1d::size_type i = 0 ; i < it->second.size() ; ++i)
      std::cout << " " << it->first << " " << it->second[i] << "\n";
  std::cout << "GlobalToLocalNodeMap " << GlobalToLocalNodeMap.size() << "\n";
  for (lvector::size_type i = 0 ; i < GlobalToLocalNodeMap.size() ; ++i)
    std::cout << " " << GlobalToLocalNodeMap[i] << "\n";
}

void FileManagerDataT::Reset()
{
  numNodes = 0;   // Number of nodes
  maxGlobalNodeID = 0;
  nodalPositions.resize(1);
  isNodeInDomain.resize(1);

  maxGlobalElemID = 0;
  totalNumElem = 0;
  numElementRegions = 0;
}

bool FileManagerDataT::AdvanceLine(std::string& inputline)
{
  current = geometry.tellg();
  getline(geometry, inputline);
  return false;
}

//for parsing format: "key=value,"
std::string FileManagerDataT::ExtractValue(const std::string& key, const std::string& inputline)
{
  const int pos1 = inputline.find(key) + key.length();
  std::string first = inputline.substr(pos1);

  const int pos2 = first.find("=") + 1;
  std::string second = first.substr(pos2);

  const int pos3 = second.find(" ");
  return pos3 > -1 ? second.substr(0, pos3) : second;
}
}
