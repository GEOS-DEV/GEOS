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
 * File: AbaqusFileManagerDataT.cpp
 * Class provides file IO data structure
 * created : SJ (11/16/2012)
 * ported from work by RRS (10/2001)
 */

#include "AbaqusFileManagerDataT.h"

namespace GPAC_IO
{
  const std::string AbaqusFileManagerDataT::nodeKey = "*NODE";
  const std::string AbaqusFileManagerDataT::elemKey = "*ELEMENT";
  const std::string AbaqusFileManagerDataT::nodeSetKey = "*NSET";
  
  const std::string AbaqusFileManagerDataT::typeKey = "TYPE";
  const std::string AbaqusFileManagerDataT::elSetKey = "ELSET";
  
  /**
   * @brief Constructor
   * @author Scott Johnson
   */
  AbaqusFileManagerDataT::AbaqusFileManagerDataT() :
      FileManagerDataT(),
      mode(UNDEF)
  {
  }

  bool AbaqusFileManagerDataT::AdvanceLine(std::string& inputline)
  {
    FileManagerDataT::AdvanceLine(inputline);

    replace(inputline.begin(), inputline.end(), ',', ' ');

    if (inputline.compare(0, nodeKey.size(), nodeKey) == 0)
    {
      if (std::streamoff(beginNodes) == 0)
        beginNodes = this->current;
      mode = NODE;
      return true;
    }
    else if (inputline.compare(0, elemKey.size(), elemKey) == 0)
    {
      if (std::streamoff(beginElems) == 0)
        beginElems = this->current;
      mode = ELEM;
      return true;
    }
    else if (inputline.compare(0, nodeSetKey.size(), nodeSetKey) == 0)
    {
      if (std::streamoff(beginNodesets) == 0)
        beginNodesets = this->current;
      mode = NSET;
      return true;
    }
    else if (inputline.compare(0, 1, "*") == 0 || inputline.size() == 0)
    {
      bool rval = false;
      if( mode != UNDEF )
      {
        rval = true;
      }
      mode = UNDEF;
      return rval;
    }
    return false;
  }

  //R.W. get maxGlobalNodeID, nodal position
  void AbaqusFileManagerDataT::AddNodeLine(std::string& inputline, const realT geometryUnits, std::map<globalIndex,R1Tensor> &tempNodalPositionsMap, std::map<globalIndex,globalIndex> &resortNodeMap)
  {
    std::istringstream linestream(inputline);
    globalIndex globalNodeNumber, newglobalNodeNumber;
    R1Tensor nodePosition;
    {
      linestream >> globalNodeNumber >> nodePosition(0) >> nodePosition(1) >> nodePosition(2);
      newglobalNodeNumber = resortNodeMap.size();
      resortNodeMap[globalNodeNumber] = newglobalNodeNumber;
      //globalNodeNumber = newglobalNodeNumber + 1;
      nodePosition*=geometryUnits; // rescale to problem units
    }
    if (newglobalNodeNumber > maxGlobalNodeID)
    {
      maxGlobalNodeID = newglobalNodeNumber;
    }

    tempNodalPositionsMap[newglobalNodeNumber] = nodePosition;

  }

  //without MPI: just to get the basic spatial partition information
  void AbaqusFileManagerDataT::AddNodalPositionLine(std::string& inputline, const realT geometryUnits)
  {
    std::istringstream linestream(inputline);
    globalIndex globalNodeNumber;
    R1Tensor nodePosition;
    {
      linestream >> globalNodeNumber >> nodePosition(0) >> nodePosition(1) >> nodePosition(2);
      nodePosition*=geometryUnits; // rescale to problem units
    }
    if (globalNodeNumber > maxGlobalNodeID)
    {
      maxGlobalNodeID = globalNodeNumber;
    }

    spatialMin.SetMin(nodePosition);
    spatialMax.SetMax(nodePosition);
  }

  void AbaqusFileManagerDataT::AddNodalPositionLineDiscrete(std::string& inputline, const realT geometryUnits)
  {
    std::istringstream linestream(inputline);
    globalIndex globalNodeNumber;
    R1Tensor nodePosition;
    {
      linestream >> globalNodeNumber >> nodePosition(0) >> nodePosition(1) >> nodePosition(2);
      nodePosition*=geometryUnits; // rescale to problem units
    }
    if (globalNodeNumber > maxGlobalNodeID)
    {
      maxGlobalNodeID = globalNodeNumber;
      nodalPositions.resize(maxGlobalNodeID + 1);
      isNodeInDomain.resize(maxGlobalNodeID + 1);
    }

    isNodeInDomain[globalNodeNumber] = 0;
    nodalPositions[globalNodeNumber] = nodePosition;

    spatialMin.SetMin(nodePosition);
    spatialMax.SetMax(nodePosition);
  }

  //with MPI: read certain number of nodes into nodalPositions with the size <= mpiNodeLimit
  void AbaqusFileManagerDataT::AddNodalPositionLine(std::string& inputline, const realT geometryUnits, std::map<globalIndex,R1Tensor> &tempNodalPositionsMap)
  {
    std::istringstream linestream(inputline);
    globalIndex globalNodeNumber;
    R1Tensor nodePosition;
    {
      linestream >> globalNodeNumber >> nodePosition(0) >> nodePosition(1) >> nodePosition(2);
      nodePosition*=geometryUnits; // rescale to problem units
    }

    tempNodalPositionsMap[globalNodeNumber] = nodePosition;
}

  void AbaqusFileManagerDataT::AddElementRegionLine(std::string& inputline)
  {
    ++numElementRegions;
    elementRegionTypes.push_back(FileManagerDataT::ExtractValue(typeKey, inputline));
    elementRegionNames.push_back(FileManagerDataT::ExtractValue(elSetKey, inputline));
    numElementsInRegion[elementRegionNames.back()];
    elemsInRegion[elementRegionNames.back()];
    elementRegionNameTypes[elementRegionNames.back()] = elementRegionTypes.back();
  }

  void AbaqusFileManagerDataT::Reset()
  {
    FileManagerDataT::Reset();

    mode = UNDEF;
  }

  void AbaqusFileManagerDataT::GoToSection(READMODE mode_)
  {
    geometry.clear();
    switch(mode_)
    {
    case NODE:
      geometry.seekg(beginNodes);
      break;
    case ELEM:
      geometry.seekg(beginElems);
      break;
    case NSET:
      geometry.seekg(beginNodesets);
      break;
    default:
      break;
    }
  }
}
