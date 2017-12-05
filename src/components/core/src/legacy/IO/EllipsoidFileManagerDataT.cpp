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
 * File: EllipsoidFileManagerDataT.cpp
 * Class provides file IO data structure
 * created : SJ (11/16/2012)
 * ported from work by RRS (10/2001)
 */

#include "EllipsoidFileManagerDataT.h"

namespace GPAC_IO
{
EllipsoidFileManagerDataT::EllipsoidFileManagerDataT(): FileManagerDataT()
{}

realT EllipsoidFileManagerDataT::ReadLine(const std::string& inputline,
                                          globalIndex& globalNodeNumber,
                                          R1Tensor& position)
{
  std::istringstream linestream(inputline);
  R1Tensor principalRadii;
  linestream >> globalNodeNumber >> position(0) >> position(1) >> position(2) >>
  principalRadii(0) >> principalRadii(1) >> principalRadii(2);
  return principalRadii.MaxVal();
}

void EllipsoidFileManagerDataT::ReadLine(const std::string& inputline,
                                         globalIndex& globalNodeNumber,
                                         R1Tensor& position,
                                         R1Tensor& principalRadii,
                                         R1TensorT<4>& rotation)
{
  std::istringstream linestream(inputline);
  linestream >> globalNodeNumber >> position(0) >> position(1) >> position(2) >>
  principalRadii(0) >> principalRadii(1) >> principalRadii(2) >>
  rotation(0) >> rotation(1) >> rotation(2) >> rotation(3);
}

bool EllipsoidFileManagerDataT::AdvanceLine(std::string& inputline)
{
  FileManagerDataT::AdvanceLine(inputline);
  replace(inputline.begin(), inputline.end(), ',', ' ');
  return (inputline.compare(0, 1, "*") == 0 || inputline.size() == 0);
}

void EllipsoidFileManagerDataT::AddNodalPositionLine(std::string& inputline, const realT geometryUnits)
{
  R1Tensor nodePosition;
  globalIndex globalNodeNumber;
  realT rad = 0.0;
  {
    rad = ReadLine(inputline, globalNodeNumber, nodePosition);
    nodePosition *= geometryUnits;   // rescale to problem units
  }
  if (globalNodeNumber > maxGlobalNodeID)
  {
    maxGlobalNodeID = globalNodeNumber;
    nodalPositions.resize(maxGlobalNodeID + 1);
    isNodeInDomain.resize(maxGlobalNodeID + 1);
  }

  isNodeInDomain[globalNodeNumber] = 0;
  nodalPositions[globalNodeNumber] = nodePosition;

  R1Tensor tpos(nodePosition);
  tpos += rad;
  spatialMax.SetMax(tpos);
  tpos -= 2*rad;
  spatialMin.SetMin(tpos);
}
}
