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
