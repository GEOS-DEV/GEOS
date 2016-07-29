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


/*
 * LinearElastic.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 

#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "IO/ticpp/HierarchicalDataNode.h"

#include "LinearElastic.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "Constitutive/Material/MaterialFactory.h"
#include <typeinfo>
#include <assert.h>

LinearElastic::LinearElastic( ):
LinearElasticIntermediate( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}

LinearElastic::~LinearElastic()
{

}

void
LinearElasticParameterData::PostReadXML( const TICPP::HierarchicalDataNode& node )
{
  realT M = 0;
  FillLinearElasticModuli(init_bulkModulus, init_shearModulus, E, Nu, Lame, M);
}
void LinearElastic::InitializeStates(const localIndex index)
{

  for (localIndex a = 0; a < m_stateData.Dimension(1); ++a)
  {
    const localIndex paramIndex = m_parameterData.size() > 1 ? index : 0;

    m_stateData[index][a].density = m_parameterData[paramIndex].init_density;
    m_stateData[index][a].BulkModulus = m_parameterData[paramIndex].init_bulkModulus;
    m_stateData[index][a].ShearModulus = m_parameterData[paramIndex].init_shearModulus;
    m_stateData[index][a].ElasticBulkModulus = m_parameterData[paramIndex].init_bulkModulus;
    m_stateData[index][a].ElasticShearModulus = m_parameterData[paramIndex].init_shearModulus;
  }

}

void
LinearElastic::StrainDrivenUpdateMember(const localIndex index0,
                                        const localIndex index1,
                                        const R2SymTensorT<3>& Ddt,
                                        const R2TensorT < 3 >& L,
                                        const R2Tensor& Rot,
                                        const realT dt)
{
  const localIndex paramIndex = m_parameterData.size() > 1 ? index0 : 0;
  const LinearElasticParameterData& matParams = m_parameterData[paramIndex];
  LinearElasticStateData& matState = m_stateData(index0, index1);

  realT& pressure = matState.pressure;
  R2SymTensorT<3>& devStress = matState.devStress;
  const realT& K = matParams.init_bulkModulus;
  const realT& G = matParams.init_shearModulus;

  const realT trDdt = Ddt.Trace();

  R2SymTensorT<3> temp;

  pressure += trDdt * K;

  temp = Ddt;
  temp.PlusIdentity(-trDdt / 3.0);
  temp *= 2.0 * G;
  devStress += temp;

  matState.RotateState(Rot);
  return;
}

void
LinearElastic::StrainDrivenUpdateMember(const localIndex index0,
                                        const localIndex index1,
                                        const R2SymTensorT<3>& Ddt,
                                        const R2TensorT < 3 >& L,
                                        const R2Tensor& Rot,
                                        const realT& volume_n,
                                        const realT& volume_np1,
                                        const realT dt)
{
  const localIndex paramIndex = m_parameterData.size() > 1 ? index0 : 0;
  const LinearElasticParameterData& matParams = m_parameterData[paramIndex];
  LinearElasticStateData& matState = m_stateData(index0, index1);

  realT& pressure = matState.pressure;
  R2SymTensorT<3>& devStress = matState.devStress;
  const realT& K = matParams.init_bulkModulus;
  const realT& G = matParams.init_shearModulus;

  const realT trDdt = Ddt.Trace();

  realT StressPowerIncrement = 0.5 * (Dot(devStress, Ddt) + pressure * trDdt) * volume_n;
  //  realT strainEnergyIncrement = - 0.5 * ( devStress.Inner() / (2*G) + pow(pressure,2)/K ) * volume_n;

  pressure += trDdt * K;

  {
    R2SymTensorT<3> temp;
    temp = Ddt;
    temp.PlusIdentity(-trDdt / 3.0);
    temp *= 2.0 * G;
    devStress += temp;
  }

  StressPowerIncrement += 0.5 * (Dot(devStress, Ddt) + pressure * trDdt) * volume_np1;
  //  strainEnergyIncrement += 0.5 * ( devStress.Inner() / (2*G) + pow(pressure,2)/K ) * volume_np1;

  matState.ElasticStrainEnergy += StressPowerIncrement; //* 0.5 * ( volume_n + volume_np1);
  matState.StressPower += StressPowerIncrement; //* 0.5 * ( volume_n + volume_np1 );

  matState.RotateState(Rot);
  return;
}

void LinearElastic::MeanPressureDevStress(const localIndex index, realT& pressure,
                                          R2SymTensor& devStress) const
{
  MeanPressureDevStressFromDerived<LinearElastic>(index, pressure, devStress);
}
void
LinearElastic::PreSetValues(const sArray1d& names)
{
  //reset mechanical constants
  for (localIndex a = 0; a < m_parameterData.size(); ++a)
  {
    m_parameterData[a].E = 0;
    m_parameterData[a].Lame = 0;
    m_parameterData[a].Nu = 0;
    m_parameterData[a].init_bulkModulus = 0;
    m_parameterData[a].init_shearModulus = 0;
  }
}

void
LinearElastic::PostSetValues(const sArray1d& names)
{
  //recalculate mechanical parameters that were not already explicitly set
  for(localIndex a = 0; a < m_parameterData.size(); a++)
  {
    realT M = 0;
    FillLinearElasticModuli(m_parameterData[a].init_bulkModulus,
                            m_parameterData[a].init_shearModulus,
                            m_parameterData[a].E,
                            m_parameterData[a].Nu,
                            m_parameterData[a].Lame,
                            M);
  }
}

/// Register class in the class factory
REGISTER_MATERIAL( LinearElastic )
