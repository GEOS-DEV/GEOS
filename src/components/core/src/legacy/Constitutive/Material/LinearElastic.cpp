// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.


/*
 * LinearElastic.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */


#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "LinearElastic.h"
#include "Constitutive/Material/MaterialFactory.h"
#include <typeinfo>
#include <assert.h>
#include "../../IO/ticpp/HierarchicalDataNode.h.old"
#include "../../IO/ticpp/HierarchicalDataNode.h.old"

LinearElastic::LinearElastic( ):
  LinearElasticIntermediate( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}

LinearElastic::~LinearElastic()
{}

void
LinearElasticParameterData::PostReadXML( const TICPP::HierarchicalDataNode& node )
{
  realT M = 0;
  FillLinearElasticModuli(init_bulkModulus, init_shearModulus, E, Nu, Lame, M);
}
void LinearElastic::InitializeStates(const localIndex index)
{

  for (localIndex a = 0 ; a < m_stateData.Dimension(1) ; ++a)
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
  //  realT strainEnergyIncrement = - 0.5 * ( devStress.Inner() / (2*G) +
  // pow(pressure,2)/K ) * volume_n;

  pressure += trDdt * K;

  {
    R2SymTensorT<3> temp;
    temp = Ddt;
    temp.PlusIdentity(-trDdt / 3.0);
    temp *= 2.0 * G;
    devStress += temp;
  }

  StressPowerIncrement += 0.5 * (Dot(devStress, Ddt) + pressure * trDdt) * volume_np1;
  //  strainEnergyIncrement += 0.5 * ( devStress.Inner() / (2*G) +
  // pow(pressure,2)/K ) * volume_np1;

  matState.ElasticStrainEnergy += StressPowerIncrement; //* 0.5 * ( volume_n +
                                                        // volume_np1);
  matState.StressPower += StressPowerIncrement; //* 0.5 * ( volume_n +
                                                // volume_np1 );

  matState.RotateState(Rot);
  return;
}

void LinearElastic::MeanPressureDevStress(const localIndex index, realT& pressure,
                                          R2SymTensor& devStress) const
{
  MeanPressureDevStressFromDerived<LinearElastic>(index, pressure, devStress);
}
void
LinearElastic::PreSetValues(const array<string>& names)
{
  //reset mechanical constants
  for (localIndex a = 0 ; a < m_parameterData.size() ; ++a)
  {
    m_parameterData[a].E = 0;
    m_parameterData[a].Lame = 0;
    m_parameterData[a].Nu = 0;
    m_parameterData[a].init_bulkModulus = 0;
    m_parameterData[a].init_shearModulus = 0;
  }
}

void
LinearElastic::PostSetValues(const array<string>& names)
{
  //recalculate mechanical parameters that were not already explicitly set
  for(localIndex a = 0 ; a < m_parameterData.size() ; a++)
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
