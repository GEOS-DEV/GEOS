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
 * MaterialBase.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */


#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "MaterialBase.h"
#include "../../IO/ticpp/HierarchicalDataNode.h.old"

MaterialBase::MaterialBase( const int paramSize, const int stateSize ):
  ConstitutiveBase(),
  m_paramSize( paramSize ),
  m_stateSize( stateSize )
{}

MaterialBase::~MaterialBase()
{}

void
MaterialBaseStateData::TotalStress(R2SymTensor& totalStress) const
{
  totalStress = devStress;
  totalStress.PlusIdentity(pressure);
}

void
MaterialBaseStateData::RotateState( const R2Tensor& Rot )
{
  R2SymTensor temp;

  devStress.PlusIdentity( pressure );
  temp.QijAjkQlk(devStress,Rot);
  pressure = temp.Trace() / 3.0;
  devStress = temp;
  devStress.PlusIdentity(-pressure);
}
void
MaterialBase::StrainDrivenUpdateMember( const localIndex index0,
                                        const localIndex index1,
                                        const R2SymTensorT < 3 >& Ddt,
                                        const R2TensorT < 3 >& L,
                                        const R2Tensor& Rot,
                                        const realT dt )
{
  throw GPException("Cannot call MaterialBase::StrainDrivenUpdateMember; must have derived method\n");
}

void
MaterialBase::StrainDrivenUpdateMember( const localIndex index0,
                                        const localIndex index1,
                                        const R2SymTensorT < 3 >& Ddt,
                                        const R2TensorT < 3 >& L,
                                        const R2Tensor& Rot,
                                        const realT& volume_n,
                                        const realT& volume_np1,
                                        const realT dt)
{
  throw GPException("Cannot call MaterialBase::StrainDrivenUpdateMember; must have derived method\n");
}

void
MaterialBase::MeanPressureDevStress( const localIndex index,
                                     realT& pressure, R2SymTensor& devStress) const
{
  throw GPException("Cannot call MaterialBase::MeanPressureDevStress; must have derived method\n");
}
