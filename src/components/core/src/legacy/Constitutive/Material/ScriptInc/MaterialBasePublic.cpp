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

//FUNCTION_BEGIN_PARSE
virtual_void
MaterialBase::StrainDrivenUpdateMember( const localIndex index0,
                                        const localIndex index1,
                                        const R2SymTensorT < 3 >& Ddt,
                                        const R2TensorT < 3 >& L,
                                        const R2Tensor& Rot,
                                        const realT dt )
{
  throw GPException("Cannot call MaterialBase::StrainDrivenUpdateMember; must have derived method\n");
}

//FUNCTION_BEGIN_PARSE
virtual_void
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

//FUNCTION_BEGIN_PARSE
virtual_void
MaterialBase::MeanPressureDevStress( const localIndex index,
                                     realT& pressure, R2SymTensor& devStress) const
{
  throw GPException("Cannot call MaterialBase::MeanPressureDevStress; must have derived method\n");
}

//FUNCTION_BEGIN_PARSE
template< typename LeafClass > void
MaterialBase::MeanPressureDevStressFromDerived( const localIndex index,
                                                realT& pressure, R2SymTensor& devStress) const
{
  pressure = 0;
  devStress = 0;

  const LeafClass& dthis = static_cast<const LeafClass&>(*this);
  const typename LeafClass::StateClass* const state = dthis.m_stateData[index];

  for( localIndex b=0 ; b<dthis.NumStateIndex1() ; ++b )
  {
    devStress += state[b].devStress;
    pressure += state[b].pressure;
  }

  if(dthis.NumStateIndex1() > 0)
  {
    pressure /= dthis.NumStateIndex1();
    devStress /= dthis.NumStateIndex1();
  }
}
