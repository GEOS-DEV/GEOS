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
 * EncapsulatedObjectBase.h
 *
 *  Created on: Jan 10, 2013
 *      Author: settgast1
 */

#ifndef ENCAPSULATEDOBJECTBASE_H_
#define ENCAPSULATEDOBJECTBASE_H_

#include "Common/Common.h"


namespace EOB_Helpers
{
template< int N, typename T>
inline void SerializeMember( const localIndex subIndex,
                             const unsigned int stride,
                             const localIndex index,
                             const T data[N],
                             array<array<T>*>& vars )
{
  for( int i=0 ; i<N ; ++i )
  {
    (*(vars[subIndex+i*stride]))[index] = data[i];
  }
}


template< int N, typename T >
inline void DeserializeMember( const localIndex varIndex,
                               const unsigned int stride,
                               const localIndex index,
                               const array<array<T>*>& vars,
                               T data[N] )
{
  for( int i=0 ; i<N ; ++i )
  {
    data[i] = (*(vars[varIndex+i*stride]))[index];
  }
}
}



template< int NINT, int NR0, int NR1, int NR2, int NR2S >
class EncapsulatedObjectBase
{
public:



  EncapsulatedObjectBase()
  {}

  ~EncapsulatedObjectBase()
  {}

  EncapsulatedObjectBase( const EncapsulatedObjectBase<NINT, NR0, NR1, NR2, NR2S>& source )
  {
    memcpy( this, &source, sizeof(EncapsulatedObjectBase< NINT, NR0, NR1, NR2, NR2S >) );
  }

  EncapsulatedObjectBase& operator=( const EncapsulatedObjectBase& source )
  {
    memcpy( this, &source, sizeof(EncapsulatedObjectBase< NINT, NR0, NR1, NR2, NR2S >) );
    return *this;
  }


  void Serialize( const localIndex subIndex,
                  const unsigned int stride,
                  const localIndex index,
                  array<array<integer>*>& intVars,
                  array<array<real64>*>& realVars,
                  array<array<R1Tensor>*>& R1Vars,
                  array<array<R2Tensor>*>& R2Vars,
                  array<array<R2SymTensor>*>& R2SymVars ) const
  {
    EOB_Helpers::SerializeMember<NINT>( subIndex, stride, index, Ints(),          intVars );
    EOB_Helpers::SerializeMember<NR0> ( subIndex, stride, index, Reals(),         realVars );
    EOB_Helpers::SerializeMember<NR1> ( subIndex, stride, index, R1Tensors(),     R1Vars );
    EOB_Helpers::SerializeMember<NR2> ( subIndex, stride, index, R2Tensors(),     R2Vars );
    EOB_Helpers::SerializeMember<NR2S>( subIndex, stride, index, R2STensors(),  R2SymVars );
  }

  void Deserialize( const localIndex subIndex,
                    const unsigned int stride,
                    const localIndex index,
                    const array<array<integer>*>& intVars,
                    const array<array<real64>*>& realVars,
                    const array<array<R1Tensor>*>& R1Vars,
                    const array<array<R2Tensor>*>& R2Vars,
                    const array<array<R2SymTensor>*>& R2SymVars )
  {
    EOB_Helpers::DeserializeMember<NINT>( subIndex, stride, index, intVars, Ints() );
    EOB_Helpers::DeserializeMember<NR0> ( subIndex, stride, index, realVars, Reals() );
    EOB_Helpers::DeserializeMember<NR1> ( subIndex, stride, index, R1Vars, R1Tensors() );
    EOB_Helpers::DeserializeMember<NR2> ( subIndex, stride, index, R2Vars, R2Tensors() );
    EOB_Helpers::DeserializeMember<NR2S>( subIndex, stride, index, R2SymVars, R2STensors() );
  }


  realT m_real[NR0];
  R1Tensor m_r1tensor[NR1];
  R2Tensor m_r2tensor[NR2];
  R2SymTensor m_r2symtensor[NR2S];
  int m_int[NINT];

  inline       int* Ints()                     { return m_int; }
  inline const int* Ints()               const { return m_int; }
  inline       realT* Reals()                  { return m_real; }
  inline const realT* Reals()            const { return m_real; }
  inline       R1Tensor* R1Tensors()           { return m_r1tensor; }
  inline const R1Tensor* R1Tensors()     const { return m_r1tensor; }
  inline       R2Tensor* R2Tensors()           { return m_r2tensor; }
  inline const R2Tensor* R2Tensors()     const { return m_r2tensor; }
  inline       R2SymTensor* R2STensors()       { return m_r2symtensor; }
  inline const R2SymTensor* R2STensors() const { return m_r2symtensor; }


  inline static int NumInt() { return NINT;}
  inline static int NumR0()  { return NR0;}
  inline static int NumR1()  { return NR1;}
  inline static int NumR2()  { return NR2;}
  inline static int NumR2S() { return NR2S;}

private:



};


#endif /* ENCAPSULATEDOBJECTBASE_H_ */
