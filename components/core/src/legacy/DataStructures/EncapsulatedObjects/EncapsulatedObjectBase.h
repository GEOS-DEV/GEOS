/*
 * EncapsulatedObjectBase.h
 *
 *  Created on: Jan 10, 2013
 *      Author: settgast1
 */

#ifndef ENCAPSULATEDOBJECTBASE_H_
#define ENCAPSULATEDOBJECTBASE_H_

#include "Common/Common.h"

#define EOB_Arrays 2

namespace EOB_Helpers
{
template< int N, typename T>
inline void SerializeMember( const localIndex subIndex,
                             const unsigned int stride,
                             const localIndex index,
                             const T* const data,
                             Array1dT<Array1dT<T>*>& vars )
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
                               const Array1dT<Array1dT<T>*>& vars,
                               T* const data )
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
                  Array1dT<iArray1d*>& intVars,
                  Array1dT<rArray1d*>& realVars,
                  Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                  Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                  Array1dT<Array1dT<R2SymTensor>*>& R2SymVars ) const
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
                    const Array1dT<iArray1d*>& intVars,
                    const Array1dT<rArray1d*>& realVars,
                    const Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                    const Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                    const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars )
  {
    EOB_Helpers::DeserializeMember<NINT>( subIndex, stride, index, intVars, Ints() );
    EOB_Helpers::DeserializeMember<NR0> ( subIndex, stride, index, realVars, Reals() );
    EOB_Helpers::DeserializeMember<NR1> ( subIndex, stride, index, R1Vars, R1Tensors() );
    EOB_Helpers::DeserializeMember<NR2> ( subIndex, stride, index, R2Vars, R2Tensors() );
    EOB_Helpers::DeserializeMember<NR2S>( subIndex, stride, index, R2SymVars, R2STensors() );
  }


#if ( EOB_Arrays == 0 || EOB_Arrays == 1)

#if ( EOB_Arrays == 0 )
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

#else
  std::array<realT,NR0> m_real;
  std::array<R1Tensor,NR1> m_r1tensor;
  std::array<R2Tensor,NR2> m_r2tensor;
  std::array<R2SymTensor,NR2S> m_r2symtensor;
  std::array<int,NINT> m_int;

  inline       int* Ints()                     { return m_int.data(); }
  inline const int* Ints()               const { return m_int.data(); }
  inline       realT* Reals()                  { return m_real.data(); }
  inline const realT* Reals()            const { return m_real.data(); }
  inline       R1Tensor* R1Tensors()           { return m_r1tensor.data(); }
  inline const R1Tensor* R1Tensors()     const { return m_r1tensor.data(); }
  inline       R2Tensor* R2Tensors()           { return m_r2tensor.data(); }
  inline const R2Tensor* R2Tensors()     const { return m_r2tensor.data(); }
  inline       R2SymTensor* R2STensors()       { return m_r2symtensor.data(); }
  inline const R2SymTensor* R2STensors() const { return m_r2symtensor.data(); }

#endif


  inline       int& Ints( const int i )                     { return m_int[i]; }
  inline const int& Ints( const int i )               const { return m_int[i]; }
  inline       realT& Reals( const int i )                  { return m_real[i]; }
  inline const realT& Reals( const int i )            const { return m_real[i]; }
  inline       R1Tensor& R1Tensors( const int i )           { return m_r1tensor[i]; }
  inline const R1Tensor& R1Tensors( const int i )     const { return m_r1tensor[i]; }
  inline       R2Tensor& R2Tensors( const int i )           { return m_r2tensor[i]; }
  inline const R2Tensor& R2Tensors( const int i )     const { return m_r2tensor[i]; }
  inline       R2SymTensor& R2STensors( const int i )       { return m_r2symtensor[i]; }
  inline const R2SymTensor& R2STensors( const int i ) const { return m_r2symtensor[i]; }

#else

  enum
  {
    intsSize = (NINT+1)*sizeof(int)/sizeof(realT),
    realsSize = NR0,
    R1TensorsSize = NR1*sizeof(R1Tensor)/sizeof(realT),
    R2TensorsSize = NR2*sizeof(R2Tensor)/sizeof(realT),
    R2STensorsSize = NR2S*sizeof(R2SymTensor)/sizeof(realT),
    dataSize = intsSize + realsSize + R1TensorsSize + R2TensorsSize + R2STensorsSize,
    int_0 = 0,
    real_0 = intsSize,
    R1_0 = real_0 + realsSize,
    R2_0 = R1_0 + R1TensorsSize,
    R2S_0 = R2_0 + R2TensorsSize
  };


  realT m_data[ dataSize ];
  inline       int* Ints()                     { return reinterpret_cast<int*>(m_data); }
  inline const int* Ints()               const { return reinterpret_cast<const int*>(m_data); }
  inline       realT* Reals()                  { return reinterpret_cast<realT*>(&(m_data[real_0])); }
  inline const realT* Reals()            const { return reinterpret_cast<const realT*>(&(m_data[real_0])); }
  inline       R1Tensor* R1Tensors()           { return reinterpret_cast<R1Tensor*>(&(m_data[R1_0])); }
  inline const R1Tensor* R1Tensors()     const { return reinterpret_cast<const R1Tensor*>(&(m_data[R1_0])); }
  inline       R2Tensor* R2Tensors()           { return reinterpret_cast<R2Tensor*>(&(m_data[R2_0])); }
  inline const R2Tensor* R2Tensors()     const { return reinterpret_cast<const R2Tensor*>(&(m_data[R2_0])); }
  inline       R2SymTensor* R2STensors()       { return reinterpret_cast<R2SymTensor*>(&(m_data[R2S_0])); }
  inline const R2SymTensor* R2STensors() const { return reinterpret_cast<const R2SymTensor*>(&(m_data[R2S_0])); }

  inline       int& Ints( const int i )                     { return Ints()[i]; }
  inline const int& Ints( const int i )               const { return Ints()[i]; }
  inline       realT& Reals( const int i )                  { return Reals()[i]; }
  inline const realT& Reals( const int i )            const { return Reals()[i]; }
  inline       R1Tensor& R1Tensors( const int i )           { return R1Tensors()[i]; }
  inline const R1Tensor& R1Tensors( const int i )     const { return R1Tensors()[i]; }
  inline       R2Tensor& R2Tensors( const int i )           { return R2Tensors()[i]; }
  inline const R2Tensor& R2Tensors( const int i )     const { return R2Tensors()[i]; }
  inline       R2SymTensor& R2STensors( const int i )       { return R2STensors()[i]; }
  inline const R2SymTensor& R2STensors( const int i ) const { return R2STensors()[i]; }

#endif

  inline static int NumInt() { return NINT;}
  inline static int NumR0()  { return NR0;}
  inline static int NumR1()  { return NR1;}
  inline static int NumR2()  { return NR2;}
  inline static int NumR2S() { return NR2S;}

private:



};


#endif /* ENCAPSULATEDOBJECTBASE_H_ */
