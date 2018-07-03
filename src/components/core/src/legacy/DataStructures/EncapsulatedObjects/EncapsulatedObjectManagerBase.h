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
 * EncapsulatedObjectManagerBase.h
 *
 *  Created on: Jan 16, 2013
 *      Author: settgast1
 */

#ifndef ENCAPSULATEDOBJECTMANAGERBASE_H_
#define ENCAPSULATEDOBJECTMANAGERBASE_H_

namespace TICPP
{
class HierarchicalDataNode;
}

#include "Common/Common.h"
#include "ArrayT/bufvector.h"

class EncapsulatedObjectManagerBase
{
public:
  EncapsulatedObjectManagerBase() {}
  virtual ~EncapsulatedObjectManagerBase() {}

  virtual void resize( const localIndex num,
                       const bool variableParams = false ) = 0;

  virtual void ReadXML( const TICPP::HierarchicalDataNode& node ) = 0;

  virtual void GetVariableNames( array<string>& intVars,
                                 array<string>& realVars,
                                 array<string>& R1TensorVars,
                                 array<string>& R2TensorVars,
                                 array<string>& R2SymTensorVars ) const = 0;


  virtual void Serialize( array<array<integer>*>& intVars,
                          array<array<real64>*>& realVars,
                          array<array<R1Tensor>*>& R1Vars,
                          array<array<R2Tensor>*>& R2Vars,
                          array<array<R2SymTensor>*>& R2SymVars ) const = 0;

  virtual void Deserialize( const array<array<integer>*>& intVars,
                            const array<array<real64>*>& realVars,
                            const array<array<R1Tensor>*>& R1Vars,
                            const array<array<R2Tensor>*>& R2Vars,
                            const array<array<R2SymTensor>*>& R2SymVars  ) = 0;

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking ) = 0;
  virtual unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking ) = 0;
  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer ) = 0;
  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer ) = 0;

protected:

  template< typename LeafClass >
  void ResizeFromDerived( const localIndex num,
                          const bool variableParams );

  template< typename LeafClass >
  void GetVariableNamesFromDerived( array<string>& intVars,
                                    array<string>& realVars,
                                    array<string>& R1TensorVars,
                                    array<string>& R2TensorVars,
                                    array<string>& R2SymTensorVars ) const;

  template< typename LeafClass >
  void SerializeFromDerived( array<array<integer>*>& intVars,
                             array<array<real64>*>& realVars,
                             array<array<R1Tensor>*>& R1Vars,
                             array<array<R2Tensor>*>& R2Vars,
                             array<array<R2SymTensor>*>& R2SymVars ) const;

  template< typename LeafClass >
  void DeserializeFromDerived( const array<array<integer>*>& intVars,
                               const array<array<real64>*>& realVars,
                               const array<array<R1Tensor>*>& R1Vars,
                               const array<array<R2Tensor>*>& R2Vars,
                               const array<array<R2SymTensor>*>& R2SymVars  );

  template< typename LeafClass, typename T_indices >
  unsigned int PackFromDerived( const T_indices& localIndices,
                                bufvector& buffer,
                                const bool doBufferPacking ) const;

  template< typename LeafClass, typename T_indices >
  unsigned int UnpackFromDerived( const T_indices& localIndices,
                                  const char*& buffer );

private:
  EncapsulatedObjectManagerBase( const EncapsulatedObjectManagerBase& );
  EncapsulatedObjectManagerBase& operator=(const EncapsulatedObjectManagerBase&);
};



template< typename LeafClass >
void EncapsulatedObjectManagerBase::ResizeFromDerived( const localIndex num,
                                                       const bool variableParams )
{
  LeafClass& dthis = static_cast<LeafClass&>(*this);

  dthis.m_stateData.resize( num );

  if( variableParams )
    dthis.m_parameterData.resize( num );

}

template< typename LeafClass >
void EncapsulatedObjectManagerBase::GetVariableNamesFromDerived( array<string>& intVars,
                                                                 array<string>& realVars,
                                                                 array<string>& R1TensorVars,
                                                                 array<string>& R2TensorVars,
                                                                 array<string>& R2SymTensorVars ) const
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);
  if( dthis.m_parameterData.size() > 1 )
  {
    LeafClass::ParameterClass::SetVariableNames( intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars);
  }

  LeafClass::StateClass::SetVariableNames( intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars);

}

template< typename LeafClass >
void EncapsulatedObjectManagerBase::SerializeFromDerived( array<array<integer>*>& intVars,
                                                          array<array<real64>*>& realVars,
                                                          array<array<R1Tensor>*>& R1Vars,
                                                          array<array<R2Tensor>*>& R2Vars,
                                                          array<array<R2SymTensor>*>& R2SymVars ) const
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);

  if( dthis.m_parameterData.size() > 1 )
  {
    typename LeafClass::ParameterArrayType::const_iterator iparam = dthis.m_parameterData.begin();
    for( localIndex a=0 ; a<dthis.NumParameterIndex0() ; ++a )
    {
      for( localIndex b=0 ; b<dthis.NumParameterIndex1() ; ++b )
      {
        iparam->Serialize( b, dthis.NumParameterIndex0(), a, intVars, realVars, R1Vars, R2Vars, R2SymVars);
        ++iparam;
      }
    }
  }

  {
    typename LeafClass::StateArrayType::const_iterator istate = dthis.m_stateData.begin();
    for( localIndex a=0 ; a<dthis.NumStateIndex0() ; ++a )
    {
      for( localIndex b=0 ; b<dthis.NumStateIndex1() ; ++b )
      {
        istate->Serialize( b, dthis.NumStateIndex1(), a, intVars, realVars, R1Vars, R2Vars, R2SymVars);
        ++istate;
      }
    }
  }
}

template< typename LeafClass >
void EncapsulatedObjectManagerBase::DeserializeFromDerived( const array<array<integer>*>& intVars,
                                                            const array<array<real64>*>& realVars,
                                                            const array<array<R1Tensor>*>& R1Vars,
                                                            const array<array<R2Tensor>*>& R2Vars,
                                                            const array<array<R2SymTensor>*>& R2SymVars  )
{
  LeafClass& dthis = static_cast<LeafClass&>(*this);

  if( dthis.m_parameterData.size() > 1 )
  {
    typename LeafClass::ParameterArrayType::iterator iparam = dthis.m_parameterData.begin();
    for( localIndex a=0 ; a<dthis.NumParameterIndex0() ; ++a )
    {
      for( localIndex b=0 ; b<dthis.NumParameterIndex1() ; ++b )
      {
        iparam->Deserialize( b, dthis.NumParameterIndex1(), a, intVars, realVars, R1Vars, R2Vars, R2SymVars);
        ++iparam;
      }
    }
  }

  {
    typename LeafClass::StateArrayType::iterator istate = dthis.m_stateData.begin();
    for( localIndex a=0 ; a<dthis.NumStateIndex0() ; ++a )
    {
      for( localIndex b=0 ; b<dthis.NumStateIndex1() ; ++b )
      {
        istate->Deserialize( b, dthis.NumStateIndex1(), a, intVars, realVars, R1Vars, R2Vars, R2SymVars);
        ++istate;
      }
    }
  }

}



template< typename LeafClass, typename T_indices >
unsigned int EncapsulatedObjectManagerBase::PackFromDerived( const T_indices& localIndices, bufvector& buffer, const bool doBufferPacking ) const
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);

  unsigned int packedSize = 0;

  for( typename T_indices::const_iterator index=localIndices.begin() ; index!=localIndices.end() ; ++index )
  {
//    packedSize += buffer.Pack(dthis.m_stateData[*index]);
  }

  if( dthis.m_parameterData.size() > 1 )
  {
    for( typename T_indices::const_iterator index=localIndices.begin() ; index!=localIndices.end() ; ++index )
    {
//      packedSize += buffer.Pack(dthis.m_parameterData[*index]);
    }
  }

  return packedSize;
}

template< typename LeafClass, typename T_indices >
unsigned int EncapsulatedObjectManagerBase::UnpackFromDerived( const T_indices& localIndices, const char*& buffer )
{
  LeafClass& dthis = static_cast<LeafClass&>(*this);

  unsigned int unpackedSize = 0;
  for( typename T_indices::const_iterator index=localIndices.begin() ; index!=localIndices.end() ; ++index )
  {
//    unpackedSize += bufvector::Unpack( buffer, dthis.m_stateData[*index] );
  }

  if( dthis.m_parameterData.size() > 1 )
  {
    for( typename T_indices::const_iterator index=localIndices.begin() ; index!=localIndices.end() ; ++index )
    {
//      unpackedSize += bufvector::Unpack( buffer, dthis.m_parameterData[*index]
// );
    }
  }
  return unpackedSize;
}


#endif /* ENCAPSULATEDOBJECTMANAGERBASE_H_ */
