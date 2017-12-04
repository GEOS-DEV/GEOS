/*
 * ConstitutiveBase.h
 *
 *  Created on: Jan 16, 2013
 *      Author: settgast1
 */

#ifndef CONSTITUTIVEBASE_H_
#define CONSTITUTIVEBASE_H_

#include "legacy/Common/Common.h"
#include "legacy/ArrayT/bufvector.h"
#include <sstream>
//#include "../IO/ticpp/HierarchicalDataNode.h.old"

namespace TICPP
{
class HierarchicalDataNode;
}

namespace ConstitutiveValueResolution
{
template<typename T>
inline int ValueType(){ throw GPException("Unrecognized type"); }
template<>
inline int ValueType<int>(){return 0;}
template<>
inline int ValueType<realT>(){return 1;}
template<>
inline int ValueType<R1Tensor>(){return 2;}
template<>
inline int ValueType<R2Tensor>(){return 3;}
template<>
inline int ValueType<R2SymTensor>(){return 4;}

template<typename T>
inline const char* PointerConversion(const T& object, const size_t offset)
{
  const char* ptr = (char*)(&object) + offset;
  return ptr;
}

template<typename T>
inline char* PointerConversion(T& object, const size_t offset)
{
  char* ptr = (char*)(&object) + offset;
  return ptr;
}
}

class ConstitutiveBase
{

public:
  ConstitutiveBase(): hasVariableParams(false) {}
  virtual ~ConstitutiveBase() {}

  inline bool VariableParameters() const { return hasVariableParams; }

  virtual void SetVariableParameters(const bool varParams, const localIndex newSize = 0) = 0;

  template< typename LeafClass >
  void SetVariableParametersFromDerived(const bool varParams, const localIndex newSize = 1);

  virtual void resize( const localIndex num ) = 0;

  virtual void ReadXML( TICPP::HierarchicalDataNode& node ) = 0;

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

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking =true ) = 0;
  virtual unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking = true) = 0;

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking =true ) = 0;
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking = true) = 0;

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer ) = 0;
  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer ) = 0;

  static void SubNames( const array<string>& names, const localIndex num, array<string>& subnames )
  {
    for( array<string>::const_iterator inames=names.begin() ; inames!=names.end() ; ++inames )
    {
      for( auto i=0 ; i<num ; ++i )
      {
        char lbl[256] = {0};
        sprintf(lbl, "%s_ip%02u", inames->c_str(), i);
        subnames.push_back( lbl );
      }

    }
  }

  virtual bool GetStateValues( const std::string& name, array<real64>& values ) const = 0;
  virtual bool GetParameterValues( const std::string& name, array<real64>& values ) const = 0;
  virtual bool SetStateValues( const std::string& name, const array<real64>& values ) = 0;
  virtual bool SetParameterValues( const std::string& name, const array<real64>& values ) = 0;

  template< typename LeafClass, typename ARRAY >
  bool GetStateValuesFromDerived( const std::string& name,
                                  ARRAY& values) const;

  template< typename LeafClass, typename ARRAY >
  bool GetParameterValuesFromDerived( const std::string& name,
                                      ARRAY& values) const;

  template< typename LeafClass, typename ARRAY >
  bool SetStateValuesFromDerived( const std::string& name,
                                  const ARRAY& values);

  template< typename LeafClass, typename ARRAY >
  bool SetParameterValuesFromDerived( const std::string& name,
                                      const ARRAY& values);


protected:
  bool hasVariableParams;

  template< typename LeafClass >
  void ReadXMLFromDerived( TICPP::HierarchicalDataNode& node );

  template< typename LeafClass >
  void ResizeFromDerived( const localIndex num );

  template< typename LeafClass >
  void InsertFromDerived( const localIndex num );

  template< typename LeafClass >
  void EraseFromDerived( const localIndex num );

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
  unsigned int PackFromDerived( const T_indices& localIndices,
                                char*& buffer,
                                const bool doBufferPacking ) const;

  template< typename LeafClass, typename T_indices >
  unsigned int UnpackFromDerived( const T_indices& localIndices,
                                  const char*& buffer );

  inline static size_t GetOffset( const std::string& name, const std::map<std::string, size_t>& offsets )
  {
    const size_t* tmp = stlMapLookupPointer( offsets, name );
    if(tmp == NULL)
      return std::numeric_limits<size_t>::max();
    else
      return *tmp;
  }

  template< typename LeafClass >
  size_t GetParameterOffsetFromDerived( const std::string& name, const int type ) const;

  template< typename LeafClass >
  size_t GetStateOffsetFromDerived( const std::string& name, const int type ) const;

  virtual void PreSetValues(const array<string>& ){}
  virtual void PostSetValues(const array<string>& ){}

public:
  void SetValues( const std::string& name, const array<real64>& values )
  {
    bool ok = SetParameterValues(name, values);
    if(!ok)
      ok = SetStateValues(name, values);
    if(!ok)
    {
      std::string err("Could not find parameter or state value: ");
      err += name;
      throw GPException(err);
    }
  }

  void SetValues( const array<string>& names, const array<array<real64> >& values )
  {
    if(names.size() != values.size())
      throw GPException("Attempt to SetValues with arrays of different lengths");
    PreSetValues(names);
    for( auto i = 0u ; i < names.size() ; i++)
      SetValues(names[i], values[i]);
    PostSetValues(names);
  }

  void GetValues( const std::string& name, array<real64>& values ) const
  {
    if(!(GetParameterValues(name, values) || GetStateValues(name, values)))
    {
      std::string err("Could not find parameter or state value: ");
      err += name;
      throw GPException(err);
    }
  }

private:
  ConstitutiveBase( const ConstitutiveBase& );
  ConstitutiveBase& operator=( const ConstitutiveBase& );

};



template< typename LeafClass >
void ConstitutiveBase::ReadXMLFromDerived( TICPP::HierarchicalDataNode& node )
{
  LeafClass& dthis = static_cast<LeafClass&>(*this);
  dthis.m_parameterData.resize(1);

  dthis.m_parameterData[0].ReadXML(node);

}

template< typename LeafClass >
void ConstitutiveBase::ResizeFromDerived( const localIndex num )
{
  LeafClass& dthis = static_cast<LeafClass&>(*this);

  dthis.m_stateData.resize( num );

  if( hasVariableParams )
    dthis.m_parameterData.resize( num );

}

template< typename LeafClass >
void ConstitutiveBase::EraseFromDerived( const localIndex num )
{
  LeafClass& dthis = static_cast<LeafClass&>(*this);

  dthis.m_stateData.Erase( num );

  if( hasVariableParams )
    dthis.m_parameterData.erase( dthis.m_parameterData.begin() + num );
}

template< typename LeafClass >
void ConstitutiveBase::InsertFromDerived( const localIndex num )
{
  LeafClass& dthis = static_cast<LeafClass&>(*this);
  //ResizeFromDerived<LeafClass>(dthis.m_stateData.size()+1, variableParams);
  {
    typename LeafClass::StateClass s;
    dthis.m_stateData.Insert( num, s);
  }
  if( hasVariableParams )
  {
    typename LeafClass::ParameterClass p;
    dthis.m_parameterData.insert( dthis.m_parameterData.begin() + num, p );
  }
}

template< typename LeafClass >
void ConstitutiveBase::GetVariableNamesFromDerived( array<string>& intVars,
                                                    array<string>& realVars,
                                                    array<string>& R1TensorVars,
                                                    array<string>& R2TensorVars,
                                                    array<string>& R2SymTensorVars ) const
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);
  if( hasVariableParams ) //dthis.m_parameterData.size() > 1 )
  {
    LeafClass::ParameterClass::GetVariableNames( intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars);
  }

  array<string> intVarsTemp;
  array<string> realVarsTemp;
  array<string> R1TensorVarsTemp;
  array<string> R2TensorVarsTemp;
  array<string> R2SymTensorVarsTemp;
  LeafClass::StateClass::GetVariableNames( intVarsTemp, realVarsTemp, R1TensorVarsTemp, R2TensorVarsTemp, R2SymTensorVarsTemp);

  SubNames( intVarsTemp, dthis.NumStateIndex1(), intVars );
  SubNames( realVarsTemp, dthis.NumStateIndex1(), realVars );
  SubNames( R1TensorVarsTemp, dthis.NumStateIndex1(), R1TensorVars );
  SubNames( R2TensorVarsTemp, dthis.NumStateIndex1(), R2TensorVars );
  SubNames( R2SymTensorVarsTemp, dthis.NumStateIndex1(), R2SymTensorVars );


}

template< typename LeafClass >
size_t ConstitutiveBase::GetParameterOffsetFromDerived( const std::string& name, const int type ) const
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);
  std::map<std::string, size_t> intOffsets;
  std::map<std::string, size_t> realOffsets;
  std::map<std::string, size_t> R1TensorOffsets;
  std::map<std::string, size_t> R2TensorOffsets;
  std::map<std::string, size_t> R2SymTensorOffsets;
  dthis.m_parameterData[0].GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets,
                                               R2TensorOffsets, R2SymTensorOffsets );
  size_t ret = std::numeric_limits<size_t>::max();
  switch(type)
  {
  case 0:
    ret = GetOffset(name, intOffsets);
    break;
  case 1:
    ret = GetOffset(name, realOffsets);
    break;
  case 2:
    ret = GetOffset(name, R1TensorOffsets);
    break;
  case 3:
    ret = GetOffset(name, R2TensorOffsets);
    break;
  case 4:
    ret = GetOffset(name, R2SymTensorOffsets);
    break;
  }
  return ret;
}

template< typename LeafClass >
size_t ConstitutiveBase::GetStateOffsetFromDerived( const std::string& name, const int type ) const
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);
  std::map<std::string, size_t> intOffsets;
  std::map<std::string, size_t> realOffsets;
  std::map<std::string, size_t> R1TensorOffsets;
  std::map<std::string, size_t> R2TensorOffsets;
  std::map<std::string, size_t> R2SymTensorOffsets;
  dthis.m_stateData(0,0).GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets,
                                             R2TensorOffsets, R2SymTensorOffsets );
  size_t ret = std::numeric_limits<size_t>::max();
  switch(type)
  {
  case 0:
    ret = GetOffset(name, intOffsets);
    break;
  case 1:
    ret = GetOffset(name, realOffsets);
    break;
  case 2:
    ret = GetOffset(name, R1TensorOffsets);
    break;
  case 3:
    ret = GetOffset(name, R2TensorOffsets);
    break;
  case 4:
    ret = GetOffset(name, R2SymTensorOffsets);
    break;
  }
  return ret;
}

template< typename LeafClass, typename ARRAY >
bool ConstitutiveBase::GetParameterValuesFromDerived( const std::string& name, ARRAY& values) const
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);
  const int type = ConstitutiveValueResolution::ValueType<typename ARRAY::value_type>();
  const size_t offset = dthis.GetParameterOffset(name, type);
  if(offset == std::numeric_limits<size_t>::max())
    return false;

  if(!hasVariableParams)
  {
    const localIndex i = 0;
    values.resize(dthis.m_stateData.Dimension(0), *((const typename ARRAY::value_type*)
                                                    ConstitutiveValueResolution::PointerConversion(dthis.m_parameterData[i], offset)));
  }
  else
  {
    values.resize(dthis.m_stateData.Dimension(0));
    for(localIndex i = 0 ; i < dthis.m_stateData.Dimension(0) ; i++)
      values[i] = *((const typename ARRAY::value_type*)(
                      ConstitutiveValueResolution::PointerConversion(dthis.m_parameterData[i], offset)));
  }
  return true;
}

template< typename LeafClass, typename ARRAY >
bool ConstitutiveBase::GetStateValuesFromDerived( const std::string& name, ARRAY& values) const
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);
  const int type = ConstitutiveValueResolution::ValueType<typename ARRAY::value_type>();
  const size_t offset = dthis.GetStateOffset(name, type);
  if (offset == std::numeric_limits<size_t>::max())
    return false;
  if (dthis.m_stateData.Dimension(1) < 1)
    throw GPException(
            "ConstitutiveBase::GetStateValuesFromDerived - Cannot get state values when the 1 dimension is zero!");

  values.resize(dthis.m_stateData.Dimension(0));
  for (localIndex i = 0 ; i < dthis.m_stateData.Dimension(0) ; i++)
  {
    values[i] = *((const typename ARRAY::value_type*)(ConstitutiveValueResolution::PointerConversion(
                                                        dthis.m_stateData(i, 0), offset)));
    for (localIndex j = 1 ; j < dthis.m_stateData.Dimension(1) ; j++)
      values[i] += *((const typename ARRAY::value_type*)(ConstitutiveValueResolution::PointerConversion(
                                                           dthis.m_stateData(i, j), offset)));
    values[i] *= 1.0 / dthis.m_stateData.Dimension(1);
  }
  return true;
}

template< typename LeafClass, typename ARRAY >
bool ConstitutiveBase::SetParameterValuesFromDerived( const std::string& name, const ARRAY& values)
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);
  const int type = ConstitutiveValueResolution::ValueType<typename ARRAY::value_type>();
  const size_t offset = dthis.GetParameterOffset(name, type);
  if(offset == std::numeric_limits<size_t>::max())
    return false;

  if(values.size() != dthis.NumStateIndex0())
    throw GPException("ConstitutiveBase::SetParameterValuesFromDerived - Value array not the same size as the state array");
  if(values.size() == 0)
    return true;
  if(values.size() > 1)
    SetVariableParameters(true, values.size());

  if(!hasVariableParams)
  {
    const localIndex i = 0;
    typename ARRAY::value_type& value = *((typename ARRAY::value_type *)
                                          ConstitutiveValueResolution::PointerConversion(dthis.m_parameterData[i], offset));
    value = values[0];
  }
  else
  {
    for(localIndex i = 0 ; i < dthis.m_stateData.Dimension(0) ; i++)
    {
      typename ARRAY::value_type& value = *((typename ARRAY::value_type *)
                                            ConstitutiveValueResolution::PointerConversion(dthis.m_parameterData[i], offset));
      value = values[i];
    }
  }
  return true;
}

template< typename LeafClass, typename ARRAY >
bool ConstitutiveBase::SetStateValuesFromDerived( const std::string& name, const ARRAY& values)
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);
  const int type = ConstitutiveValueResolution::ValueType<typename ARRAY::value_type>();
  const size_t offset = dthis.GetStateOffset(name, type);
  if(offset == std::numeric_limits<size_t>::max())
    return false;

  if (values.size() != dthis.NumStateIndex0())
    throw GPException(
            "ConstitutiveBase::SetStateValuesFromDerived - Value array not the same size as the state array");
  if (values.size() == 0)
    return true;

  for (localIndex i = 0 ; i < dthis.m_stateData.Dimension(0) ; i++)
  {
    for (localIndex j = 0 ; j < dthis.m_stateData.Dimension(1) ; j++)
    {
      typename ARRAY::value_type& tmp = *((typename ARRAY::value_type *)ConstitutiveValueResolution::PointerConversion(
                                            dthis.m_stateData(i, j), offset));
      tmp = values[i];
    }
  }
  return true;
}

template< typename LeafClass >
void ConstitutiveBase::SetVariableParametersFromDerived(const bool varParams, const localIndex newSize)
{
  if(varParams != hasVariableParams )
  {
    LeafClass& dthis = static_cast<LeafClass&>(*this);
    hasVariableParams = varParams;
    if(varParams)
    {
      if(newSize < 1)
        return;

      dthis.m_parameterData.resize(newSize);
      for(localIndex a = 1 ; a < newSize ; a++)
        memcpy(&dthis.m_parameterData[a], &dthis.m_parameterData[0], sizeof(typename LeafClass::ParameterClass));
    }
    else
    {
      dthis.m_parameterData.resize(1);
    }
  }
}


template< typename LeafClass >
void ConstitutiveBase::SerializeFromDerived( array<array<integer>*>& intVars,
                                             array<array<real64>*>& realVars,
                                             array<array<R1Tensor>*>& R1Vars,
                                             array<array<R2Tensor>*>& R2Vars,
                                             array<array<R2SymTensor>*>& R2SymVars ) const
{

  const LeafClass& dthis = static_cast<const LeafClass&>(*this);

  localIndex intVarCounts = 0;
  localIndex realVarCounts = 0;
  localIndex R1TensorVarCounts = 0;
  localIndex R2TensorVarCounts = 0;
  localIndex R2SymTensorVarCounts = 0;

  if( hasVariableParams )//dthis.m_parameterData.size() > 1 )
  {
    for( localIndex a=0 ; a<dthis.m_parameterData.size() ; ++a )
    {
      intVarCounts = 0;
      realVarCounts = 0;
      R1TensorVarCounts = 0;
      R2TensorVarCounts = 0;
      R2SymTensorVarCounts = 0;
      dthis.m_parameterData[a].Serialize( a, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                                          intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts);
    }
  }

  {
    typename LeafClass::StateArrayType::const_iterator istate = dthis.m_stateData.begin();
    for( localIndex a=0 ; a<dthis.NumStateIndex0() ; ++a )
    {
      for( localIndex b=0 ; b<dthis.NumStateIndex1() ; ++b )
      {
        localIndex intVarCountsTmp = intVarCounts + b;
        localIndex realVarCountsTmp = realVarCounts + b;
        localIndex R1TensorVarCountsTmp = R1TensorVarCounts + b;
        localIndex R2TensorVarCountsTmp = R2TensorVarCounts + b;
        localIndex R2SymTensorVarCountsTmp = R2SymTensorVarCounts + b;

        istate->Serialize( b, dthis.NumStateIndex1(), a, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                           intVarCountsTmp, realVarCountsTmp, R1TensorVarCountsTmp, R2TensorVarCountsTmp, R2SymTensorVarCountsTmp);
        ++istate;
      }
    }
  }
}

template< typename LeafClass >
void ConstitutiveBase::DeserializeFromDerived( const array<array<integer>*>& intVars,
                                               const array<array<real64>*>& realVars,
                                               const array<array<R1Tensor>*>& R1Vars,
                                               const array<array<R2Tensor>*>& R2Vars,
                                               const array<array<R2SymTensor>*>& R2SymVars  )
{
  LeafClass& dthis = static_cast<LeafClass&>(*this);

  localIndex intVarCounts = 0;
  localIndex realVarCounts = 0;
  localIndex R1TensorVarCounts = 0;
  localIndex R2TensorVarCounts = 0;
  localIndex R2SymTensorVarCounts = 0;

  if( hasVariableParams ) //dthis.m_parameterData.size() > 1 )
  {
    for( localIndex a=0 ; a<dthis.m_parameterData.size() ; ++a )
    {
      intVarCounts = 0;
      realVarCounts = 0;
      R1TensorVarCounts = 0;
      R2TensorVarCounts = 0;
      R2SymTensorVarCounts = 0;
      dthis.m_parameterData[a].Deserialize( a, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                                            intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts);
    }
  }

  {
    typename LeafClass::StateArrayType::iterator istate = dthis.m_stateData.begin();
    for( localIndex a=0 ; a<dthis.NumStateIndex0() ; ++a )
    {
      for( localIndex b=0 ; b<dthis.NumStateIndex1() ; ++b )
      {
        localIndex intVarCountsTmp = intVarCounts + b;
        localIndex realVarCountsTmp = realVarCounts + b;
        localIndex R1TensorVarCountsTmp = R1TensorVarCounts + b;
        localIndex R2TensorVarCountsTmp = R2TensorVarCounts + b;
        localIndex R2SymTensorVarCountsTmp = R2SymTensorVarCounts + b;

        istate->Deserialize( b, dthis.NumStateIndex1(), a, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                             intVarCountsTmp, realVarCountsTmp, R1TensorVarCountsTmp, R2TensorVarCountsTmp, R2SymTensorVarCountsTmp);
        ++istate;
      }
    }
  }

}



template< typename LeafClass, typename T_indices >
unsigned int ConstitutiveBase::PackFromDerived( const T_indices& localIndices, bufvector& buffer, const bool doBufferPacking ) const
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);

  unsigned int packedSize = 0;

  for( typename T_indices::const_iterator index=localIndices.begin() ; index!=localIndices.end() ; ++index )
  {
    const typename LeafClass::StateClass* const pstate = &(dthis.m_stateData.data()[(*index)*(dthis.NumStateIndex1())]);
    for( localIndex b=0 ; b<dthis.NumStateIndex1() ; ++b )
    {
      if( doBufferPacking )
      {
        packedSize += buffer.PackSerialObject(pstate[b]);
      }
      else
      {
        packedSize += sizeof(typename LeafClass::StateClass);
      }
    }
  }

  if( hasVariableParams )//dthis.m_parameterData.size() > 1 )
  {
    for( typename T_indices::const_iterator index=localIndices.begin() ; index!=localIndices.end() ; ++index )
    {
      if( doBufferPacking )
      {
        packedSize += buffer.PackSerialObject(dthis.m_parameterData[*index]);
      }
      else
      {
        packedSize += sizeof(typename LeafClass::ParameterClass);
      }
    }
  }

  return packedSize;
}


template< typename LeafClass, typename T_indices >
unsigned int ConstitutiveBase::PackFromDerived( const T_indices& localIndices, char*& buffer, const bool doBufferPacking ) const
{
  const LeafClass& dthis = static_cast<const LeafClass&>(*this);

  unsigned int packedSize = 0;

  for( typename T_indices::const_iterator index=localIndices.begin() ; index!=localIndices.end() ; ++index )
  {
    const typename LeafClass::StateClass* const pstate = &(dthis.m_stateData.data()[(*index)*(dthis.NumStateIndex1())]);
    for( localIndex b=0 ; b<dthis.NumStateIndex1() ; ++b )
    {
      if( doBufferPacking )
      {
        packedSize += bufvector::PackSerialObject(buffer, pstate[b]);
      }
      else
      {
        packedSize += sizeof(typename LeafClass::StateClass);
      }
    }
  }

  if( hasVariableParams )//dthis.m_parameterData.size() > 1 )
  {
    for( typename T_indices::const_iterator index=localIndices.begin() ; index!=localIndices.end() ; ++index )
    {
      if( doBufferPacking )
      {
        packedSize += bufvector::PackSerialObject(buffer, dthis.m_parameterData[*index]);
      }
      else
      {
        packedSize += sizeof(typename LeafClass::ParameterClass);
      }
    }
  }

  return packedSize;
}

template< typename LeafClass, typename T_indices >
unsigned int ConstitutiveBase::UnpackFromDerived( const T_indices& localIndices, const char*& buffer )
{
  LeafClass& dthis = static_cast<LeafClass&>(*this);

  unsigned int unpackedSize = 0;
  for( typename T_indices::const_iterator index=localIndices.begin() ; index!=localIndices.end() ; ++index )
  {
    typename LeafClass::StateClass* const pstate = &(dthis.m_stateData.data()[(*index)*(dthis.NumStateIndex1())]);
    for( localIndex b=0 ; b<dthis.NumStateIndex1() ; ++b )
    {
      unpackedSize += bufvector::UnpackSerialObject( buffer, pstate[b] );
    }
  }

  if( hasVariableParams ) //dthis.m_parameterData.size() > 1 )
  {
    for( typename T_indices::const_iterator index=localIndices.begin() ; index!=localIndices.end() ; ++index )
    {
      unpackedSize += bufvector::UnpackSerialObject( buffer, dthis.m_parameterData[*index] );
    }
  }
  return unpackedSize;
}



#endif /* CONSTITUTIVEBASE_H_ */
