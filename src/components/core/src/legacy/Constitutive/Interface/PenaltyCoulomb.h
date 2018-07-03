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


#ifndef PENALTYCOULOMB_H_
#define PENALTYCOULOMB_H_

#include "Utilities/GeometryUtilities.h"
#include "PenaltyCoulombIntermediate.h"

/*
 * PenaltyCoulomb.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */



//**********************************************************************************************************************
//**********************************************************************************************************************


class PenaltyCoulombParameterData : public PenaltyCoulombIntermediateParameterData
{

public:

  typedef PenaltyCoulombIntermediateParameterData base;


  PenaltyCoulombParameterData():
    base()
  {}

  PenaltyCoulombParameterData( const PenaltyCoulombParameterData& source):
    base( source )
  {}

  ~PenaltyCoulombParameterData() {}
  friend class ConstitutiveBase;
  friend class PenaltyCoulomb;

  static void GetVariableCounts( localIndex& intVarCounts,
                                 localIndex& realVarCounts,
                                 localIndex& R1TensorVarCounts,
                                 localIndex& R2TensorVarCounts,
                                 localIndex& R2SymTensorVarCounts )
  {
    base::GetVariableCounts( intVarCounts,
                             realVarCounts,
                             R1TensorVarCounts,
                             R2TensorVarCounts,
                             R2SymTensorVarCounts );

  }

  static void GetVariableNames( array<string>& intNames,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>& R2TensorNames,
                                array<string>& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                   std::map<std::string, size_t>& realOffsets,
                                   std::map<std::string, size_t>& R1TensorOffsets,
                                   std::map<std::string, size_t>& R2TensorOffsets,
                                   std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                  std::map<std::string, realT>& realValues,
                                  std::map<std::string, R1Tensor>& R1TensorValues,
                                  std::map<std::string, R2Tensor>& R2TensorValues,
                                  std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
  }

  void Serialize(const localIndex index,
                 array<array<integer>*>& intVars,
                 array<array<real64>*>& realVars,
                 array<array<R1Tensor>*>& R1Vars,
                 array<array<R2Tensor>*>& R2Vars,
                 array<array<R2SymTensor>*>& R2SymVars,
                 localIndex& intVarCounts,
                 localIndex& realVarCounts,
                 localIndex& R1TensorVarCounts,
                 localIndex& R2TensorVarCounts,
                 localIndex& R2SymTensorVarCounts  ) const
  {
    base::Serialize(index, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
  }


  void  Deserialize( const localIndex index,
                     const array<array<integer>*>& intVars,
                     const array<array<real64>*>& realVars,
                     const array<array<R1Tensor>*>& R1Vars,
                     const array<array<R2Tensor>*>& R2Vars,
                     const array<array<R2SymTensor>*>& R2SymVars,
                     localIndex& intVarCounts,
                     localIndex& realVarCounts,
                     localIndex& R1TensorVarCounts,
                     localIndex& R2TensorVarCounts,
                     localIndex& R2SymTensorVarCounts )
  {
    base::Deserialize(index, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                      intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
  }
  inline PenaltyCoulombParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    return *this;
  }

  inline PenaltyCoulombParameterData&
  operator=(const PenaltyCoulombParameterData& datum)
  {
    base::operator=(datum);
    return *this;
  }

  inline PenaltyCoulombParameterData&
  operator+=(const PenaltyCoulombParameterData& datum)
  {
    base::operator+=(datum);
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, PenaltyCoulombParameterData& p0, PenaltyCoulombParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);

  }

  void MapFromRegion(const PenaltyCoulombParameterData& p0, const PenaltyCoulombParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    PenaltyCoulombIntermediateParameterData::ReadXML( node );
    PostReadXML( node );

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class PenaltyCoulombStateData : public PenaltyCoulombIntermediateStateData
{

public:

  typedef PenaltyCoulombIntermediateStateData base;


  PenaltyCoulombStateData():
    base()
  {}

  PenaltyCoulombStateData( const PenaltyCoulombStateData& source):
    base( source )
  {}

  ~PenaltyCoulombStateData() {}
  friend class ConstitutiveBase;
  friend class PenaltyCoulomb;

  static void GetVariableCounts( localIndex& intVarCounts,
                                 localIndex& realVarCounts,
                                 localIndex& R1TensorVarCounts,
                                 localIndex& R2TensorVarCounts,
                                 localIndex& R2SymTensorVarCounts )
  {
    base::GetVariableCounts( intVarCounts,
                             realVarCounts,
                             R1TensorVarCounts,
                             R2TensorVarCounts,
                             R2SymTensorVarCounts );

  }

  static void GetVariableNames( array<string>& intNames,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>& R2TensorNames,
                                array<string>& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                   std::map<std::string, size_t>& realOffsets,
                                   std::map<std::string, size_t>& R1TensorOffsets,
                                   std::map<std::string, size_t>& R2TensorOffsets,
                                   std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                  std::map<std::string, realT>& realValues,
                                  std::map<std::string, R1Tensor>& R1TensorValues,
                                  std::map<std::string, R2Tensor>& R2TensorValues,
                                  std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
  }

  void Serialize(const localIndex index,
                 const unsigned int stride,
                 const localIndex elemNum,
                 array<array<integer>*>& intVars,
                 array<array<real64>*>& realVars,
                 array<array<R1Tensor>*>& R1Vars,
                 array<array<R2Tensor>*>& R2Vars,
                 array<array<R2SymTensor>*>& R2SymVars,
                 localIndex& intVarCounts,
                 localIndex& realVarCounts,
                 localIndex& R1TensorVarCounts,
                 localIndex& R2TensorVarCounts,
                 localIndex& R2SymTensorVarCounts  ) const
  {
    base::Serialize(index, stride, elemNum, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
  }


  void  Deserialize( const localIndex index,
                     const unsigned int stride,
                     const localIndex elemNum,
                     const array<array<integer>*>& intVars,
                     const array<array<real64>*>& realVars,
                     const array<array<R1Tensor>*>& R1Vars,
                     const array<array<R2Tensor>*>& R2Vars,
                     const array<array<R2SymTensor>*>& R2SymVars,
                     localIndex& intVarCounts,
                     localIndex& realVarCounts,
                     localIndex& R1TensorVarCounts,
                     localIndex& R2TensorVarCounts,
                     localIndex& R2SymTensorVarCounts )
  {
    base::Deserialize(index, stride, elemNum, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                      intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
  }
  inline PenaltyCoulombStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    return *this;
  }

  inline PenaltyCoulombStateData&
  operator=(const PenaltyCoulombStateData& datum)
  {
    base::operator=(datum);
    return *this;
  }

  inline PenaltyCoulombStateData&
  operator+=(const PenaltyCoulombStateData& datum)
  {
    base::operator+=(datum);
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, PenaltyCoulombStateData& p0, PenaltyCoulombStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);

  }

  void MapFromRegion(const PenaltyCoulombStateData& p0, const PenaltyCoulombStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************
class PenaltyCoulomb : public PenaltyCoulombIntermediate
{
public:

  typedef PenaltyCoulombParameterData ParameterClass;
  typedef PenaltyCoulombStateData     StateClass;

  typedef array<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return "PenaltyCoulomb"; }

  PenaltyCoulomb();
  virtual ~PenaltyCoulomb();

  inline void MapToRegion(const realT fctNormal, const realT fct0, const realT fct1,
                          const localIndex from0, const localIndex from1,
                          StateClass& s0, StateClass& s1)
  {
    StateData(from0, from1)->MapToRegion(fctNormal, fct0, fct1, s0, s1);
    //ParameterData(from)->MapToRegion(fctNormal, fct0, fct1, p0, p1);
  }

  inline void MapFromRegion(const realT fct0, const realT fct1,
                            const StateClass& s0, const StateClass& s1,
                            const localIndex to0, const localIndex to1)
  {
    StateData(to0, to1)->MapFromRegion(s0, s1, fct0, fct1);
    //ParameterData(to)->MapFromRegion(p0, p1, fct0, fct1);
  }

  inline void MapToRegion(const realT fctNormal, const realT fct0, const realT fct1,
                          const localIndex from,
                          ParameterClass& p0, ParameterClass& p1)
  {
    //StateData(from0, from1)->MapToRegion(fctNormal, fct0, fct1, s0, s1);
    ParameterData(from)->MapToRegion(fctNormal, fct0, fct1, p0, p1);
  }

  inline void MapFromRegion(const realT fct0, const realT fct1,
                            const ParameterClass& p0, const ParameterClass& p1,
                            const localIndex to)
  {
    //StateData(to0, to1)->MapFromRegion(s0, s1, fct0, fct1);
    ParameterData(to)->MapFromRegion(p0, p1, fct0, fct1);
  }


  virtual void ZeroStates()
  {
    for(localIndex j = 0 ; j < m_stateData.Dimension(1) ; j++)
    {
      for(localIndex i = 0 ; i < m_stateData.Dimension(0) ; i++)
      {
        m_stateData(i,j) *= 0.0;
      }
    }
  }

  virtual void SetVariableParameters(const bool varParams, const localIndex newSize = 0)
  { SetVariableParametersFromDerived<PenaltyCoulomb>(varParams, newSize); }

  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<PenaltyCoulomb>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<PenaltyCoulomb>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<PenaltyCoulomb>( num0 );
  }

  virtual void insert( const localIndex num )
  { InsertFromDerived<PenaltyCoulomb>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<PenaltyCoulomb>( num ); }

  void GetVariableNames( array<string>& intVars, array<string>& realVars, array<string>& R1TensorVars, array<string>& R2TensorVars,
                         array<string>& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<PenaltyCoulomb>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<PenaltyCoulomb>(name, type); }

  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<PenaltyCoulomb>(name, type ); }

  bool GetStateValues( const std::string& name, array<real64>& values ) const
  { return GetStateValuesFromDerived<PenaltyCoulomb>(name, values); }

  bool GetParameterValues( const std::string& name, array<real64>& values ) const
  { return GetParameterValuesFromDerived<PenaltyCoulomb>(name, values); }

  bool SetStateValues( const std::string& name, const array<real64>& values )
  { return SetStateValuesFromDerived<PenaltyCoulomb>(name, values); }

  bool SetParameterValues( const std::string& name, const array<real64>& values )
  { return SetParameterValuesFromDerived<PenaltyCoulomb>(name, values); }

  virtual void Serialize( array<array<integer>*>& intVars, array<array<real64>*>& realVars, array<array<R1Tensor>*>& R1Vars, array<array<R2Tensor>*>& R2Vars,
                          array<array<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<PenaltyCoulomb>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const array<array<integer>*>& intVars, const array<array<real64>*>& realVars, const array<array<R1Tensor>*>& R1Vars,
                            const array<array<R2Tensor>*>& R2Vars, const array<array<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<PenaltyCoulomb>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<PenaltyCoulomb>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<PenaltyCoulomb>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<PenaltyCoulomb>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<PenaltyCoulomb>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<PenaltyCoulomb>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<PenaltyCoulomb>( localIndices, buffer ); }

  const StateClass* StateData( const localIndex index0, const localIndex index1 ) const
  { return &( m_stateData(index0,index1) );  }
  StateClass* StateData( const localIndex index0, const localIndex index1 )
  { return &( m_stateData(index0,index1) );  }

  const ParameterClass* ParameterData( const localIndex index ) const
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  ParameterClass* ParameterData( const localIndex index )
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }



private:
  PenaltyCoulomb(const PenaltyCoulomb&);
  PenaltyCoulomb& operator=(const PenaltyCoulomb&);


};
#endif /* PENALTYCOULOMB_H_ */
