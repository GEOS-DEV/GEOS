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


#ifndef LINEARELASTICINTERMEDIATE_H_
#define LINEARELASTICINTERMEDIATE_H_

#include "Utilities/GeometryUtilities.h"
#include "MaterialBase.h"

/*
 * LinearElasticIntermediate.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */



//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearElasticIntermediateParameterData : public MaterialBaseParameterData
{

public:

  typedef MaterialBaseParameterData base;
  realT init_bulkModulus;


  LinearElasticIntermediateParameterData():
    base(),
    init_bulkModulus(0)
  {}

  LinearElasticIntermediateParameterData( const LinearElasticIntermediateParameterData& source):
    base( source ),
    init_bulkModulus(source.init_bulkModulus)
  {}

  ~LinearElasticIntermediateParameterData() {}
  friend class ConstitutiveBase;
  friend class LinearElasticIntermediate;

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
    realVarCounts += 1;

  }

  static void GetVariableNames( array<string>& intNames,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>& R2TensorNames,
                                array<string>& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("BulkModulus");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                   std::map<std::string, size_t>& realOffsets,
                                   std::map<std::string, size_t>& R1TensorOffsets,
                                   std::map<std::string, size_t>& R2TensorOffsets,
                                   std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["BulkModulus"] = (char*)(&init_bulkModulus) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                  std::map<std::string, realT>& realValues,
                                  std::map<std::string, R1Tensor>& R1TensorValues,
                                  std::map<std::string, R2Tensor>& R2TensorValues,
                                  std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["BulkModulus"] = init_bulkModulus;
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
    (*(realVars[realVarCounts]))[index] = init_bulkModulus; realVarCounts++;
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
    init_bulkModulus = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline LinearElasticIntermediateParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    init_bulkModulus *= factor;
    return *this;
  }

  inline LinearElasticIntermediateParameterData&
  operator=(const LinearElasticIntermediateParameterData& datum)
  {
    base::operator=(datum);
    init_bulkModulus = datum.init_bulkModulus;
    return *this;
  }

  inline LinearElasticIntermediateParameterData&
  operator+=(const LinearElasticIntermediateParameterData& datum)
  {
    base::operator+=(datum);
    init_bulkModulus += datum.init_bulkModulus;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, LinearElasticIntermediateParameterData& p0, LinearElasticIntermediateParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, init_bulkModulus, p0.init_bulkModulus, p1.init_bulkModulus);

  }

  void MapFromRegion(const LinearElasticIntermediateParameterData& p0, const LinearElasticIntermediateParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.init_bulkModulus, p1.init_bulkModulus, fct0, fct1, init_bulkModulus);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    MaterialBaseParameterData::ReadXML( node );
    init_bulkModulus = node.GetAttributeOrDefault("BulkModulus", 0.0);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearElasticIntermediateStateData : public MaterialBaseStateData
{

public:

  typedef MaterialBaseStateData base;
  realT density;
  realT ElasticBulkModulus;
  realT ElasticShearModulus;


  LinearElasticIntermediateStateData():
    base(),
    density(0),
    ElasticBulkModulus(0),
    ElasticShearModulus(0)
  {}

  LinearElasticIntermediateStateData( const LinearElasticIntermediateStateData& source):
    base( source ),
    density(source.density),
    ElasticBulkModulus(source.ElasticBulkModulus),
    ElasticShearModulus(source.ElasticShearModulus)
  {}

  ~LinearElasticIntermediateStateData() {}
  friend class ConstitutiveBase;
  friend class LinearElasticIntermediate;

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
    realVarCounts += 3;

  }

  static void GetVariableNames( array<string>& intNames,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>& R2TensorNames,
                                array<string>& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("density");
    realNames.push_back("ElasticBulkModulus");
    realNames.push_back("ElasticShearModulus");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                   std::map<std::string, size_t>& realOffsets,
                                   std::map<std::string, size_t>& R1TensorOffsets,
                                   std::map<std::string, size_t>& R2TensorOffsets,
                                   std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["density"] = (char*)(&density) - (char*)this;
    realOffsets["ElasticBulkModulus"] = (char*)(&ElasticBulkModulus) - (char*)this;
    realOffsets["ElasticShearModulus"] = (char*)(&ElasticShearModulus) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                  std::map<std::string, realT>& realValues,
                                  std::map<std::string, R1Tensor>& R1TensorValues,
                                  std::map<std::string, R2Tensor>& R2TensorValues,
                                  std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["density"] = density;
    realValues["ElasticBulkModulus"] = ElasticBulkModulus;
    realValues["ElasticShearModulus"] = ElasticShearModulus;
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
    (*(realVars[realVarCounts]))[elemNum] = density; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = ElasticBulkModulus; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = ElasticShearModulus; realVarCounts += stride;
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
    density = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    ElasticBulkModulus = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    ElasticShearModulus = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
  }
  inline LinearElasticIntermediateStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    density *= factor;
    ElasticBulkModulus *= factor;
    ElasticShearModulus *= factor;
    return *this;
  }

  inline LinearElasticIntermediateStateData&
  operator=(const LinearElasticIntermediateStateData& datum)
  {
    base::operator=(datum);
    density = datum.density;
    ElasticBulkModulus = datum.ElasticBulkModulus;
    ElasticShearModulus = datum.ElasticShearModulus;
    return *this;
  }

  inline LinearElasticIntermediateStateData&
  operator+=(const LinearElasticIntermediateStateData& datum)
  {
    base::operator+=(datum);
    density += datum.density;
    ElasticBulkModulus += datum.ElasticBulkModulus;
    ElasticShearModulus += datum.ElasticShearModulus;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, LinearElasticIntermediateStateData& p0, LinearElasticIntermediateStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, density, p0.density, p1.density);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ElasticBulkModulus, p0.ElasticBulkModulus, p1.ElasticBulkModulus);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ElasticShearModulus, p0.ElasticShearModulus, p1.ElasticShearModulus);

  }

  void MapFromRegion(const LinearElasticIntermediateStateData& p0, const LinearElasticIntermediateStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.density, p1.density, fct0, fct1, density);
    GeometryUtilities::MapFromRegion(p0.ElasticBulkModulus, p1.ElasticBulkModulus, fct0, fct1, ElasticBulkModulus);
    GeometryUtilities::MapFromRegion(p0.ElasticShearModulus, p1.ElasticShearModulus, fct0, fct1, ElasticShearModulus);

  }



};


//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearElasticIntermediate : public MaterialBase
{
public:

  typedef LinearElasticIntermediateParameterData ParameterClass;
  typedef LinearElasticIntermediateStateData     StateClass;


  LinearElasticIntermediate( const int paramSize, const int stateSize );

  virtual ~LinearElasticIntermediate();

  virtual void ReadXML( TICPP::HierarchicalDataNode& node ) = 0;

  virtual void resize( const localIndex num ) = 0;

  virtual void resize( const localIndex num0,
                       const localIndex num1 ) = 0;

  virtual void insert( const localIndex num ) = 0;

  virtual void erase( const localIndex num ) = 0;

  virtual const LinearElasticIntermediateStateData* StateData( const localIndex index0,
                                                               const localIndex index1 ) const = 0;
  virtual       LinearElasticIntermediateStateData* StateData( const localIndex index0,
                                                               const localIndex index1 )  = 0;

  virtual const LinearElasticIntermediateParameterData* ParameterData( const localIndex index ) const = 0;
  virtual       LinearElasticIntermediateParameterData* ParameterData( const localIndex index ) = 0;

  inline void IncrementPtr( const LinearElasticIntermediateStateData* ptr ) const
  {
    ptr = reinterpret_cast<const LinearElasticIntermediateStateData*>( reinterpret_cast<const char*>(ptr) + m_stateSize );
  }

  inline void IncrementPtr( const LinearElasticIntermediateParameterData* ptr ) const
  {
    ptr = reinterpret_cast<const LinearElasticIntermediateParameterData*>( reinterpret_cast<const char*>(ptr) + m_paramSize );
  }

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

  virtual void ZeroStates() = 0;

  virtual localIndex NumStateIndex0() const = 0;
  virtual localIndex NumStateIndex1() const = 0;

  virtual localIndex NumParameterIndex0() const = 0;
  virtual localIndex NumParameterIndex1() const = 0;



private:
  LinearElasticIntermediate();
  LinearElasticIntermediate( const LinearElasticIntermediate& );
  LinearElasticIntermediate& operator=( const LinearElasticIntermediate& );


};
#endif /* LINEARELASTICINTERMEDIATE_H_ */
