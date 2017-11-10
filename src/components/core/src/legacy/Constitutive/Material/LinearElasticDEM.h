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

#ifndef LINEARELASTICDEM_H_
#define LINEARELASTICDEM_H_

#include "Utilities/GeometryUtilities.h"
#include "MaterialBase.h"

/*
 * LinearElasticDEM.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearElasticDEMParameterData : public MaterialBaseParameterData
{

public:

  typedef MaterialBaseParameterData base;
  realT yieldStrength;
  realT cor;
  realT velHalf;
  realT cement;
  realT surfaceEnergy;


  LinearElasticDEMParameterData():
    base(),
    yieldStrength(0),
    cor(0),
    velHalf(0),
    cement(0),
    surfaceEnergy(0)
  {}

  LinearElasticDEMParameterData( const LinearElasticDEMParameterData& source):
    base( source ),
    yieldStrength(source.yieldStrength),
    cor(source.cor),
    velHalf(source.velHalf),
    cement(source.cement),
    surfaceEnergy(source.surfaceEnergy)
  {}

  ~LinearElasticDEMParameterData() {}
  friend class ConstitutiveBase;
  friend class LinearElasticDEM;

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
    realVarCounts += 5;

  }

  static void GetVariableNames( array<string>& intNames,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>& R2TensorNames,
                                array<string>& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("yieldStrength");
    realNames.push_back("coefficientOfRestitution");
    realNames.push_back("velocityOfCoROfHalf");
    realNames.push_back("cementBondStrength");
    realNames.push_back("surfaceEnergy");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["yieldStrength"] = (char*)(&yieldStrength) - (char*)this;
    realOffsets["coefficientOfRestitution"] = (char*)(&cor) - (char*)this;
    realOffsets["velocityOfCoROfHalf"] = (char*)(&velHalf) - (char*)this;
    realOffsets["cementBondStrength"] = (char*)(&cement) - (char*)this;
    realOffsets["surfaceEnergy"] = (char*)(&surfaceEnergy) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["yieldStrength"] = yieldStrength;
    realValues["coefficientOfRestitution"] = cor;
    realValues["velocityOfCoROfHalf"] = velHalf;
    realValues["cementBondStrength"] = cement;
    realValues["surfaceEnergy"] = surfaceEnergy;
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
    (*(realVars[realVarCounts]))[index] = yieldStrength; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = cor; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = velHalf; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = cement; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = surfaceEnergy; realVarCounts++;
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
    yieldStrength = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    cor = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    velHalf = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    cement = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    surfaceEnergy = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline LinearElasticDEMParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    yieldStrength *= factor;
    cor *= factor;
    velHalf *= factor;
    cement *= factor;
    surfaceEnergy *= factor;
    return *this;
  }

  inline LinearElasticDEMParameterData&
  operator=(const LinearElasticDEMParameterData& datum)
  {
    base::operator=(datum);
    yieldStrength = datum.yieldStrength;
    cor = datum.cor;
    velHalf = datum.velHalf;
    cement = datum.cement;
    surfaceEnergy = datum.surfaceEnergy;
    return *this;
  }

  inline LinearElasticDEMParameterData&
  operator+=(const LinearElasticDEMParameterData& datum)
  {
    base::operator+=(datum);
    yieldStrength += datum.yieldStrength;
    cor += datum.cor;
    velHalf += datum.velHalf;
    cement += datum.cement;
    surfaceEnergy += datum.surfaceEnergy;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, LinearElasticDEMParameterData& p0, LinearElasticDEMParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, yieldStrength, p0.yieldStrength, p1.yieldStrength);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, cor, p0.cor, p1.cor);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, velHalf, p0.velHalf, p1.velHalf);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, cement, p0.cement, p1.cement);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, surfaceEnergy, p0.surfaceEnergy, p1.surfaceEnergy);

  }

  void MapFromRegion(const LinearElasticDEMParameterData& p0, const LinearElasticDEMParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.yieldStrength, p1.yieldStrength, fct0, fct1, yieldStrength);
    GeometryUtilities::MapFromRegion(p0.cor, p1.cor, fct0, fct1, cor);
    GeometryUtilities::MapFromRegion(p0.velHalf, p1.velHalf, fct0, fct1, velHalf);
    GeometryUtilities::MapFromRegion(p0.cement, p1.cement, fct0, fct1, cement);
    GeometryUtilities::MapFromRegion(p0.surfaceEnergy, p1.surfaceEnergy, fct0, fct1, surfaceEnergy);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    MaterialBaseParameterData::ReadXML( node );
    yieldStrength = node.GetAttributeOrDefault("yieldStrength", 0.0);
    cor = node.GetAttributeOrDefault("coefficientOfRestitution", 0.0);
    velHalf = node.GetAttributeOrDefault("velocityOfCoROfHalf", 0.0);
    cement = node.GetAttributeOrDefault("cementBondStrength", 0.0);
    surfaceEnergy = node.GetAttributeOrDefault("surfaceEnergy", 0.0);
    PostReadXML( node );

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearElasticDEMStateData : public MaterialBaseStateData
{

public:

  typedef MaterialBaseStateData base;


  LinearElasticDEMStateData():
    base()
  {}

  LinearElasticDEMStateData( const LinearElasticDEMStateData& source):
    base( source )
  {}

  ~LinearElasticDEMStateData() {}
  friend class ConstitutiveBase;
  friend class LinearElasticDEM;

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
  inline LinearElasticDEMStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    return *this;
  }

  inline LinearElasticDEMStateData&
  operator=(const LinearElasticDEMStateData& datum)
  {
    base::operator=(datum);
    return *this;
  }

  inline LinearElasticDEMStateData&
  operator+=(const LinearElasticDEMStateData& datum)
  {
    base::operator+=(datum);
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, LinearElasticDEMStateData& p0, LinearElasticDEMStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);

  }

  void MapFromRegion(const LinearElasticDEMStateData& p0, const LinearElasticDEMStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************
class LinearElasticDEM: public MaterialBase
{
public:

  typedef LinearElasticDEMParameterData ParameterClass;
  typedef LinearElasticDEMStateData     StateClass;

  typedef array<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return "LinearElasticDEM"; }

  LinearElasticDEM();
  virtual ~LinearElasticDEM();

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
    for(localIndex j = 0; j < m_stateData.Dimension(1); j++)
    {
      for(localIndex i = 0; i < m_stateData.Dimension(0); i++)
      {
        m_stateData(i,j) *= 0.0;
      }
    }
  }
  
  virtual void SetVariableParameters(const bool varParams, const localIndex newSize = 0)
  { SetVariableParametersFromDerived<LinearElasticDEM>(varParams, newSize); }

  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<LinearElasticDEM>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<LinearElasticDEM>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<LinearElasticDEM>( num0 );
  }
  
  virtual void insert( const localIndex num )
  { InsertFromDerived<LinearElasticDEM>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<LinearElasticDEM>( num ); }
 
  void GetVariableNames( array<string>& intVars, array<string>& realVars, array<string>& R1TensorVars, array<string>& R2TensorVars, array<string>& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<LinearElasticDEM>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<LinearElasticDEM>(name, type); }
  
  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<LinearElasticDEM>(name, type ); }
  
  bool GetStateValues( const std::string& name, array<real64>& values ) const
  { return GetStateValuesFromDerived<LinearElasticDEM>(name, values); }

  bool GetParameterValues( const std::string& name, array<real64>& values ) const
  { return GetParameterValuesFromDerived<LinearElasticDEM>(name, values); }

  bool SetStateValues( const std::string& name, const array<real64>& values )
  { return SetStateValuesFromDerived<LinearElasticDEM>(name, values); }

  bool SetParameterValues( const std::string& name, const array<real64>& values )
  { return SetParameterValuesFromDerived<LinearElasticDEM>(name, values); }

  virtual void Serialize( array<array<integer>*>& intVars, array<array<real64>*>& realVars, array<array<R1Tensor>*>& R1Vars, array<array<R2Tensor>*>& R2Vars, array<array<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<LinearElasticDEM>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const array<array<integer>*>& intVars, const array<array<real64>*>& realVars, const array<array<R1Tensor>*>& R1Vars, const array<array<R2Tensor>*>& R2Vars, const array<array<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<LinearElasticDEM>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticDEM>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticDEM>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticDEM>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticDEM>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<LinearElasticDEM>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<LinearElasticDEM>( localIndices, buffer ); }

  const StateClass* StateData( const localIndex index0, const localIndex index1 ) const
  { return &( m_stateData(index0,index1) );  }
  StateClass* StateData( const localIndex index0, const localIndex index1 )
  { return &( m_stateData(index0,index1) );  }

  const ParameterClass* ParameterData( const localIndex index ) const
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  ParameterClass* ParameterData( const localIndex index )
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  


private:
  LinearElasticDEM(const LinearElasticDEM&);
  LinearElasticDEM& operator=(const LinearElasticDEM&);
  

};
#endif /* LINEARELASTICDEM_H_ */
