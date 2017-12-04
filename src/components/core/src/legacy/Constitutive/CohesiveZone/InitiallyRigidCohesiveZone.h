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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INITIALLYRIGIDCOHESIVEZONE_H_
#define INITIALLYRIGIDCOHESIVEZONE_H_

#include "Utilities/GeometryUtilities.h"
#include "CohesiveZoneBase.h"

/*
 * InitiallyRigidCohesiveZone.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */



//**********************************************************************************************************************
//**********************************************************************************************************************


class InitiallyRigidCohesiveZoneParameterData : public CohesiveZoneBaseParameterData
{

public:

  typedef CohesiveZoneBaseParameterData base;
  realT failStress;
  realT failGap;


  InitiallyRigidCohesiveZoneParameterData():
    base(),
    failStress(0),
    failGap(0)
  {}

  InitiallyRigidCohesiveZoneParameterData( const InitiallyRigidCohesiveZoneParameterData& source):
    base( source ),
    failStress(source.failStress),
    failGap(source.failGap)
  {}

  ~InitiallyRigidCohesiveZoneParameterData() {}
  friend class ConstitutiveBase;
  friend class InitiallyRigidCohesiveZone;

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
    realVarCounts += 2;

  }

  static void GetVariableNames( array<string>& intNames,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>& R2TensorNames,
                                array<string>& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("failstress");
    realNames.push_back("failgap");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                   std::map<std::string, size_t>& realOffsets,
                                   std::map<std::string, size_t>& R1TensorOffsets,
                                   std::map<std::string, size_t>& R2TensorOffsets,
                                   std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["failstress"] = (char*)(&failStress) - (char*)this;
    realOffsets["failgap"] = (char*)(&failGap) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                  std::map<std::string, realT>& realValues,
                                  std::map<std::string, R1Tensor>& R1TensorValues,
                                  std::map<std::string, R2Tensor>& R2TensorValues,
                                  std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["failstress"] = failStress;
    realValues["failgap"] = failGap;
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
    (*(realVars[realVarCounts]))[index] = failStress; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = failGap; realVarCounts++;
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
    failStress = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    failGap = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline InitiallyRigidCohesiveZoneParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    failStress *= factor;
    failGap *= factor;
    return *this;
  }

  inline InitiallyRigidCohesiveZoneParameterData&
  operator=(const InitiallyRigidCohesiveZoneParameterData& datum)
  {
    base::operator=(datum);
    failStress = datum.failStress;
    failGap = datum.failGap;
    return *this;
  }

  inline InitiallyRigidCohesiveZoneParameterData&
  operator+=(const InitiallyRigidCohesiveZoneParameterData& datum)
  {
    base::operator+=(datum);
    failStress += datum.failStress;
    failGap += datum.failGap;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, InitiallyRigidCohesiveZoneParameterData& p0, InitiallyRigidCohesiveZoneParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, failStress, p0.failStress, p1.failStress);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, failGap, p0.failGap, p1.failGap);

  }

  void MapFromRegion(const InitiallyRigidCohesiveZoneParameterData& p0, const InitiallyRigidCohesiveZoneParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.failStress, p1.failStress, fct0, fct1, failStress);
    GeometryUtilities::MapFromRegion(p0.failGap, p1.failGap, fct0, fct1, failGap);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    CohesiveZoneBaseParameterData::ReadXML( node );
    failStress = node.GetAttributeOrDefault("failstress", 0.0);
    failGap = node.GetAttributeOrDefault("failgap", 0.0);
    PostReadXML( node );

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class InitiallyRigidCohesiveZoneStateData : public CohesiveZoneBaseStateData
{

public:

  typedef CohesiveZoneBaseStateData base;
  realT maxTraction;
  R1Tensor maxGap;


  InitiallyRigidCohesiveZoneStateData():
    base(),
    maxTraction(0),
    maxGap(0)
  {}

  InitiallyRigidCohesiveZoneStateData( const InitiallyRigidCohesiveZoneStateData& source):
    base( source ),
    maxTraction(source.maxTraction),
    maxGap(source.maxGap)
  {}

  ~InitiallyRigidCohesiveZoneStateData() {}
  friend class ConstitutiveBase;
  friend class InitiallyRigidCohesiveZone;

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
    R1TensorVarCounts += 1;

  }

  static void GetVariableNames( array<string>& intNames,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>& R2TensorNames,
                                array<string>& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("maxtraction");
    R1TensorNames.push_back("maxgap");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                   std::map<std::string, size_t>& realOffsets,
                                   std::map<std::string, size_t>& R1TensorOffsets,
                                   std::map<std::string, size_t>& R2TensorOffsets,
                                   std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["maxtraction"] = (char*)(&maxTraction) - (char*)this;
    R1TensorOffsets["maxgap"] = (char*)(&maxGap) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                  std::map<std::string, realT>& realValues,
                                  std::map<std::string, R1Tensor>& R1TensorValues,
                                  std::map<std::string, R2Tensor>& R2TensorValues,
                                  std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["maxtraction"] = maxTraction;
    R1TensorValues["maxgap"] = maxGap;
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
    (*(realVars[realVarCounts]))[elemNum] = maxTraction; realVarCounts += stride;
    (*(R1Vars[R1TensorVarCounts]))[elemNum] = maxGap; R1TensorVarCounts += stride;
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
    maxTraction = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    maxGap = (*(R1Vars[R1TensorVarCounts]))[elemNum]; R1TensorVarCounts += stride;
  }
  inline InitiallyRigidCohesiveZoneStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    maxTraction *= factor;
    maxGap *= factor;
    return *this;
  }

  inline InitiallyRigidCohesiveZoneStateData&
  operator=(const InitiallyRigidCohesiveZoneStateData& datum)
  {
    base::operator=(datum);
    maxTraction = datum.maxTraction;
    maxGap = datum.maxGap;
    return *this;
  }

  inline InitiallyRigidCohesiveZoneStateData&
  operator+=(const InitiallyRigidCohesiveZoneStateData& datum)
  {
    base::operator+=(datum);
    maxTraction += datum.maxTraction;
    maxGap += datum.maxGap;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, InitiallyRigidCohesiveZoneStateData& p0, InitiallyRigidCohesiveZoneStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, maxTraction, p0.maxTraction, p1.maxTraction);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, maxGap, p0.maxGap, p1.maxGap);

  }

  void MapFromRegion(const InitiallyRigidCohesiveZoneStateData& p0, const InitiallyRigidCohesiveZoneStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.maxTraction, p1.maxTraction, fct0, fct1, maxTraction);
    GeometryUtilities::MapFromRegion(p0.maxGap, p1.maxGap, fct0, fct1, maxGap);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************
class InitiallyRigidCohesiveZone : public CohesiveZoneBase
{
public:

  typedef InitiallyRigidCohesiveZoneParameterData ParameterClass;
  typedef InitiallyRigidCohesiveZoneStateData     StateClass;

  typedef array<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return "InitiallyRigidCohesiveZone"; }

  InitiallyRigidCohesiveZone();
  virtual ~InitiallyRigidCohesiveZone();

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
  { SetVariableParametersFromDerived<InitiallyRigidCohesiveZone>(varParams, newSize); }

  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<InitiallyRigidCohesiveZone>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<InitiallyRigidCohesiveZone>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<InitiallyRigidCohesiveZone>( num0 );
  }

  virtual void insert( const localIndex num )
  { InsertFromDerived<InitiallyRigidCohesiveZone>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<InitiallyRigidCohesiveZone>( num ); }

  void GetVariableNames( array<string>& intVars, array<string>& realVars, array<string>& R1TensorVars, array<string>& R2TensorVars,
                         array<string>& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<InitiallyRigidCohesiveZone>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<InitiallyRigidCohesiveZone>(name, type); }

  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<InitiallyRigidCohesiveZone>(name, type ); }

  bool GetStateValues( const std::string& name, array<real64>& values ) const
  { return GetStateValuesFromDerived<InitiallyRigidCohesiveZone>(name, values); }

  bool GetParameterValues( const std::string& name, array<real64>& values ) const
  { return GetParameterValuesFromDerived<InitiallyRigidCohesiveZone>(name, values); }

  bool SetStateValues( const std::string& name, const array<real64>& values )
  { return SetStateValuesFromDerived<InitiallyRigidCohesiveZone>(name, values); }

  bool SetParameterValues( const std::string& name, const array<real64>& values )
  { return SetParameterValuesFromDerived<InitiallyRigidCohesiveZone>(name, values); }

  virtual void Serialize( array<array<integer>*>& intVars, array<array<real64>*>& realVars, array<array<R1Tensor>*>& R1Vars, array<array<R2Tensor>*>& R2Vars,
                          array<array<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<InitiallyRigidCohesiveZone>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const array<array<integer>*>& intVars, const array<array<real64>*>& realVars, const array<array<R1Tensor>*>& R1Vars,
                            const array<array<R2Tensor>*>& R2Vars, const array<array<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<InitiallyRigidCohesiveZone>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<InitiallyRigidCohesiveZone>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<InitiallyRigidCohesiveZone>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<InitiallyRigidCohesiveZone>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<InitiallyRigidCohesiveZone>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<InitiallyRigidCohesiveZone>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<InitiallyRigidCohesiveZone>( localIndices, buffer ); }

  const StateClass* StateData( const localIndex index0, const localIndex index1 ) const
  { return &( m_stateData(index0,index1) );  }
  StateClass* StateData( const localIndex index0, const localIndex index1 )
  { return &( m_stateData(index0,index1) );  }

  const ParameterClass* ParameterData( const localIndex index ) const
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  ParameterClass* ParameterData( const localIndex index )
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }

  virtual int
  UpdateCohesiveZone( const localIndex index,
                      const R1Tensor& gap,
                      const R1Tensor& N,
                      const std::pair< ElementRegionT*, localIndex >& elem0,
                      const std::pair< ElementRegionT*, localIndex >& elem1,
                      R1Tensor& traction,
                      R2Tensor& stiffness );



private:
  InitiallyRigidCohesiveZone(const InitiallyRigidCohesiveZone&);
  InitiallyRigidCohesiveZone& operator=(const InitiallyRigidCohesiveZone&);


};
#endif /* INITIALLYRIGIDCOHESIVEZONE_H_ */
