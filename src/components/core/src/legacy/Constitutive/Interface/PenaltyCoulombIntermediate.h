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

#ifndef PENALTYCOULOMBINTERMEDIATE_H_
#define PENALTYCOULOMBINTERMEDIATE_H_

#include "Utilities/GeometryUtilities.h"
#include "InterfaceBase.h"

/*
 * PenaltyCoulombIntermediate.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class PenaltyCoulombIntermediateParameterData : public InterfaceBaseParameterData
{

public:

  typedef InterfaceBaseParameterData base;
  realT aperture;
  realT normalApproachYield;
  realT ksoften;
  realT stressYield;
  realT stressSoften;
  realT kshear;
  realT ktildeAperture;
  realT kyield;
  realT normalApproachSoften;


  PenaltyCoulombIntermediateParameterData():
    base(),
    aperture(0),
    normalApproachYield(0),
    ksoften(0),
    stressYield(0),
    stressSoften(0),
    kshear(0),
    ktildeAperture(0),
    kyield(0),
    normalApproachSoften(0)
  {}

  PenaltyCoulombIntermediateParameterData( const PenaltyCoulombIntermediateParameterData& source):
    base( source ),
    aperture(source.aperture),
    normalApproachYield(source.normalApproachYield),
    ksoften(source.ksoften),
    stressYield(source.stressYield),
    stressSoften(source.stressSoften),
    kshear(source.kshear),
    ktildeAperture(source.ktildeAperture),
    kyield(source.kyield),
    normalApproachSoften(source.normalApproachSoften)
  {}

  ~PenaltyCoulombIntermediateParameterData() {}
  friend class ConstitutiveBase;
  friend class PenaltyCoulombIntermediate;

  realT
  Stiffness(const realT normalApproach) const;

  realT
  Stress(const realT normalApproach) const;

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
    realVarCounts += 9;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("aperture");
    realNames.push_back("normalApproachYield");
    realNames.push_back("arealStiffnessSoften");
    realNames.push_back("stressYield");
    realNames.push_back("stressSoften");
    realNames.push_back("arealStiffnessShear");
    realNames.push_back("ktildeAperture");
    realNames.push_back("arealStiffnessYield");
    realNames.push_back("normalApproachSoften");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["aperture"] = (char*)(&aperture) - (char*)this;
    realOffsets["normalApproachYield"] = (char*)(&normalApproachYield) - (char*)this;
    realOffsets["arealStiffnessSoften"] = (char*)(&ksoften) - (char*)this;
    realOffsets["stressYield"] = (char*)(&stressYield) - (char*)this;
    realOffsets["stressSoften"] = (char*)(&stressSoften) - (char*)this;
    realOffsets["arealStiffnessShear"] = (char*)(&kshear) - (char*)this;
    realOffsets["ktildeAperture"] = (char*)(&ktildeAperture) - (char*)this;
    realOffsets["arealStiffnessYield"] = (char*)(&kyield) - (char*)this;
    realOffsets["normalApproachSoften"] = (char*)(&normalApproachSoften) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["aperture"] = aperture;
    realValues["normalApproachYield"] = normalApproachYield;
    realValues["arealStiffnessSoften"] = ksoften;
    realValues["stressYield"] = stressYield;
    realValues["stressSoften"] = stressSoften;
    realValues["arealStiffnessShear"] = kshear;
    realValues["ktildeAperture"] = ktildeAperture;
    realValues["arealStiffnessYield"] = kyield;
    realValues["normalApproachSoften"] = normalApproachSoften;
  }

  void Serialize(const localIndex index,
                  Array1dT<iArray1d*>& intVars,
                  Array1dT<rArray1d*>& realVars,
                  Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                  Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                  Array1dT<Array1dT<R2SymTensor>*>& R2SymVars,
                  localIndex& intVarCounts,
                  localIndex& realVarCounts,
                  localIndex& R1TensorVarCounts,
                  localIndex& R2TensorVarCounts,
                  localIndex& R2SymTensorVarCounts  ) const
  {
    base::Serialize(index, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
    (*(realVars[realVarCounts]))[index] = aperture; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = normalApproachYield; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = ksoften; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = stressYield; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = stressSoften; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = kshear; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = ktildeAperture; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = kyield; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = normalApproachSoften; realVarCounts++;
  }


  void  Deserialize( const localIndex index,
                     const Array1dT<iArray1d*>& intVars,
                     const Array1dT<rArray1d*>& realVars,
                     const Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                     const Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                     const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars,
                  localIndex& intVarCounts,
                  localIndex& realVarCounts,
                  localIndex& R1TensorVarCounts,
                  localIndex& R2TensorVarCounts,
                  localIndex& R2SymTensorVarCounts )
  {
    base::Deserialize(index, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
    aperture = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    normalApproachYield = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    ksoften = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    stressYield = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    stressSoften = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    kshear = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    ktildeAperture = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    kyield = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    normalApproachSoften = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline PenaltyCoulombIntermediateParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    aperture *= factor;
    normalApproachYield *= factor;
    ksoften *= factor;
    stressYield *= factor;
    stressSoften *= factor;
    kshear *= factor;
    ktildeAperture *= factor;
    kyield *= factor;
    normalApproachSoften *= factor;
    return *this;
  }

  inline PenaltyCoulombIntermediateParameterData&
  operator=(const PenaltyCoulombIntermediateParameterData& datum)
  {
    base::operator=(datum);
    aperture = datum.aperture;
    normalApproachYield = datum.normalApproachYield;
    ksoften = datum.ksoften;
    stressYield = datum.stressYield;
    stressSoften = datum.stressSoften;
    kshear = datum.kshear;
    ktildeAperture = datum.ktildeAperture;
    kyield = datum.kyield;
    normalApproachSoften = datum.normalApproachSoften;
    return *this;
  }

  inline PenaltyCoulombIntermediateParameterData&
  operator+=(const PenaltyCoulombIntermediateParameterData& datum)
  {
    base::operator+=(datum);
    aperture += datum.aperture;
    normalApproachYield += datum.normalApproachYield;
    ksoften += datum.ksoften;
    stressYield += datum.stressYield;
    stressSoften += datum.stressSoften;
    kshear += datum.kshear;
    ktildeAperture += datum.ktildeAperture;
    kyield += datum.kyield;
    normalApproachSoften += datum.normalApproachSoften;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, PenaltyCoulombIntermediateParameterData& p0, PenaltyCoulombIntermediateParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, aperture, p0.aperture, p1.aperture);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, normalApproachYield, p0.normalApproachYield, p1.normalApproachYield);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ksoften, p0.ksoften, p1.ksoften);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressYield, p0.stressYield, p1.stressYield);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressSoften, p0.stressSoften, p1.stressSoften);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, kshear, p0.kshear, p1.kshear);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ktildeAperture, p0.ktildeAperture, p1.ktildeAperture);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, kyield, p0.kyield, p1.kyield);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, normalApproachSoften, p0.normalApproachSoften, p1.normalApproachSoften);

  }

  void MapFromRegion(const PenaltyCoulombIntermediateParameterData& p0, const PenaltyCoulombIntermediateParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.aperture, p1.aperture, fct0, fct1, aperture);
    GeometryUtilities::MapFromRegion(p0.normalApproachYield, p1.normalApproachYield, fct0, fct1, normalApproachYield);
    GeometryUtilities::MapFromRegion(p0.ksoften, p1.ksoften, fct0, fct1, ksoften);
    GeometryUtilities::MapFromRegion(p0.stressYield, p1.stressYield, fct0, fct1, stressYield);
    GeometryUtilities::MapFromRegion(p0.stressSoften, p1.stressSoften, fct0, fct1, stressSoften);
    GeometryUtilities::MapFromRegion(p0.kshear, p1.kshear, fct0, fct1, kshear);
    GeometryUtilities::MapFromRegion(p0.ktildeAperture, p1.ktildeAperture, fct0, fct1, ktildeAperture);
    GeometryUtilities::MapFromRegion(p0.kyield, p1.kyield, fct0, fct1, kyield);
    GeometryUtilities::MapFromRegion(p0.normalApproachSoften, p1.normalApproachSoften, fct0, fct1, normalApproachSoften);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    InterfaceBaseParameterData::ReadXML( node );
    aperture = node.GetAttributeOrDefault("aperture", 0.0);
    normalApproachYield = node.GetAttributeOrDefault("normalApproachYield", 0.0);
    ksoften = node.GetAttributeOrDefault("arealStiffnessSoften", 0.0);
    stressYield = node.GetAttributeOrDefault("stressYield", 0.0);
    stressSoften = node.GetAttributeOrDefault("stressSoften", 0.0);
    kshear = node.GetAttributeOrDefault("arealStiffnessShear", 0.0);
    ktildeAperture = node.GetAttributeOrDefault("ktildeAperture", 0.0);
    kyield = node.GetAttributeOrDefault("arealStiffnessYield", 0.0);
    normalApproachSoften = node.GetAttributeOrDefault("normalApproachSoften", 0.0);

  }

protected:
  void
  Initialize();

  virtual void
  PostReadXML( TICPP::HierarchicalDataNode& );

  virtual void
  PostSetValue();



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class PenaltyCoulombIntermediateStateData : public InterfaceBaseStateData
{

public:

  typedef InterfaceBaseStateData base;


  PenaltyCoulombIntermediateStateData():
    base()
  {}

  PenaltyCoulombIntermediateStateData( const PenaltyCoulombIntermediateStateData& source):
    base( source )
  {}

  ~PenaltyCoulombIntermediateStateData() {}
  friend class ConstitutiveBase;
  friend class PenaltyCoulombIntermediate;

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

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
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
                  Array1dT<iArray1d*>& intVars,
                  Array1dT<rArray1d*>& realVars,
                  Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                  Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                  Array1dT<Array1dT<R2SymTensor>*>& R2SymVars,
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
                     const Array1dT<iArray1d*>& intVars,
                     const Array1dT<rArray1d*>& realVars,
                     const Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                     const Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                     const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars,
                  localIndex& intVarCounts,
                  localIndex& realVarCounts,
                  localIndex& R1TensorVarCounts,
                  localIndex& R2TensorVarCounts,
                  localIndex& R2SymTensorVarCounts )
  {
    base::Deserialize(index, stride, elemNum, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
  }
  inline PenaltyCoulombIntermediateStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    return *this;
  }

  inline PenaltyCoulombIntermediateStateData&
  operator=(const PenaltyCoulombIntermediateStateData& datum)
  {
    base::operator=(datum);
    return *this;
  }

  inline PenaltyCoulombIntermediateStateData&
  operator+=(const PenaltyCoulombIntermediateStateData& datum)
  {
    base::operator+=(datum);
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, PenaltyCoulombIntermediateStateData& p0, PenaltyCoulombIntermediateStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);

  }

  void MapFromRegion(const PenaltyCoulombIntermediateStateData& p0, const PenaltyCoulombIntermediateStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);

  }



};


//**********************************************************************************************************************
//**********************************************************************************************************************


class PenaltyCoulombIntermediate: public InterfaceBase
{
public:
  
  typedef PenaltyCoulombIntermediateParameterData ParameterClass;
  typedef PenaltyCoulombIntermediateStateData     StateClass;
  

  PenaltyCoulombIntermediate( const int paramSize, const int stateSize );

  virtual ~PenaltyCoulombIntermediate();
  
  virtual void ReadXML( TICPP::HierarchicalDataNode& node ) = 0;

  virtual void resize( const localIndex num ) = 0;
  
  virtual void resize( const localIndex num0,
                       const localIndex num1 ) = 0;
  
  virtual void insert( const localIndex num ) = 0;

  virtual void erase( const localIndex num ) = 0;
  
  virtual const PenaltyCoulombIntermediateStateData* StateData( const localIndex index0,
                                                  const localIndex index1 ) const = 0;
  virtual       PenaltyCoulombIntermediateStateData* StateData( const localIndex index0,
                                                  const localIndex index1 )  = 0;

  virtual const PenaltyCoulombIntermediateParameterData* ParameterData( const localIndex index ) const = 0;
  virtual       PenaltyCoulombIntermediateParameterData* ParameterData( const localIndex index ) = 0;
  
  inline void IncrementPtr( const PenaltyCoulombIntermediateStateData* ptr ) const
  {
    ptr = reinterpret_cast<const PenaltyCoulombIntermediateStateData*>( reinterpret_cast<const char*>(ptr) + m_stateSize );
  }

  inline void IncrementPtr( const PenaltyCoulombIntermediateParameterData* ptr ) const
  {
    ptr = reinterpret_cast<const PenaltyCoulombIntermediateParameterData*>( reinterpret_cast<const char*>(ptr) + m_paramSize );
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

  virtual void
  StrainDrivenUpdate( const localIndex index );

  virtual realT
  StiffnessProjected(const localIndex index);


  
private:
  PenaltyCoulombIntermediate();
  PenaltyCoulombIntermediate( const PenaltyCoulombIntermediate& );
  PenaltyCoulombIntermediate& operator=( const PenaltyCoulombIntermediate& );
  
  
};
#endif /* PENALTYCOULOMBINTERMEDIATE_H_ */
