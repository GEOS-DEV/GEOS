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

#ifndef LINEARIZED_H_
#define LINEARIZED_H_

#include "Utilities/GeometryUtilities.h"
#include "HertzianIntermediate.h"

/*
 * Linearized.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearizedParameterData : public HertzianIntermediateParameterData
{

public:

  typedef HertzianIntermediateParameterData base;
  realT interfaceEnergy;


  LinearizedParameterData():
    base(),
    interfaceEnergy(0)
  {}

  LinearizedParameterData( const LinearizedParameterData& source):
    base( source ),
    interfaceEnergy(source.interfaceEnergy)
  {}

  ~LinearizedParameterData() {}
  friend class ConstitutiveBase;
  friend class Linearized;

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
    realNames.push_back("interfaceEnergy");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["interfaceEnergy"] = (char*)(&interfaceEnergy) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["interfaceEnergy"] = interfaceEnergy;
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
    (*(realVars[realVarCounts]))[index] = interfaceEnergy; realVarCounts++;
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
    interfaceEnergy = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline LinearizedParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    interfaceEnergy *= factor;
    return *this;
  }

  inline LinearizedParameterData&
  operator=(const LinearizedParameterData& datum)
  {
    base::operator=(datum);
    interfaceEnergy = datum.interfaceEnergy;
    return *this;
  }

  inline LinearizedParameterData&
  operator+=(const LinearizedParameterData& datum)
  {
    base::operator+=(datum);
    interfaceEnergy += datum.interfaceEnergy;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, LinearizedParameterData& p0, LinearizedParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, interfaceEnergy, p0.interfaceEnergy, p1.interfaceEnergy);

  }

  void MapFromRegion(const LinearizedParameterData& p0, const LinearizedParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.interfaceEnergy, p1.interfaceEnergy, fct0, fct1, interfaceEnergy);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    HertzianIntermediateParameterData::ReadXML( node );
    interfaceEnergy = node.GetAttributeOrDefault("interfaceEnergy", 0.0);
    PostReadXML( node );

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearizedStateData : public HertzianIntermediateStateData
{

public:

  typedef HertzianIntermediateStateData base;
  int isCohesive;
  int loading;
  realT a0;
  realT at;
  realT ac;
  realT dcoh0;
  realT dupre;
  realT pullOffForce;
  realT Tk;
  realT coefficientOfRestitution2;
  realT yield;
  realT halfRestVel;
  realT ecrit;
  realT fcrit;
  realT linearStiffnessElastic;
  realT vark2;


  LinearizedStateData():
    base(),
    isCohesive(0),
    loading(0),
    a0(0),
    at(0),
    ac(0),
    dcoh0(0),
    dupre(0),
    pullOffForce(0),
    Tk(0),
    coefficientOfRestitution2(0),
    yield(0),
    halfRestVel(0),
    ecrit(0),
    fcrit(0),
    linearStiffnessElastic(0),
    vark2(0)
  {}

  LinearizedStateData( const LinearizedStateData& source):
    base( source ),
    isCohesive(source.isCohesive),
    loading(source.loading),
    a0(source.a0),
    at(source.at),
    ac(source.ac),
    dcoh0(source.dcoh0),
    dupre(source.dupre),
    pullOffForce(source.pullOffForce),
    Tk(source.Tk),
    coefficientOfRestitution2(source.coefficientOfRestitution2),
    yield(source.yield),
    halfRestVel(source.halfRestVel),
    ecrit(source.ecrit),
    fcrit(source.fcrit),
    linearStiffnessElastic(source.linearStiffnessElastic),
    vark2(source.vark2)
  {}

  ~LinearizedStateData() {}
  friend class ConstitutiveBase;
  friend class Linearized;

  virtual void
  Update(const realT curvature1,
                               const realT curvature2);

  virtual void
  Initialize(const realT curvature1, const realT curvature2,
                                   const realT poissons1, const realT poissons2,
                                   const realT youngs1, const realT youngs2,
                                   const realT mass1, const realT mass2,
                                   const realT rest1, const realT rest2,
                                   const realT yield1, const realT yield2,
                                   const realT velHalf1, const realT velHalf2,
                                   const realT surfaceEnergy1, const realT surfaceEnergy2,
                                   const realT cement1, const realT cement2);

  realT
  PullOffForce(const realT cement1,
                                    const realT cement2) const;

  realT
  EffectiveCoefficientOfRestitutionSquared(const realT Eeff,
                                                                const realT rest1, const realT rest2,
                                                                const realT poissons1, const realT poissons2,
                                                                const realT youngs1, const realT youngs2) const;

  realT
  EffectiveYieldStrength(const realT yield1,
                                              const realT yield2) const;

  realT
  HalfCoefficientOfRestitutionVelocity(const realT Eeff,
                                                            const realT velHalf1, const realT velHalf2,
                                                            const realT poissons1, const realT poissons2,
                                                            const realT youngs1, const realT youngs2) const;

  realT
  min_dt(const realT stiffness,
                              const realT tmass) const;

  realT
  ThetaK(const realT Tk_old,
                              const R1Tensor& shearDirection,
                              const R1Tensor& tangentialForces,
                              const realT mu_,
                              const realT fn_mag,
                              const realT dfn_mag,
                              int& loading_) const;

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
    intVarCounts += 2;
    realVarCounts += 14;

  }

  static void GetVariableNames( array<string>& intNames,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>& R2TensorNames,
                                array<string>& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    intNames.push_back("isCohesive");
    intNames.push_back("loading");
    realNames.push_back("a0");
    realNames.push_back("at");
    realNames.push_back("ac");
    realNames.push_back("cohesionDistance");
    realNames.push_back("dupre");
    realNames.push_back("pullOffForce");
    realNames.push_back("Tk");
    realNames.push_back("coefficientOfRestitution2");
    realNames.push_back("yield");
    realNames.push_back("halfRestVel");
    realNames.push_back("ecrit");
    realNames.push_back("fcrit");
    realNames.push_back("linearStiffnessElastic");
    realNames.push_back("vark2");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    intOffsets["isCohesive"] = (char*)(&isCohesive) - (char*)this;
    intOffsets["loading"] = (char*)(&loading) - (char*)this;
    realOffsets["a0"] = (char*)(&a0) - (char*)this;
    realOffsets["at"] = (char*)(&at) - (char*)this;
    realOffsets["ac"] = (char*)(&ac) - (char*)this;
    realOffsets["cohesionDistance"] = (char*)(&dcoh0) - (char*)this;
    realOffsets["dupre"] = (char*)(&dupre) - (char*)this;
    realOffsets["pullOffForce"] = (char*)(&pullOffForce) - (char*)this;
    realOffsets["Tk"] = (char*)(&Tk) - (char*)this;
    realOffsets["coefficientOfRestitution2"] = (char*)(&coefficientOfRestitution2) - (char*)this;
    realOffsets["yield"] = (char*)(&yield) - (char*)this;
    realOffsets["halfRestVel"] = (char*)(&halfRestVel) - (char*)this;
    realOffsets["ecrit"] = (char*)(&ecrit) - (char*)this;
    realOffsets["fcrit"] = (char*)(&fcrit) - (char*)this;
    realOffsets["linearStiffnessElastic"] = (char*)(&linearStiffnessElastic) - (char*)this;
    realOffsets["vark2"] = (char*)(&vark2) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    intValues["isCohesive"] = isCohesive;
    intValues["loading"] = loading;
    realValues["a0"] = a0;
    realValues["at"] = at;
    realValues["ac"] = ac;
    realValues["cohesionDistance"] = dcoh0;
    realValues["dupre"] = dupre;
    realValues["pullOffForce"] = pullOffForce;
    realValues["Tk"] = Tk;
    realValues["coefficientOfRestitution2"] = coefficientOfRestitution2;
    realValues["yield"] = yield;
    realValues["halfRestVel"] = halfRestVel;
    realValues["ecrit"] = ecrit;
    realValues["fcrit"] = fcrit;
    realValues["linearStiffnessElastic"] = linearStiffnessElastic;
    realValues["vark2"] = vark2;
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
    (*(intVars[intVarCounts]))[elemNum] = isCohesive; intVarCounts += stride;
    (*(intVars[intVarCounts]))[elemNum] = loading; intVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = a0; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = at; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = ac; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dcoh0; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dupre; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = pullOffForce; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = Tk; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = coefficientOfRestitution2; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = yield; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = halfRestVel; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = ecrit; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = fcrit; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = linearStiffnessElastic; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = vark2; realVarCounts += stride;
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
    isCohesive = (*(intVars[intVarCounts]))[elemNum]; intVarCounts += stride;
    loading = (*(intVars[intVarCounts]))[elemNum]; intVarCounts += stride;
    a0 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    at = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    ac = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dcoh0 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dupre = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    pullOffForce = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    Tk = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    coefficientOfRestitution2 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    yield = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    halfRestVel = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    ecrit = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    fcrit = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    linearStiffnessElastic = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    vark2 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
  }
  inline LinearizedStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    isCohesive *= factor;
    loading *= factor;
    a0 *= factor;
    at *= factor;
    ac *= factor;
    dcoh0 *= factor;
    dupre *= factor;
    pullOffForce *= factor;
    Tk *= factor;
    coefficientOfRestitution2 *= factor;
    yield *= factor;
    halfRestVel *= factor;
    ecrit *= factor;
    fcrit *= factor;
    linearStiffnessElastic *= factor;
    vark2 *= factor;
    return *this;
  }

  inline LinearizedStateData&
  operator=(const LinearizedStateData& datum)
  {
    base::operator=(datum);
    isCohesive = datum.isCohesive;
    loading = datum.loading;
    a0 = datum.a0;
    at = datum.at;
    ac = datum.ac;
    dcoh0 = datum.dcoh0;
    dupre = datum.dupre;
    pullOffForce = datum.pullOffForce;
    Tk = datum.Tk;
    coefficientOfRestitution2 = datum.coefficientOfRestitution2;
    yield = datum.yield;
    halfRestVel = datum.halfRestVel;
    ecrit = datum.ecrit;
    fcrit = datum.fcrit;
    linearStiffnessElastic = datum.linearStiffnessElastic;
    vark2 = datum.vark2;
    return *this;
  }

  inline LinearizedStateData&
  operator+=(const LinearizedStateData& datum)
  {
    base::operator+=(datum);
    isCohesive += datum.isCohesive;
    loading += datum.loading;
    a0 += datum.a0;
    at += datum.at;
    ac += datum.ac;
    dcoh0 += datum.dcoh0;
    dupre += datum.dupre;
    pullOffForce += datum.pullOffForce;
    Tk += datum.Tk;
    coefficientOfRestitution2 += datum.coefficientOfRestitution2;
    yield += datum.yield;
    halfRestVel += datum.halfRestVel;
    ecrit += datum.ecrit;
    fcrit += datum.fcrit;
    linearStiffnessElastic += datum.linearStiffnessElastic;
    vark2 += datum.vark2;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, LinearizedStateData& p0, LinearizedStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, isCohesive, p0.isCohesive, p1.isCohesive);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, loading, p0.loading, p1.loading);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, a0, p0.a0, p1.a0);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, at, p0.at, p1.at);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ac, p0.ac, p1.ac);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dcoh0, p0.dcoh0, p1.dcoh0);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dupre, p0.dupre, p1.dupre);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, pullOffForce, p0.pullOffForce, p1.pullOffForce);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, Tk, p0.Tk, p1.Tk);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, coefficientOfRestitution2, p0.coefficientOfRestitution2, p1.coefficientOfRestitution2);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, yield, p0.yield, p1.yield);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, halfRestVel, p0.halfRestVel, p1.halfRestVel);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ecrit, p0.ecrit, p1.ecrit);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, fcrit, p0.fcrit, p1.fcrit);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, linearStiffnessElastic, p0.linearStiffnessElastic, p1.linearStiffnessElastic);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, vark2, p0.vark2, p1.vark2);

  }

  void MapFromRegion(const LinearizedStateData& p0, const LinearizedStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.isCohesive, p1.isCohesive, fct0, fct1, isCohesive);
    GeometryUtilities::MapFromRegion(p0.loading, p1.loading, fct0, fct1, loading);
    GeometryUtilities::MapFromRegion(p0.a0, p1.a0, fct0, fct1, a0);
    GeometryUtilities::MapFromRegion(p0.at, p1.at, fct0, fct1, at);
    GeometryUtilities::MapFromRegion(p0.ac, p1.ac, fct0, fct1, ac);
    GeometryUtilities::MapFromRegion(p0.dcoh0, p1.dcoh0, fct0, fct1, dcoh0);
    GeometryUtilities::MapFromRegion(p0.dupre, p1.dupre, fct0, fct1, dupre);
    GeometryUtilities::MapFromRegion(p0.pullOffForce, p1.pullOffForce, fct0, fct1, pullOffForce);
    GeometryUtilities::MapFromRegion(p0.Tk, p1.Tk, fct0, fct1, Tk);
    GeometryUtilities::MapFromRegion(p0.coefficientOfRestitution2, p1.coefficientOfRestitution2, fct0, fct1, coefficientOfRestitution2);
    GeometryUtilities::MapFromRegion(p0.yield, p1.yield, fct0, fct1, yield);
    GeometryUtilities::MapFromRegion(p0.halfRestVel, p1.halfRestVel, fct0, fct1, halfRestVel);
    GeometryUtilities::MapFromRegion(p0.ecrit, p1.ecrit, fct0, fct1, ecrit);
    GeometryUtilities::MapFromRegion(p0.fcrit, p1.fcrit, fct0, fct1, fcrit);
    GeometryUtilities::MapFromRegion(p0.linearStiffnessElastic, p1.linearStiffnessElastic, fct0, fct1, linearStiffnessElastic);
    GeometryUtilities::MapFromRegion(p0.vark2, p1.vark2, fct0, fct1, vark2);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************
class Linearized: public HertzianIntermediate
{
public:

  typedef LinearizedParameterData ParameterClass;
  typedef LinearizedStateData     StateClass;

  typedef array<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return "Linearized"; }

  Linearized();
  virtual ~Linearized();

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
  { SetVariableParametersFromDerived<Linearized>(varParams, newSize); }

  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<Linearized>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<Linearized>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<Linearized>( num0 );
  }
  
  virtual void insert( const localIndex num )
  { InsertFromDerived<Linearized>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<Linearized>( num ); }
 
  void GetVariableNames( array<string>& intVars, array<string>& realVars, array<string>& R1TensorVars, array<string>& R2TensorVars, array<string>& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<Linearized>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<Linearized>(name, type); }
  
  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<Linearized>(name, type ); }
  
  bool GetStateValues( const std::string& name, array<real64>& values ) const
  { return GetStateValuesFromDerived<Linearized>(name, values); }

  bool GetParameterValues( const std::string& name, array<real64>& values ) const
  { return GetParameterValuesFromDerived<Linearized>(name, values); }

  bool SetStateValues( const std::string& name, const array<real64>& values )
  { return SetStateValuesFromDerived<Linearized>(name, values); }

  bool SetParameterValues( const std::string& name, const array<real64>& values )
  { return SetParameterValuesFromDerived<Linearized>(name, values); }

  virtual void Serialize( array<array<integer>*>& intVars, array<array<real64>*>& realVars, array<array<R1Tensor>*>& R1Vars, array<array<R2Tensor>*>& R2Vars, array<array<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<Linearized>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const array<array<integer>*>& intVars, const array<array<real64>*>& realVars, const array<array<R1Tensor>*>& R1Vars, const array<array<R2Tensor>*>& R2Vars, const array<array<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<Linearized>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<Linearized>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<Linearized>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<Linearized>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<Linearized>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<Linearized>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<Linearized>( localIndices, buffer ); }

  const StateClass* StateData( const localIndex index0, const localIndex index1 ) const
  { return &( m_stateData(index0,index1) );  }
  StateClass* StateData( const localIndex index0, const localIndex index1 )
  { return &( m_stateData(index0,index1) );  }

  const ParameterClass* ParameterData( const localIndex index ) const
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  ParameterClass* ParameterData( const localIndex index )
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  

protected:
  virtual realT
  NormalStiffness(const InterfaceBaseParameterData& ,
                              InterfaceBaseStateData& matStateBase,
                              const realT normalApproach,
                              const bool setForces) const;

  realT
  PlasticUnloadingSlope(const realT Eeff, const realT Eeff2, const realT tmass,
                                    const int vark2, const realT linearStiffnessElastic,
                                    const realT at, const realT ac, const realT velHalf) const;


private:
  Linearized(const Linearized&);
  Linearized& operator=(const Linearized&);
  

};
#endif /* LINEARIZED_H_ */
