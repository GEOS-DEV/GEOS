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

#ifndef HERTZIANINTERMEDIATE_H_
#define HERTZIANINTERMEDIATE_H_

#include "Utilities/GeometryUtilities.h"
#include "InterfaceBase.h"

/*
 * HertzianIntermediate.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */



//**********************************************************************************************************************
//**********************************************************************************************************************


class HertzianIntermediateParameterData : public InterfaceBaseParameterData
{

public:

  typedef InterfaceBaseParameterData base;


  HertzianIntermediateParameterData():
    base()
  {}

  HertzianIntermediateParameterData( const HertzianIntermediateParameterData& source):
    base( source )
  {}

  ~HertzianIntermediateParameterData() {}
  friend class ConstitutiveBase;
  friend class HertzianIntermediate;

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
  inline HertzianIntermediateParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    return *this;
  }

  inline HertzianIntermediateParameterData&
  operator=(const HertzianIntermediateParameterData& datum)
  {
    base::operator=(datum);
    return *this;
  }

  inline HertzianIntermediateParameterData&
  operator+=(const HertzianIntermediateParameterData& datum)
  {
    base::operator+=(datum);
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, HertzianIntermediateParameterData& p0, HertzianIntermediateParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);

  }

  void MapFromRegion(const HertzianIntermediateParameterData& p0, const HertzianIntermediateParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    InterfaceBaseParameterData::ReadXML( node );

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class HertzianIntermediateStateData : public InterfaceBaseStateData
{

public:

  typedef InterfaceBaseStateData base;
  realT mu;
  realT dissipatedViscousEnergy;
  realT radius;
  realT youngs;
  realT mass;
  realT hertzCf;


  HertzianIntermediateStateData():
    base(),
    mu(0),
    dissipatedViscousEnergy(0),
    radius(0),
    youngs(0),
    mass(0),
    hertzCf(0)
  {}

  HertzianIntermediateStateData( const HertzianIntermediateStateData& source):
    base( source ),
    mu(source.mu),
    dissipatedViscousEnergy(source.dissipatedViscousEnergy),
    radius(source.radius),
    youngs(source.youngs),
    mass(source.mass),
    hertzCf(source.hertzCf)
  {}

  ~HertzianIntermediateStateData() {}
  friend class ConstitutiveBase;
  friend class HertzianIntermediate;

  virtual void
  Update(const realT curvature1,
         const realT curvature2);

  virtual void
  Initialize(const realT curvature1, const realT curvature2,
             const realT poissons1, const realT poissons2,
             const realT youngs1, const realT youngs2,
             const realT mass1, const realT mass2,
             const realT, const realT,
             const realT, const realT,
             const realT, const realT,
             const realT, const realT,
             const realT, const realT );

  realT
  EffectiveRadius(const realT curvature1,
                  const realT curvature2) const;

  realT
  EffectiveYoungsModulus(const realT poissons1,
                         const realT poissons2,
                         const realT youngs1,
                         const realT youngs2) const;

  realT
  EffectiveMass(const realT mass1,
                const realT mass2) const;

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
    realVarCounts += 6;

  }

  static void GetVariableNames( array<string>& intNames,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>& R2TensorNames,
                                array<string>& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("currentFrictionCoefficient");
    realNames.push_back("dissipatedViscousEnergy");
    realNames.push_back("radius");
    realNames.push_back("youngsModulus");
    realNames.push_back("mass");
    realNames.push_back("hertzianCoefficient");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                   std::map<std::string, size_t>& realOffsets,
                                   std::map<std::string, size_t>& R1TensorOffsets,
                                   std::map<std::string, size_t>& R2TensorOffsets,
                                   std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["currentFrictionCoefficient"] = (char*)(&mu) - (char*)this;
    realOffsets["dissipatedViscousEnergy"] = (char*)(&dissipatedViscousEnergy) - (char*)this;
    realOffsets["radius"] = (char*)(&radius) - (char*)this;
    realOffsets["youngsModulus"] = (char*)(&youngs) - (char*)this;
    realOffsets["mass"] = (char*)(&mass) - (char*)this;
    realOffsets["hertzianCoefficient"] = (char*)(&hertzCf) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                  std::map<std::string, realT>& realValues,
                                  std::map<std::string, R1Tensor>& R1TensorValues,
                                  std::map<std::string, R2Tensor>& R2TensorValues,
                                  std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["currentFrictionCoefficient"] = mu;
    realValues["dissipatedViscousEnergy"] = dissipatedViscousEnergy;
    realValues["radius"] = radius;
    realValues["youngsModulus"] = youngs;
    realValues["mass"] = mass;
    realValues["hertzianCoefficient"] = hertzCf;
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
    (*(realVars[realVarCounts]))[elemNum] = mu; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dissipatedViscousEnergy; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = radius; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = youngs; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = mass; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = hertzCf; realVarCounts += stride;
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
    mu = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dissipatedViscousEnergy = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    radius = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    youngs = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    mass = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    hertzCf = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
  }
  inline HertzianIntermediateStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    mu *= factor;
    dissipatedViscousEnergy *= factor;
    radius *= factor;
    youngs *= factor;
    mass *= factor;
    hertzCf *= factor;
    return *this;
  }

  inline HertzianIntermediateStateData&
  operator=(const HertzianIntermediateStateData& datum)
  {
    base::operator=(datum);
    mu = datum.mu;
    dissipatedViscousEnergy = datum.dissipatedViscousEnergy;
    radius = datum.radius;
    youngs = datum.youngs;
    mass = datum.mass;
    hertzCf = datum.hertzCf;
    return *this;
  }

  inline HertzianIntermediateStateData&
  operator+=(const HertzianIntermediateStateData& datum)
  {
    base::operator+=(datum);
    mu += datum.mu;
    dissipatedViscousEnergy += datum.dissipatedViscousEnergy;
    radius += datum.radius;
    youngs += datum.youngs;
    mass += datum.mass;
    hertzCf += datum.hertzCf;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, HertzianIntermediateStateData& p0, HertzianIntermediateStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mu, p0.mu, p1.mu);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dissipatedViscousEnergy, p0.dissipatedViscousEnergy, p1.dissipatedViscousEnergy);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, radius, p0.radius, p1.radius);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, youngs, p0.youngs, p1.youngs);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mass, p0.mass, p1.mass);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, hertzCf, p0.hertzCf, p1.hertzCf);

  }

  void MapFromRegion(const HertzianIntermediateStateData& p0, const HertzianIntermediateStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.mu, p1.mu, fct0, fct1, mu);
    GeometryUtilities::MapFromRegion(p0.dissipatedViscousEnergy, p1.dissipatedViscousEnergy, fct0, fct1, dissipatedViscousEnergy);
    GeometryUtilities::MapFromRegion(p0.radius, p1.radius, fct0, fct1, radius);
    GeometryUtilities::MapFromRegion(p0.youngs, p1.youngs, fct0, fct1, youngs);
    GeometryUtilities::MapFromRegion(p0.mass, p1.mass, fct0, fct1, mass);
    GeometryUtilities::MapFromRegion(p0.hertzCf, p1.hertzCf, fct0, fct1, hertzCf);

  }



};


//**********************************************************************************************************************
//**********************************************************************************************************************


class HertzianIntermediate : public InterfaceBase
{
public:

  typedef HertzianIntermediateParameterData ParameterClass;
  typedef HertzianIntermediateStateData     StateClass;


  HertzianIntermediate( const int paramSize, const int stateSize );

  virtual ~HertzianIntermediate();

  virtual void ReadXML( TICPP::HierarchicalDataNode& node ) = 0;

  virtual void resize( const localIndex num ) = 0;

  virtual void resize( const localIndex num0,
                       const localIndex num1 ) = 0;

  virtual void insert( const localIndex num ) = 0;

  virtual void erase( const localIndex num ) = 0;

  virtual const HertzianIntermediateStateData* StateData( const localIndex index0,
                                                          const localIndex index1 ) const = 0;
  virtual       HertzianIntermediateStateData* StateData( const localIndex index0,
                                                          const localIndex index1 )  = 0;

  virtual const HertzianIntermediateParameterData* ParameterData( const localIndex index ) const = 0;
  virtual       HertzianIntermediateParameterData* ParameterData( const localIndex index ) = 0;

  inline void IncrementPtr( const HertzianIntermediateStateData* ptr ) const
  {
    ptr = reinterpret_cast<const HertzianIntermediateStateData*>( reinterpret_cast<const char*>(ptr) + m_stateSize );
  }

  inline void IncrementPtr( const HertzianIntermediateParameterData* ptr ) const
  {
    ptr = reinterpret_cast<const HertzianIntermediateParameterData*>( reinterpret_cast<const char*>(ptr) + m_paramSize );
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
  UpdateProperties(const localIndex index, std::map<std::string, realT>& p0, std::map<std::string, realT>& p1);


protected:
  virtual realT
  NormalStiffness(const InterfaceBaseParameterData&,
                  InterfaceBaseStateData& matStateBase,
                  const realT normalApproach,
                  const bool setForces) const;

  virtual void
  UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                  InterfaceBaseStateData& matStateBase) const;

  virtual realT
  ShearStrength(const InterfaceBaseParameterData&,
                InterfaceBaseStateData& matStateBase) const;


private:
  HertzianIntermediate();
  HertzianIntermediate( const HertzianIntermediate& );
  HertzianIntermediate& operator=( const HertzianIntermediate& );


};
#endif /* HERTZIANINTERMEDIATE_H_ */
