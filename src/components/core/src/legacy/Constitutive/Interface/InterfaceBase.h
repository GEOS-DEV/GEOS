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

#ifndef INTERFACEBASE_H_
#define INTERFACEBASE_H_

#include "Utilities/GeometryUtilities.h"
#include "Constitutive/ConstitutiveBase.h"

/*
 * InterfaceBase.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */



//**********************************************************************************************************************
//**********************************************************************************************************************


class InterfaceBaseParameterData
{

public:
  realT ston;
  realT mu0;
  realT kappa0;
  realT dkappadsFct;
  realT dkappadsExp;
  realT dkappadnFct;
  realT dkappadnExp;


  InterfaceBaseParameterData():
    ston(0),
    mu0(0),
    kappa0(0),
    dkappadsFct(0),
    dkappadsExp(0),
    dkappadnFct(0),
    dkappadnExp(0)
  {}

  InterfaceBaseParameterData( const InterfaceBaseParameterData& source):
    ston(source.ston),
    mu0(source.mu0),
    kappa0(source.kappa0),
    dkappadsFct(source.dkappadsFct),
    dkappadsExp(source.dkappadsExp),
    dkappadnFct(source.dkappadnFct),
    dkappadnExp(source.dkappadnExp)
  {}

  ~InterfaceBaseParameterData() {}

  static void GetVariableCounts( localIndex&,
                                 localIndex& realVarCounts,
                                 localIndex&,
                                 localIndex&,
                                 localIndex&  )
  {
    realVarCounts = 7;

  }

  static void GetVariableNames( array<string>&,
                                array<string>& realNames,
                                array<string>&,
                                array<string>&,
                                array<string>&  )
  {
    realNames.push_back("shearToNormalStiffnessRatio");
    realNames.push_back("frictionCoefficient");
    realNames.push_back("permeabilityInitial");
    realNames.push_back("permeabilityShearCoefficient");
    realNames.push_back("permeabilityShearExponent");
    realNames.push_back("permeabilityNormalCoefficient");
    realNames.push_back("permeabilityNormalExponent");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>&,
                                   std::map<std::string, size_t>& realOffsets,
                                   std::map<std::string, size_t>&,
                                   std::map<std::string, size_t>&,
                                   std::map<std::string, size_t>&  ) const
  {
    realOffsets["shearToNormalStiffnessRatio"] = (char*)(&ston) - (char*)this;
    realOffsets["frictionCoefficient"] = (char*)(&mu0) - (char*)this;
    realOffsets["permeabilityInitial"] = (char*)(&kappa0) - (char*)this;
    realOffsets["permeabilityShearCoefficient"] = (char*)(&dkappadsFct) - (char*)this;
    realOffsets["permeabilityShearExponent"] = (char*)(&dkappadsExp) - (char*)this;
    realOffsets["permeabilityNormalCoefficient"] = (char*)(&dkappadnFct) - (char*)this;
    realOffsets["permeabilityNormalExponent"] = (char*)(&dkappadnExp) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>&,
                                  std::map<std::string, realT>& realValues,
                                  std::map<std::string, R1Tensor>&,
                                  std::map<std::string, R2Tensor>&,
                                  std::map<std::string, R2SymTensor>&  )
  {
    realValues["shearToNormalStiffnessRatio"] = ston;
    realValues["frictionCoefficient"] = mu0;
    realValues["permeabilityInitial"] = kappa0;
    realValues["permeabilityShearCoefficient"] = dkappadsFct;
    realValues["permeabilityShearExponent"] = dkappadsExp;
    realValues["permeabilityNormalCoefficient"] = dkappadnFct;
    realValues["permeabilityNormalExponent"] = dkappadnExp;
  }

  void Serialize(const localIndex index,
                 array<array<integer>*>&,
                 array<array<real64>*>& realVars,
                 array<array<R1Tensor>*>&,
                 array<array<R2Tensor>*>&,
                 array<array<R2SymTensor>*>&,
                 localIndex&,
                 localIndex& realVarCounts,
                 localIndex&,
                 localIndex&,
                 localIndex&   ) const
  {
    (*(realVars[realVarCounts]))[index] = ston; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = mu0; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = kappa0; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = dkappadsFct; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = dkappadsExp; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = dkappadnFct; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = dkappadnExp; realVarCounts++;
  }


  void  Deserialize( const localIndex index,
                     const array<array<integer>*>&,
                     const array<array<real64>*>& realVars,
                     const array<array<R1Tensor>*>&,
                     const array<array<R2Tensor>*>&,
                     const array<array<R2SymTensor>*>&,
                     localIndex&,
                     localIndex& realVarCounts,
                     localIndex&,
                     localIndex&,
                     localIndex&   )
  {
    ston = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    mu0 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    kappa0 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    dkappadsFct = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    dkappadsExp = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    dkappadnFct = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    dkappadnExp = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline InterfaceBaseParameterData&
  operator*=(const realT factor)
  {
    ston *= factor;
    mu0 *= factor;
    kappa0 *= factor;
    dkappadsFct *= factor;
    dkappadsExp *= factor;
    dkappadnFct *= factor;
    dkappadnExp *= factor;
    return *this;
  }

  inline InterfaceBaseParameterData&
  operator=(const InterfaceBaseParameterData& datum)
  {
    ston = datum.ston;
    mu0 = datum.mu0;
    kappa0 = datum.kappa0;
    dkappadsFct = datum.dkappadsFct;
    dkappadsExp = datum.dkappadsExp;
    dkappadnFct = datum.dkappadnFct;
    dkappadnExp = datum.dkappadnExp;
    return *this;
  }

  inline InterfaceBaseParameterData&
  operator+=(const InterfaceBaseParameterData& datum)
  {
    ston += datum.ston;
    mu0 += datum.mu0;
    kappa0 += datum.kappa0;
    dkappadsFct += datum.dkappadsFct;
    dkappadsExp += datum.dkappadsExp;
    dkappadnFct += datum.dkappadnFct;
    dkappadnExp += datum.dkappadnExp;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, InterfaceBaseParameterData& p0, InterfaceBaseParameterData& p1)
  {
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ston, p0.ston, p1.ston);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mu0, p0.mu0, p1.mu0);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, kappa0, p0.kappa0, p1.kappa0);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dkappadsFct, p0.dkappadsFct, p1.dkappadsFct);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dkappadsExp, p0.dkappadsExp, p1.dkappadsExp);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dkappadnFct, p0.dkappadnFct, p1.dkappadnFct);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dkappadnExp, p0.dkappadnExp, p1.dkappadnExp);

  }

  void MapFromRegion(const InterfaceBaseParameterData& p0, const InterfaceBaseParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    GeometryUtilities::MapFromRegion(p0.ston, p1.ston, fct0, fct1, ston);
    GeometryUtilities::MapFromRegion(p0.mu0, p1.mu0, fct0, fct1, mu0);
    GeometryUtilities::MapFromRegion(p0.kappa0, p1.kappa0, fct0, fct1, kappa0);
    GeometryUtilities::MapFromRegion(p0.dkappadsFct, p1.dkappadsFct, fct0, fct1, dkappadsFct);
    GeometryUtilities::MapFromRegion(p0.dkappadsExp, p1.dkappadsExp, fct0, fct1, dkappadsExp);
    GeometryUtilities::MapFromRegion(p0.dkappadnFct, p1.dkappadnFct, fct0, fct1, dkappadnFct);
    GeometryUtilities::MapFromRegion(p0.dkappadnExp, p1.dkappadnExp, fct0, fct1, dkappadnExp);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    ston = node.GetAttributeOrDefault("shearToNormalStiffnessRatio", 0.0);
    mu0 = node.GetAttributeOrDefault("frictionCoefficient", 0.0);
    kappa0 = node.GetAttributeOrDefault("permeabilityInitial", 0.0);
    dkappadsFct = node.GetAttributeOrDefault("permeabilityShearCoefficient", 0.0);
    dkappadsExp = node.GetAttributeOrDefault("permeabilityShearExponent", 0.0);
    dkappadnFct = node.GetAttributeOrDefault("permeabilityNormalCoefficient", 0.0);
    dkappadnExp = node.GetAttributeOrDefault("permeabilityNormalExponent", 0.0);

  }
  virtual void PostReadXML( TICPP::HierarchicalDataNode& ) {}



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class InterfaceBaseStateData
{

public:
  realT ElasticStrainEnergy;
  realT DissipatedEnergy;
  realT stress;
  realT stressShear;
  realT dxndt;
  realT dxsdt;
  realT xs;
  realT normalApproach;
  realT dt;
  realT kappa;
  R1Tensor stressShearVector;


  InterfaceBaseStateData():
    ElasticStrainEnergy(0),
    DissipatedEnergy(0),
    stress(0),
    stressShear(0),
    dxndt(0),
    dxsdt(0),
    xs(0),
    normalApproach(0),
    dt(0),
    kappa(0),
    stressShearVector(0)
  {}

  InterfaceBaseStateData( const InterfaceBaseStateData& source):
    ElasticStrainEnergy(source.ElasticStrainEnergy),
    DissipatedEnergy(source.DissipatedEnergy),
    stress(source.stress),
    stressShear(source.stressShear),
    dxndt(source.dxndt),
    dxsdt(source.dxsdt),
    xs(source.xs),
    normalApproach(source.normalApproach),
    dt(source.dt),
    kappa(source.kappa),
    stressShearVector(source.stressShearVector)
  {}

  ~InterfaceBaseStateData() {}

  void
  UpdateOrientation(const realT normalApproachIn,
                    const realT dtIn,
                    const R1Tensor& normal,
                    const R1Tensor& velocity,
                    R1Tensor& dShearSlip,
                    R1Tensor& shearSlip);

  static void GetVariableCounts( localIndex&,
                                 localIndex& realVarCounts,
                                 localIndex& R1TensorVarCounts,
                                 localIndex&,
                                 localIndex&  )
  {
    realVarCounts = 10;
    R1TensorVarCounts = 1;

  }

  static void GetVariableNames( array<string>&,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>&,
                                array<string>&  )
  {
    realNames.push_back("elasticStrainEnergy");
    realNames.push_back("dissipatedEnergy");
    realNames.push_back("stress");
    realNames.push_back("shearStress");
    realNames.push_back("normalVelocity");
    realNames.push_back("shearVelocity");
    realNames.push_back("shear");
    realNames.push_back("normalApproach");
    realNames.push_back("timestep");
    realNames.push_back("permeability");
    R1TensorNames.push_back("shearStressVector");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>&,
                                   std::map<std::string, size_t>& realOffsets,
                                   std::map<std::string, size_t>& R1TensorOffsets,
                                   std::map<std::string, size_t>&,
                                   std::map<std::string, size_t>&  ) const
  {
    realOffsets["elasticStrainEnergy"] = (char*)(&ElasticStrainEnergy) - (char*)this;
    realOffsets["dissipatedEnergy"] = (char*)(&DissipatedEnergy) - (char*)this;
    realOffsets["stress"] = (char*)(&stress) - (char*)this;
    realOffsets["shearStress"] = (char*)(&stressShear) - (char*)this;
    realOffsets["normalVelocity"] = (char*)(&dxndt) - (char*)this;
    realOffsets["shearVelocity"] = (char*)(&dxsdt) - (char*)this;
    realOffsets["shear"] = (char*)(&xs) - (char*)this;
    realOffsets["normalApproach"] = (char*)(&normalApproach) - (char*)this;
    realOffsets["timestep"] = (char*)(&dt) - (char*)this;
    realOffsets["permeability"] = (char*)(&kappa) - (char*)this;
    R1TensorOffsets["shearStressVector"] = (char*)(&stressShearVector) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>&,
                                  std::map<std::string, realT>& realValues,
                                  std::map<std::string, R1Tensor>& R1TensorValues,
                                  std::map<std::string, R2Tensor>&,
                                  std::map<std::string, R2SymTensor>&  )
  {
    realValues["elasticStrainEnergy"] = ElasticStrainEnergy;
    realValues["dissipatedEnergy"] = DissipatedEnergy;
    realValues["stress"] = stress;
    realValues["shearStress"] = stressShear;
    realValues["normalVelocity"] = dxndt;
    realValues["shearVelocity"] = dxsdt;
    realValues["shear"] = xs;
    realValues["normalApproach"] = normalApproach;
    realValues["timestep"] = dt;
    realValues["permeability"] = kappa;
    R1TensorValues["shearStressVector"] = stressShearVector;
  }

  void Serialize(const localIndex index,
                 const unsigned int stride,
                 const localIndex elemNum,
                 array<array<integer>*>&,
                 array<array<real64>*>& realVars,
                 array<array<R1Tensor>*>& R1Vars,
                 array<array<R2Tensor>*>&,
                 array<array<R2SymTensor>*>&,
                 localIndex&,
                 localIndex& realVarCounts,
                 localIndex& R1TensorVarCounts,
                 localIndex&,
                 localIndex&   ) const
  {
    (void)index;
    (*(realVars[realVarCounts]))[elemNum] = ElasticStrainEnergy; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = DissipatedEnergy; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = stress; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = stressShear; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dxndt; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dxsdt; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = xs; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = normalApproach; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dt; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = kappa; realVarCounts += stride;
    (*(R1Vars[R1TensorVarCounts]))[elemNum] = stressShearVector; R1TensorVarCounts += stride;
  }


  void  Deserialize( const localIndex /*index*/,
                     const unsigned int stride,
                     const localIndex elemNum,
                     const array<array<integer>*>&,
                     const array<array<real64>*>& realVars,
                     const array<array<R1Tensor>*>& R1Vars,
                     const array<array<R2Tensor>*>&,
                     const array<array<R2SymTensor>*>&,
                     localIndex&,
                     localIndex& realVarCounts,
                     localIndex& R1TensorVarCounts,
                     localIndex&,
                     localIndex&   )
  {
    ElasticStrainEnergy = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    DissipatedEnergy = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    stress = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    stressShear = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dxndt = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dxsdt = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    xs = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    normalApproach = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dt = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    kappa = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    stressShearVector = (*(R1Vars[R1TensorVarCounts]))[elemNum]; R1TensorVarCounts += stride;
  }
  inline InterfaceBaseStateData&
  operator*=(const realT factor)
  {
    ElasticStrainEnergy *= factor;
    DissipatedEnergy *= factor;
    stress *= factor;
    stressShear *= factor;
    dxndt *= factor;
    dxsdt *= factor;
    xs *= factor;
    normalApproach *= factor;
    dt *= factor;
    kappa *= factor;
    stressShearVector *= factor;
    return *this;
  }

  inline InterfaceBaseStateData&
  operator=(const InterfaceBaseStateData& datum)
  {
    ElasticStrainEnergy = datum.ElasticStrainEnergy;
    DissipatedEnergy = datum.DissipatedEnergy;
    stress = datum.stress;
    stressShear = datum.stressShear;
    dxndt = datum.dxndt;
    dxsdt = datum.dxsdt;
    xs = datum.xs;
    normalApproach = datum.normalApproach;
    dt = datum.dt;
    kappa = datum.kappa;
    stressShearVector = datum.stressShearVector;
    return *this;
  }

  inline InterfaceBaseStateData&
  operator+=(const InterfaceBaseStateData& datum)
  {
    ElasticStrainEnergy += datum.ElasticStrainEnergy;
    DissipatedEnergy += datum.DissipatedEnergy;
    stress += datum.stress;
    stressShear += datum.stressShear;
    dxndt += datum.dxndt;
    dxsdt += datum.dxsdt;
    xs += datum.xs;
    normalApproach += datum.normalApproach;
    dt += datum.dt;
    kappa += datum.kappa;
    stressShearVector += datum.stressShearVector;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, InterfaceBaseStateData& p0, InterfaceBaseStateData& p1)
  {
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ElasticStrainEnergy, p0.ElasticStrainEnergy, p1.ElasticStrainEnergy);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, DissipatedEnergy, p0.DissipatedEnergy, p1.DissipatedEnergy);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stress, p0.stress, p1.stress);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressShear, p0.stressShear, p1.stressShear);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dxndt, p0.dxndt, p1.dxndt);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dxsdt, p0.dxsdt, p1.dxsdt);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, xs, p0.xs, p1.xs);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, normalApproach, p0.normalApproach, p1.normalApproach);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dt, p0.dt, p1.dt);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, kappa, p0.kappa, p1.kappa);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressShearVector, p0.stressShearVector, p1.stressShearVector);

  }

  void MapFromRegion(const InterfaceBaseStateData& p0, const InterfaceBaseStateData& p1, const realT fct0,
                     const realT fct1)
  {
    GeometryUtilities::MapFromRegion(p0.ElasticStrainEnergy, p1.ElasticStrainEnergy, fct0, fct1, ElasticStrainEnergy);
    GeometryUtilities::MapFromRegion(p0.DissipatedEnergy, p1.DissipatedEnergy, fct0, fct1, DissipatedEnergy);
    GeometryUtilities::MapFromRegion(p0.stress, p1.stress, fct0, fct1, stress);
    GeometryUtilities::MapFromRegion(p0.stressShear, p1.stressShear, fct0, fct1, stressShear);
    GeometryUtilities::MapFromRegion(p0.dxndt, p1.dxndt, fct0, fct1, dxndt);
    GeometryUtilities::MapFromRegion(p0.dxsdt, p1.dxsdt, fct0, fct1, dxsdt);
    GeometryUtilities::MapFromRegion(p0.xs, p1.xs, fct0, fct1, xs);
    GeometryUtilities::MapFromRegion(p0.normalApproach, p1.normalApproach, fct0, fct1, normalApproach);
    GeometryUtilities::MapFromRegion(p0.dt, p1.dt, fct0, fct1, dt);
    GeometryUtilities::MapFromRegion(p0.kappa, p1.kappa, fct0, fct1, kappa);
    GeometryUtilities::MapFromRegion(p0.stressShearVector, p1.stressShearVector, fct0, fct1, stressShearVector);

  }



};


//**********************************************************************************************************************
//**********************************************************************************************************************


class InterfaceBase : public ConstitutiveBase
{
public:
  const int m_paramSize;
  const int m_stateSize;


  typedef InterfaceBaseParameterData ParameterClass;
  typedef InterfaceBaseStateData     StateClass;

  inline std::string BaseName() { return "Interface"; }

  InterfaceBase( const int paramSize, const int stateSize );

  virtual ~InterfaceBase();

  virtual void ReadXML( TICPP::HierarchicalDataNode& node ) = 0;

  virtual void resize( const localIndex num ) = 0;

  virtual void resize( const localIndex num0,
                       const localIndex num1 ) = 0;

  virtual void insert( const localIndex num ) = 0;

  virtual void erase( const localIndex num ) = 0;

  virtual void InitializeStates( const localIndex ){}

  virtual const InterfaceBaseStateData* StateData( const localIndex index0,
                                                   const localIndex index1 ) const = 0;
  virtual       InterfaceBaseStateData* StateData( const localIndex index0,
                                                   const localIndex index1 )  = 0;

  virtual const InterfaceBaseParameterData* ParameterData( const localIndex index ) const = 0;
  virtual       InterfaceBaseParameterData* ParameterData( const localIndex index ) = 0;

  inline void IncrementPtr( const InterfaceBaseStateData* ptr ) const
  {
    ptr = reinterpret_cast<const InterfaceBaseStateData*>( reinterpret_cast<const char*>(ptr) + m_stateSize );
  }

  inline void IncrementPtr( const InterfaceBaseParameterData* ptr ) const
  {
    ptr = reinterpret_cast<const InterfaceBaseParameterData*>( reinterpret_cast<const char*>(ptr) + m_paramSize );
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
  Initialize( const localIndex index,
              const realT stressNormal,
              const realT stressShear);

  virtual void
  UpdateProperties(const localIndex, std::map<std::string, realT>&, std::map<std::string, realT>& );

  virtual realT
  SetPermeabilityTerm(const localIndex index);

  virtual realT
  StiffnessProjected( const localIndex index );

  virtual void
  StrainDrivenUpdate( const localIndex index );


protected:
  virtual realT ShearStrength(const InterfaceBaseParameterData& matParams,
                              InterfaceBaseStateData& matState) const;

  virtual void UpdateFriction(const InterfaceBaseParameterData&,
                              InterfaceBaseStateData&) const;

  virtual realT
  NormalStiffness(const InterfaceBaseParameterData&,
                  InterfaceBaseStateData&,
                  const realT,
                  const bool ) const;

  virtual void
  ThresholdToFailureSurface(const InterfaceBaseParameterData& matParams,
                            InterfaceBaseStateData& matState,
                            const realT dxs);


private:
  InterfaceBase();
  InterfaceBase( const InterfaceBase& );
  InterfaceBase& operator=( const InterfaceBase& );


};
#endif /* INTERFACEBASE_H_ */
