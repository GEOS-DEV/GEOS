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

#ifndef COHESIVEZONEBASE_H_
#define COHESIVEZONEBASE_H_

#include "Utilities/GeometryUtilities.h"
#include "Constitutive/ConstitutiveBase.h"

/*
 * CohesiveZoneBase.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */



//**********************************************************************************************************************
//**********************************************************************************************************************


class CohesiveZoneBaseParameterData
{

public:
  int unloadFlag;


  CohesiveZoneBaseParameterData():
    unloadFlag(0)
  {}

  CohesiveZoneBaseParameterData( const CohesiveZoneBaseParameterData& source):
    unloadFlag(source.unloadFlag)
  {}

  ~CohesiveZoneBaseParameterData() {}

  static void GetVariableCounts( localIndex&intVarCounts,
                                 localIndex&,
                                 localIndex&,
                                 localIndex&,
                                 localIndex&  )
  {
    intVarCounts = 1;

  }

  static void GetVariableNames( array<string>&intNames,
                                array<string>&,
                                array<string>&,
                                array<string>&,
                                array<string>&  )
  {
    intNames.push_back("unloadFlag");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>&intOffsets,
                                   std::map<std::string, size_t>&,
                                   std::map<std::string, size_t>&,
                                   std::map<std::string, size_t>&,
                                   std::map<std::string, size_t>&  ) const
  {
    intOffsets["unloadFlag"] = (char*)(&unloadFlag) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>&intValues,
                                  std::map<std::string, realT>&,
                                  std::map<std::string, R1Tensor>&,
                                  std::map<std::string, R2Tensor>&,
                                  std::map<std::string, R2SymTensor>&  )
  {
    intValues["unloadFlag"] = unloadFlag;
  }

  void Serialize(const localIndex index,
                 array<array<integer>*>& intVars,
                 array<array<real64>*>&,
                 array<array<R1Tensor>*>&,
                 array<array<R2Tensor>*>&,
                 array<array<R2SymTensor>*>&,
                 localIndex& intVarCounts,
                 localIndex&,
                 localIndex&,
                 localIndex&,
                 localIndex&   ) const
  {
    (*(intVars[intVarCounts]))[index] = unloadFlag; intVarCounts++;
  }


  void  Deserialize( const localIndex index,
                     const array<array<integer>*>& intVars,
                     const array<array<real64>*>&,
                     const array<array<R1Tensor>*>&,
                     const array<array<R2Tensor>*>&,
                     const array<array<R2SymTensor>*>&,
                     localIndex& intVarCounts,
                     localIndex&,
                     localIndex&,
                     localIndex&,
                     localIndex&   )
  {
    unloadFlag = (*(intVars[intVarCounts]))[index]; intVarCounts++;
  }
  inline CohesiveZoneBaseParameterData&
  operator*=(const realT factor)
  {
    unloadFlag *= factor;
    return *this;
  }

  inline CohesiveZoneBaseParameterData&
  operator=(const CohesiveZoneBaseParameterData& datum)
  {
    unloadFlag = datum.unloadFlag;
    return *this;
  }

  inline CohesiveZoneBaseParameterData&
  operator+=(const CohesiveZoneBaseParameterData& datum)
  {
    unloadFlag += datum.unloadFlag;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, CohesiveZoneBaseParameterData& p0, CohesiveZoneBaseParameterData& p1)
  {
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, unloadFlag, p0.unloadFlag, p1.unloadFlag);

  }

  void MapFromRegion(const CohesiveZoneBaseParameterData& p0, const CohesiveZoneBaseParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    GeometryUtilities::MapFromRegion(p0.unloadFlag, p1.unloadFlag, fct0, fct1, unloadFlag);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {}
  virtual void PostReadXML( TICPP::HierarchicalDataNode& ) {}



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class CohesiveZoneBaseStateData
{

public:
  realT separationCoeff;
  R1Tensor traction;
  R2Tensor stiffness;


  CohesiveZoneBaseStateData():
    separationCoeff(0),
    traction(0)
  {}

  CohesiveZoneBaseStateData( const CohesiveZoneBaseStateData& source):
    separationCoeff(source.separationCoeff),
    traction(source.traction),
    stiffness(source.stiffness)
  {}

  ~CohesiveZoneBaseStateData() {}

  static void GetVariableCounts( localIndex&,
                                 localIndex& realVarCounts,
                                 localIndex& R1TensorVarCounts,
                                 localIndex& R2TensorVarCounts,
                                 localIndex&  )
  {
    realVarCounts = 1;
    R1TensorVarCounts = 1;
    R2TensorVarCounts = 1;

  }

  static void GetVariableNames( array<string>&,
                                array<string>& realNames,
                                array<string>& R1TensorNames,
                                array<string>& R2TensorNames,
                                array<string>&  )
  {
    realNames.push_back("separationCoeff");
    R1TensorNames.push_back("traction");
    R2TensorNames.push_back("stiffness");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>&,
                                   std::map<std::string, size_t>& realOffsets,
                                   std::map<std::string, size_t>& R1TensorOffsets,
                                   std::map<std::string, size_t>& R2TensorOffsets,
                                   std::map<std::string, size_t>&  ) const
  {
    realOffsets["separationCoeff"] = (char*)(&separationCoeff) - (char*)this;
    R1TensorOffsets["traction"] = (char*)(&traction) - (char*)this;
    R2TensorOffsets["stiffness"] = (char*)(&stiffness) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>&,
                                  std::map<std::string, realT>& realValues,
                                  std::map<std::string, R1Tensor>& R1TensorValues,
                                  std::map<std::string, R2Tensor>& R2TensorValues,
                                  std::map<std::string, R2SymTensor>&  )
  {
    realValues["separationCoeff"] = separationCoeff;
    R1TensorValues["traction"] = traction;
    R2TensorValues["stiffness"] = stiffness;
  }

  void Serialize(const localIndex index,
                 const unsigned int stride,
                 const localIndex elemNum,
                 array<array<integer>*>&,
                 array<array<real64>*>& realVars,
                 array<array<R1Tensor>*>& R1Vars,
                 array<array<R2Tensor>*>& R2Vars,
                 array<array<R2SymTensor>*>&,
                 localIndex&,
                 localIndex& realVarCounts,
                 localIndex& R1TensorVarCounts,
                 localIndex& R2TensorVarCounts,
                 localIndex&   ) const
  {
    (*(realVars[realVarCounts]))[elemNum] = separationCoeff; realVarCounts += stride;
    (*(R1Vars[R1TensorVarCounts]))[elemNum] = traction; R1TensorVarCounts += stride;
    (*(R2Vars[R2TensorVarCounts]))[elemNum] = stiffness; R2TensorVarCounts += stride;
  }


  void  Deserialize( const localIndex index,
                     const unsigned int stride,
                     const localIndex elemNum,
                     const array<array<integer>*>&,
                     const array<array<real64>*>& realVars,
                     const array<array<R1Tensor>*>& R1Vars,
                     const array<array<R2Tensor>*>& R2Vars,
                     const array<array<R2SymTensor>*>&,
                     localIndex&,
                     localIndex& realVarCounts,
                     localIndex& R1TensorVarCounts,
                     localIndex& R2TensorVarCounts,
                     localIndex&   )
  {
    separationCoeff = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    traction = (*(R1Vars[R1TensorVarCounts]))[elemNum]; R1TensorVarCounts += stride;
    stiffness = (*(R2Vars[R2TensorVarCounts]))[elemNum]; R2TensorVarCounts += stride;
  }
  inline CohesiveZoneBaseStateData&
  operator*=(const realT factor)
  {
    separationCoeff *= factor;
    traction *= factor;
    stiffness *= factor;
    return *this;
  }

  inline CohesiveZoneBaseStateData&
  operator=(const CohesiveZoneBaseStateData& datum)
  {
    separationCoeff = datum.separationCoeff;
    traction = datum.traction;
    stiffness = datum.stiffness;
    return *this;
  }

  inline CohesiveZoneBaseStateData&
  operator+=(const CohesiveZoneBaseStateData& datum)
  {
    separationCoeff += datum.separationCoeff;
    traction += datum.traction;
    stiffness += datum.stiffness;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, CohesiveZoneBaseStateData& p0, CohesiveZoneBaseStateData& p1)
  {
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, separationCoeff, p0.separationCoeff, p1.separationCoeff);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, traction, p0.traction, p1.traction);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stiffness, p0.stiffness, p1.stiffness);

  }

  void MapFromRegion(const CohesiveZoneBaseStateData& p0, const CohesiveZoneBaseStateData& p1, const realT fct0,
                     const realT fct1)
  {
    GeometryUtilities::MapFromRegion(p0.separationCoeff, p1.separationCoeff, fct0, fct1, separationCoeff);
    GeometryUtilities::MapFromRegion(p0.traction, p1.traction, fct0, fct1, traction);
    GeometryUtilities::MapFromRegion(p0.stiffness, p1.stiffness, fct0, fct1, stiffness);

  }



};


//**********************************************************************************************************************
//**********************************************************************************************************************


class CohesiveZoneBase : public ConstitutiveBase
{
public:
  const int m_paramSize;
  const int m_stateSize;


  typedef CohesiveZoneBaseParameterData ParameterClass;
  typedef CohesiveZoneBaseStateData     StateClass;

  inline std::string BaseName() { return "CohesiveZone"; }

  CohesiveZoneBase( const int paramSize, const int stateSize );

  virtual ~CohesiveZoneBase();

  virtual void ReadXML( TICPP::HierarchicalDataNode& node ) = 0;

  virtual void resize( const localIndex num ) = 0;

  virtual void resize( const localIndex num0,
                       const localIndex num1 ) = 0;

  virtual void insert( const localIndex num ) = 0;

  virtual void erase( const localIndex num ) = 0;

  virtual void InitializeStates( const localIndex index ){}

  virtual const CohesiveZoneBaseStateData* StateData( const localIndex index0,
                                                      const localIndex index1 ) const = 0;
  virtual       CohesiveZoneBaseStateData* StateData( const localIndex index0,
                                                      const localIndex index1 )  = 0;

  virtual const CohesiveZoneBaseParameterData* ParameterData( const localIndex index ) const = 0;
  virtual       CohesiveZoneBaseParameterData* ParameterData( const localIndex index ) = 0;

  inline void IncrementPtr( const CohesiveZoneBaseStateData* ptr ) const
  {
    ptr = reinterpret_cast<const CohesiveZoneBaseStateData*>( reinterpret_cast<const char*>(ptr) + m_stateSize );
  }

  inline void IncrementPtr( const CohesiveZoneBaseParameterData* ptr ) const
  {
    ptr = reinterpret_cast<const CohesiveZoneBaseParameterData*>( reinterpret_cast<const char*>(ptr) + m_paramSize );
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

  virtual int
  UpdateCohesiveZone( const localIndex index,
                      const R1Tensor& gap,
                      const R1Tensor& N,
                      const std::pair< ElementRegionT*, localIndex >& elem0,
                      const std::pair< ElementRegionT*, localIndex >& elem1,
                      R1Tensor& traction,
                      R2Tensor& stiffness );



private:
  CohesiveZoneBase();
  CohesiveZoneBase( const CohesiveZoneBase& );
  CohesiveZoneBase& operator=( const CohesiveZoneBase& );


};
#endif /* COHESIVEZONEBASE_H_ */
