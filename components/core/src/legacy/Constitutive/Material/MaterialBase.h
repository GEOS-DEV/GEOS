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

#ifndef MATERIALBASE_H_
#define MATERIALBASE_H_

#include "Utilities/GeometryUtilities.h"
#include "Constitutive/ConstitutiveBase.h"

/*
 * MaterialBase.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class MaterialBaseParameterData 
{

public:
  realT init_density;
  realT E;
  realT Nu;
  realT init_shearModulus;
  realT Lame;


  MaterialBaseParameterData():
    init_density(0),
    E(0),
    Nu(0),
    init_shearModulus(0),
    Lame(0)
  {}

  MaterialBaseParameterData( const MaterialBaseParameterData& source):
    init_density(source.init_density),
    E(source.E),
    Nu(source.Nu),
    init_shearModulus(source.init_shearModulus),
    Lame(source.Lame)
  {}

  ~MaterialBaseParameterData() {}

  static void GetVariableCounts( localIndex&,
                                 localIndex& realVarCounts,
                                 localIndex& ,
                                 localIndex& ,
                                 localIndex&  )
  {
    realVarCounts = 5;

  }

  static void GetVariableNames( sArray1d&,
                                 sArray1d& realNames,
                                 sArray1d& ,
                                 sArray1d& ,
                                 sArray1d&  )
  {
    realNames.push_back("Density");
    realNames.push_back("E");
    realNames.push_back("Nu");
    realNames.push_back("ShearModulus");
    realNames.push_back("Lame");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>&,
                                 std::map<std::string, size_t>& realOffsets,
                                 std::map<std::string, size_t>& ,
                                 std::map<std::string, size_t>& ,
                                 std::map<std::string, size_t>&  ) const
  {
    realOffsets["Density"] = (char*)(&init_density) - (char*)this;
    realOffsets["E"] = (char*)(&E) - (char*)this;
    realOffsets["Nu"] = (char*)(&Nu) - (char*)this;
    realOffsets["ShearModulus"] = (char*)(&init_shearModulus) - (char*)this;
    realOffsets["Lame"] = (char*)(&Lame) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>&,
                                 std::map<std::string, realT>& realValues,
                                 std::map<std::string, R1Tensor>& ,
                                 std::map<std::string, R2Tensor>& ,
                                 std::map<std::string, R2SymTensor>&  )
  {
    realValues["Density"] = init_density;
    realValues["E"] = E;
    realValues["Nu"] = Nu;
    realValues["ShearModulus"] = init_shearModulus;
    realValues["Lame"] = Lame;
  }

  void Serialize(const localIndex index,
                  Array1dT<iArray1d*>& ,
                  Array1dT<rArray1d*>& realVars,
                  Array1dT<Array1dT<R1Tensor>*>& ,
                  Array1dT<Array1dT<R2Tensor>*>& ,
                  Array1dT<Array1dT<R2SymTensor>*>& ,
                  localIndex& ,
                  localIndex& realVarCounts,
                  localIndex& ,
                  localIndex& ,
                  localIndex&   ) const
  {
    (*(realVars[realVarCounts]))[index] = init_density; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = E; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = Nu; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = init_shearModulus; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = Lame; realVarCounts++;
  }


  void  Deserialize( const localIndex index,
                     const Array1dT<iArray1d*>& ,
                  const Array1dT<rArray1d*>& realVars,
                  const Array1dT<Array1dT<R1Tensor>*>& ,
                  const Array1dT<Array1dT<R2Tensor>*>& ,
                  const Array1dT<Array1dT<R2SymTensor>*>& ,
                  localIndex& ,
                  localIndex& realVarCounts,
                  localIndex& ,
                  localIndex& ,
                  localIndex&   )
  {
    init_density = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    E = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    Nu = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    init_shearModulus = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    Lame = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline MaterialBaseParameterData&
  operator*=(const realT factor)
  {
    init_density *= factor;
    E *= factor;
    Nu *= factor;
    init_shearModulus *= factor;
    Lame *= factor;
    return *this;
  }

  inline MaterialBaseParameterData&
  operator=(const MaterialBaseParameterData& datum)
  {
    init_density = datum.init_density;
    E = datum.E;
    Nu = datum.Nu;
    init_shearModulus = datum.init_shearModulus;
    Lame = datum.Lame;
    return *this;
  }

  inline MaterialBaseParameterData&
  operator+=(const MaterialBaseParameterData& datum)
  {
    init_density += datum.init_density;
    E += datum.E;
    Nu += datum.Nu;
    init_shearModulus += datum.init_shearModulus;
    Lame += datum.Lame;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, MaterialBaseParameterData& p0, MaterialBaseParameterData& p1)
  {
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, init_density, p0.init_density, p1.init_density);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, E, p0.E, p1.E);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, Nu, p0.Nu, p1.Nu);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, init_shearModulus, p0.init_shearModulus, p1.init_shearModulus);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, Lame, p0.Lame, p1.Lame);

  }

  void MapFromRegion(const MaterialBaseParameterData& p0, const MaterialBaseParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    GeometryUtilities::MapFromRegion(p0.init_density, p1.init_density, fct0, fct1, init_density);
    GeometryUtilities::MapFromRegion(p0.E, p1.E, fct0, fct1, E);
    GeometryUtilities::MapFromRegion(p0.Nu, p1.Nu, fct0, fct1, Nu);
    GeometryUtilities::MapFromRegion(p0.init_shearModulus, p1.init_shearModulus, fct0, fct1, init_shearModulus);
    GeometryUtilities::MapFromRegion(p0.Lame, p1.Lame, fct0, fct1, Lame);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    init_density = node.GetAttributeOrDefault("Density", 0.0);
    E = node.GetAttributeOrDefault("E", 0.0);
    Nu = node.GetAttributeOrDefault("Nu", 0.0);
    init_shearModulus = node.GetAttributeOrDefault("ShearModulus", 0.0);
    Lame = node.GetAttributeOrDefault("Lame", 0.0);

  }
  virtual void PostReadXML( const TICPP::HierarchicalDataNode& ) {}


  inline virtual void IncreaseDamageParameter(realT tmp)
  {}


};

//**********************************************************************************************************************
//**********************************************************************************************************************


class MaterialBaseStateData 
{

public:
  realT StressPower;
  realT DissipatedEnergy;
  realT ElasticStrainEnergy;
  realT pressure;
  //realT density;
  realT BulkModulus;
  realT ShearModulus;
  R2SymTensor devStress;


  MaterialBaseStateData():
    StressPower(0),
    DissipatedEnergy(0),
    ElasticStrainEnergy(0),
    pressure(0),
   // density(0),
    BulkModulus(0),
    ShearModulus(0),
    devStress(0)
  {}

  MaterialBaseStateData( const MaterialBaseStateData& source):
    StressPower(source.StressPower),
    DissipatedEnergy(source.DissipatedEnergy),
    ElasticStrainEnergy(source.ElasticStrainEnergy),
    pressure(source.pressure),
  //  density(source.density),
    BulkModulus(source.BulkModulus),
    ShearModulus(source.ShearModulus),
    devStress(source.devStress)
  {}

  ~MaterialBaseStateData() {}

  void
  TotalStress(R2SymTensor& totalStress) const;

  void
  RotateState( const R2Tensor& Rot );


  virtual void SetDensity(realT rho){};

  virtual realT GetDensity() const{return 0.0; };
  virtual void SetSpecificInternalEnergy(realT rho){};
  virtual realT GetSpecificInternalEnergy() const{return 0;};

  static void GetVariableCounts( localIndex&,
                                 localIndex& realVarCounts,
                                 localIndex& ,
                                 localIndex& ,
                                 localIndex& R2SymTensorVarCounts )
  {
    realVarCounts = 6;
    R2SymTensorVarCounts = 1;

  }

  static void GetVariableNames( sArray1d&,
                                 sArray1d& realNames,
                                 sArray1d& ,
                                 sArray1d& ,
                                 sArray1d& R2SymTensorNames )
  {
    realNames.push_back("StressPower");
    realNames.push_back("DissipatedEnergy");
    realNames.push_back("ElasticStrainEnergy");
    realNames.push_back("pressure");
   // realNames.push_back("density");
    realNames.push_back("BulkModulusCurrent");
    realNames.push_back("ShearModulusCurrent");
    R2SymTensorNames.push_back("devStress");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>&,
                                 std::map<std::string, size_t>& realOffsets,
                                 std::map<std::string, size_t>& ,
                                 std::map<std::string, size_t>& ,
                                 std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    realOffsets["StressPower"] = (char*)(&StressPower) - (char*)this;
    realOffsets["DissipatedEnergy"] = (char*)(&DissipatedEnergy) - (char*)this;
    realOffsets["ElasticStrainEnergy"] = (char*)(&ElasticStrainEnergy) - (char*)this;
    realOffsets["pressure"] = (char*)(&pressure) - (char*)this;
  //  realOffsets["density"] = (char*)(&pressure) - (char*)this;
    realOffsets["BulkModulusCurrent"] = (char*)(&BulkModulus) - (char*)this;
    realOffsets["ShearModulusCurrent"] = (char*)(&ShearModulus) - (char*)this;
    R2SymTensorOffsets["devStress"] = (char*)(&devStress) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>&,
                                 std::map<std::string, realT>& realValues,
                                 std::map<std::string, R1Tensor>& ,
                                 std::map<std::string, R2Tensor>& ,
                                 std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    realValues["StressPower"] = StressPower;
    realValues["DissipatedEnergy"] = DissipatedEnergy;
    realValues["ElasticStrainEnergy"] = ElasticStrainEnergy;
    realValues["pressure"] = pressure;
    //realOffsets["density"] = density;
    realValues["BulkModulusCurrent"] = BulkModulus;
    realValues["ShearModulusCurrent"] = ShearModulus;
    R2SymTensorValues["devStress"] = devStress;
  }

  void Serialize(const localIndex index,
                  const unsigned int stride,
                  const localIndex elemNum,
                  Array1dT<iArray1d*>& ,
                  Array1dT<rArray1d*>& realVars,
                  Array1dT<Array1dT<R1Tensor>*>& ,
                  Array1dT<Array1dT<R2Tensor>*>& ,
                  Array1dT<Array1dT<R2SymTensor>*>& R2SymVars,
                  localIndex& ,
                  localIndex& realVarCounts,
                  localIndex& ,
                  localIndex& ,
                  localIndex& R2SymTensorVarCounts  ) const
  {
    (*(realVars[realVarCounts]))[elemNum] = StressPower; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = DissipatedEnergy; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = ElasticStrainEnergy; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = pressure; realVarCounts += stride;
  //  (*(realVars[realVarCounts]))[elemNum] = density; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = BulkModulus; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = ShearModulus; realVarCounts += stride;
    (*(R2SymVars[R2SymTensorVarCounts]))[elemNum] = devStress; R2SymTensorVarCounts += stride;
  }


  void  Deserialize( const localIndex index,
                  const unsigned int stride,
                  const localIndex elemNum,
                     const Array1dT<iArray1d*>& ,
                  const Array1dT<rArray1d*>& realVars,
                  const Array1dT<Array1dT<R1Tensor>*>& ,
                  const Array1dT<Array1dT<R2Tensor>*>& ,
                  const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars,
                  localIndex& ,
                  localIndex& realVarCounts,
                  localIndex& ,
                  localIndex& ,
                  localIndex& R2SymTensorVarCounts  )
  {
    StressPower = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    DissipatedEnergy = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    ElasticStrainEnergy = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    pressure = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
   // density = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    BulkModulus = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    ShearModulus = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    devStress = (*(R2SymVars[R2SymTensorVarCounts]))[elemNum]; R2SymTensorVarCounts += stride;
  }
  inline MaterialBaseStateData&
  operator*=(const realT factor)
  {
    StressPower *= factor;
    DissipatedEnergy *= factor;
    ElasticStrainEnergy *= factor;
    pressure *= factor;
   // density *= factor;
    BulkModulus *= factor;
    ShearModulus *= factor;
    devStress *= factor;
    return *this;
  }

  inline MaterialBaseStateData&
  operator=(const MaterialBaseStateData& datum)
  {
    StressPower = datum.StressPower;
    DissipatedEnergy = datum.DissipatedEnergy;
    ElasticStrainEnergy = datum.ElasticStrainEnergy;
    pressure = datum.pressure;
   // density = datum.density;
    BulkModulus = datum.BulkModulus;
    ShearModulus = datum.ShearModulus;
    devStress = datum.devStress;
    return *this;
  }

  inline MaterialBaseStateData&
  operator+=(const MaterialBaseStateData& datum)
  {
    StressPower += datum.StressPower;
    DissipatedEnergy += datum.DissipatedEnergy;
    ElasticStrainEnergy += datum.ElasticStrainEnergy;
    pressure += datum.pressure;
    // density += datum.density;
    BulkModulus += datum.BulkModulus;
    ShearModulus += datum.ShearModulus;
    devStress += datum.devStress;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, MaterialBaseStateData& p0, MaterialBaseStateData& p1)
  {
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, StressPower, p0.StressPower, p1.StressPower);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, DissipatedEnergy, p0.DissipatedEnergy, p1.DissipatedEnergy);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ElasticStrainEnergy, p0.ElasticStrainEnergy, p1.ElasticStrainEnergy);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, pressure, p0.pressure, p1.pressure);
   // GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, density, p0.density, p1.density);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, BulkModulus, p0.BulkModulus, p1.BulkModulus);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ShearModulus, p0.ShearModulus, p1.ShearModulus);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, devStress, p0.devStress, p1.devStress);

  }

  void MapFromRegion(const MaterialBaseStateData& p0, const MaterialBaseStateData& p1, const realT fct0,
                     const realT fct1)
  {
    GeometryUtilities::MapFromRegion(p0.StressPower, p1.StressPower, fct0, fct1, StressPower);
    GeometryUtilities::MapFromRegion(p0.DissipatedEnergy, p1.DissipatedEnergy, fct0, fct1, DissipatedEnergy);
    GeometryUtilities::MapFromRegion(p0.ElasticStrainEnergy, p1.ElasticStrainEnergy, fct0, fct1, ElasticStrainEnergy);
    GeometryUtilities::MapFromRegion(p0.pressure, p1.pressure, fct0, fct1, pressure);
    // GeometryUtilities::MapFromRegion(p0.density, p1.density, fct0, fct1, density);
    GeometryUtilities::MapFromRegion(p0.BulkModulus, p1.BulkModulus, fct0, fct1, BulkModulus);
    GeometryUtilities::MapFromRegion(p0.ShearModulus, p1.ShearModulus, fct0, fct1, ShearModulus);
    GeometryUtilities::MapFromRegion(p0.devStress, p1.devStress, fct0, fct1, devStress);

  }



};


//**********************************************************************************************************************
//**********************************************************************************************************************


class MaterialBase: public ConstitutiveBase
{
public:
  const int m_paramSize;
  const int m_stateSize;

  
  typedef MaterialBaseParameterData ParameterClass;
  typedef MaterialBaseStateData     StateClass;
  
  inline std::string BaseName() { return "Material"; }

  MaterialBase( const int paramSize, const int stateSize );

  virtual ~MaterialBase();
  
  virtual void ReadXML( TICPP::HierarchicalDataNode& node ) = 0;

  virtual void resize( const localIndex num ) = 0;
  
  virtual void resize( const localIndex num0,
                       const localIndex num1 ) = 0;
  
  virtual void insert( const localIndex num ) = 0;

  virtual void erase( const localIndex num ) = 0;
  
  virtual void InitializeStates( const localIndex index ){}

  virtual const MaterialBaseStateData* StateData( const localIndex index0,
                                                  const localIndex index1 ) const = 0;
  virtual       MaterialBaseStateData* StateData( const localIndex index0,
                                                  const localIndex index1 )  = 0;

  virtual const MaterialBaseParameterData* ParameterData( const localIndex index ) const = 0;
  virtual       MaterialBaseParameterData* ParameterData( const localIndex index ) = 0;
  
  inline void IncrementPtr( const MaterialBaseStateData* ptr ) const
  {
    ptr = reinterpret_cast<const MaterialBaseStateData*>( reinterpret_cast<const char*>(ptr) + m_stateSize );
  }

  inline void IncrementPtr( const MaterialBaseParameterData* ptr ) const
  {
    ptr = reinterpret_cast<const MaterialBaseParameterData*>( reinterpret_cast<const char*>(ptr) + m_paramSize );
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
  StrainDrivenUpdateMember( const localIndex index0,
                              const localIndex index1,
                              const R2SymTensorT < 3 >& Ddt,
                              const R2TensorT < 3 >& L,
                              const R2Tensor& Rot,
                              const realT dt );

  virtual void
  StrainDrivenUpdateMember( const localIndex index0,
                              const localIndex index1,
                              const R2SymTensorT < 3 >& Ddt,
                              const R2TensorT < 3 >& L,
                              const R2Tensor& Rot,
                              const realT& volume_n,
                              const realT& volume_np1,
                              const realT dt);

  virtual void
  MeanPressureDevStress( const localIndex index,
                                       realT& pressure, R2SymTensor& devStress) const;

  template< typename LeafClass > void
  MeanPressureDevStressFromDerived( const localIndex index,
                                    realT& pressure, R2SymTensor& devStress) const
  {
    pressure = 0;
    devStress = 0;

    const LeafClass& dthis = static_cast<const LeafClass&>(*this);
    const typename LeafClass::StateClass* const state = dthis.m_stateData[index];

    for( localIndex b=0 ; b<dthis.NumStateIndex1() ; ++b )
    {
      devStress += state[b].devStress;
      pressure += state[b].pressure;
    }

    if(dthis.NumStateIndex1() > 0)
    {
      pressure /= dthis.NumStateIndex1();
      devStress /= dthis.NumStateIndex1();
    }
  }



  virtual bool NeedsDensity(){return false;};
  virtual bool NeedsSpecificInternalEnergy(){return false;};


  
private:
  MaterialBase();
  MaterialBase( const MaterialBase& );
  MaterialBase& operator=( const MaterialBase& );
  
  
};
#endif /* MATERIALBASE_H_ */
