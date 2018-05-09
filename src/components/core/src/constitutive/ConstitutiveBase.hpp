/*
 * ConstitutiveBase.hpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#ifndef CONSTITUTIVEBASE_HPP_
#define CONSTITUTIVEBASE_HPP_

#include "common/DataTypes.hpp"
#include "ObjectCatalog.hpp"
#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"
#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const stateData("StateData");
string const parameterData("ParameterData");
}
}

namespace constitutive
{


class ConstitutiveBase : public dataRepository::ManagedGroup
{
public:


  ConstitutiveBase( std::string const & name,
                    ManagedGroup * const parent  );

  virtual ~ConstitutiveBase() override;


  virtual void SetParamStatePointers( void *& ) = 0;


  typedef void (*UpdateFunctionPointer)( R2SymTensor const & D,
                                        R2Tensor const & Rot,
                                        localIndex const i,
                                        void * dataPtrs,
                                        integer const systemAssembleFlag);

  virtual UpdateFunctionPointer GetStateUpdateFunctionPointer( ) = 0;

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const = 0;


  virtual R2SymTensor StateUpdatePoint( R2SymTensor const & D,
                                        R2Tensor const & Rot,
                                        localIndex const i,
                                        integer const systemAssembleFlag ) = 0;


  virtual void FillDocumentationNode() override = 0;

  virtual void resize( localIndex ) override;

  void SetVariableParameters();

  virtual void GetStiffness( realT c[6][6]) const = 0;


  using CatalogInterface = cxx_utilities::CatalogInterface< ConstitutiveBase, std::string const &, ManagedGroup * const >;
  static typename CatalogInterface::CatalogType& GetCatalog();

  struct ViewKeyStruct
  {} viewKeys;

  struct GroupKeyStruct
  {
    dataRepository::GroupKey StateData      = { "StateData" };
    dataRepository::GroupKey ParameterData  = { "ParameterData" };
  } groupKeys;

  ManagedGroup * GetParameterData()             { return this->GetGroup(groupKeys.ParameterData); }
  ManagedGroup const * GetParameterData() const { return this->GetGroup(groupKeys.ParameterData); }

  ManagedGroup * GetStateData()             { return this->GetGroup(groupKeys.StateData); }
  ManagedGroup const * GetStateData() const { return this->GetGroup(groupKeys.StateData); }

};



}
}
#endif
