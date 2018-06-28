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

  struct viewKeyStruct
  {} m_ConstitutiveBaseViewKeys;

  struct groupKeyStruct
  {
    dataRepository::GroupKey StateData      = { "StateData" };
    dataRepository::GroupKey ParameterData  = { "ParameterData" };
  } m_ConstitutiveBaseGroupKeys;

  virtual viewKeyStruct       & viewKeys()        { return m_ConstitutiveBaseViewKeys; }
  virtual viewKeyStruct const & viewKeys() const  { return m_ConstitutiveBaseViewKeys; }

  virtual groupKeyStruct       & groupKeys()       { return m_ConstitutiveBaseGroupKeys; }
  virtual groupKeyStruct const & groupKeys() const { return m_ConstitutiveBaseGroupKeys; }

  ManagedGroup * GetParameterData()             { return this->GetGroup(m_ConstitutiveBaseGroupKeys.ParameterData); }
  ManagedGroup const * GetParameterData() const { return this->GetGroup(m_ConstitutiveBaseGroupKeys.ParameterData); }

  ManagedGroup * GetStateData()             { return this->GetGroup(m_ConstitutiveBaseGroupKeys.StateData); }
  ManagedGroup const * GetStateData() const { return this->GetGroup(m_ConstitutiveBaseGroupKeys.StateData); }

};



}
}
#endif
