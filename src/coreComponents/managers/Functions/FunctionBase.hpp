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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file FunctionBase.hpp
 */

#ifndef FUNCTIONBASE_HPP_
#define FUNCTIONBASE_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const inputVarNames("inputVarNames");
string const inputVarTypes("inputVarTypes");

}
}

/**
 * @class FunctionBase
 *
 * An object for interfacing with arbitrary N-dimensional functions.
 */
class FunctionBase : public dataRepository::ManagedGroup
{
public:
  /// Main constructor
  FunctionBase( const std::string& name,
                dataRepository::ManagedGroup * const parent );

  /// Destructor
  virtual ~FunctionBase() override;

  /// Documentation assignment
  virtual void FillDocumentationNode() override;

  /// Catalog name interface
  static string CatalogName() { return "FunctionBase"; }

  /// After reading the xml, call the function initialization
  virtual void ReadXML_PostProcess() override { InitializeFunction(); }

  /// Function initialization
  virtual void InitializeFunction(){}

  /// Test to see if the function is a 1D function of time
  integer isFunctionOfTime() const;

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  virtual void Evaluate( dataRepository::ManagedGroup const * const group,
                         real64 const time,
                         set<localIndex> const & set,
                         real64_array & result ) const = 0;

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   */
  virtual real64 Evaluate( real64 const * const input ) const = 0;

  // Setup catalog
  using CatalogInterface = cxx_utilities::CatalogInterface< FunctionBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  /*
   * @brief This generates statistics by applying a function to an object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @return An array holding the min, average, max values of the results
   */
  real64_array EvaluateStats( dataRepository::ManagedGroup const * const group,
                              real64 const time,
                              set<localIndex> const & set) const;

protected:
  template< typename LEAF >
  void EvaluateT( dataRepository::ManagedGroup const * const group,
                  real64 const time,
                  set<localIndex> const & set,
                  real64_array & result ) const;

};

/// Method to apply an function with an arbitrary type of output
template< typename LEAF >
void FunctionBase::EvaluateT( dataRepository::ManagedGroup const * const group,
                              real64 const time,
                              set<localIndex> const & set,
                              real64_array & result ) const
{
  real64 const * input_ptrs[4];

  string_array const & inputVarNames = this->getReference<string_array>( dataRepository::keys::inputVarNames );
  string_array const & inputVarTypes = this->getReference<string_array>( dataRepository::keys::inputVarTypes );

  localIndex const numVars = integer_conversion<localIndex>(inputVarNames.size());
  localIndex varSize[4];
  for( auto varIndex=0 ; varIndex<numVars ; ++varIndex )
  {
    string const & varName = inputVarNames[varIndex];
    string const & varTypeName = inputVarTypes[varIndex];

    if( varName=="time")
    {
      input_ptrs[varIndex] = &time;
      varSize[varIndex] = 1;
    }
    else
    {
      dataRepository::ViewWrapperBase const & vwb = *(group->getWrapperBase( varName ));
      std::type_index typeIndex = std::type_index(vwb.get_typeid());
      rtTypes::ApplyTypeLambda2( rtTypes::typeID(typeIndex), [&]( auto container_type, auto var_type ) -> void
        {
          using containerType = decltype(container_type);
          using varType = decltype(var_type);
          dataRepository::ViewWrapper<containerType> const & view =
            dynamic_cast< dataRepository::ViewWrapper<containerType> const & >(vwb);

          input_ptrs[varIndex] = reinterpret_cast<double const*>(view.dataPtr());
          varSize[varIndex] = sizeof(varType) / sizeof(double);
        });
    }
  }

  integer count=0;
  for( auto const & index : set )
  {
    double input[4];
    int c = 0;
    for( int a=0 ; a<numVars ; ++a )
    {
      for( int b=0 ; b<varSize[a] ; ++b )
      {
        input[c] = input_ptrs[a][index*varSize[a]+b];
        ++c;
      }
    }

    // TODO: Check this line to make sure it is correct
    result[count] = static_cast<LEAF const *>(this)->Evaluate(input);
    ++count;
  }

}
} /* namespace geosx */

#endif /* FUNCTIONBASE_HPP_ */
