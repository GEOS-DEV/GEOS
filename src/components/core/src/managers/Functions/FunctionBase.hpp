/*
 * FunctionBase.hpp
 *
 *  Created on: July 6, 2017
 *      Author: sherman
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

class FunctionBase : public dataRepository::ManagedGroup
{
public:
  FunctionBase( const std::string& name,
                dataRepository::ManagedGroup * const parent );

  virtual ~FunctionBase();

  virtual void FillDocumentationNode() override;

  static string CatalogName() { return "FunctionBase"; }
  virtual void InitializeFunction(){}

  integer isFunctionOfTime() const;

  virtual void Evaluate( dataRepository::ManagedGroup const * const group,
                         real64 const time,
                         lSet const & set,
                         real64_array & result ) const = 0;


  virtual real64 Evaluate( real64 const * const input ) const = 0;

  // Setup catalog
  using CatalogInterface = cxx_utilities::CatalogInterface< FunctionBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

protected:
  template< typename LEAF >
  void EvaluateT( dataRepository::ManagedGroup const * const group,
                  real64 const time,
                  lSet const & set,
                  real64_array & result ) const;

};


template< typename LEAF >
void FunctionBase::EvaluateT( dataRepository::ManagedGroup const * const group,
                              real64 const time,
                              lSet const & set,
                              real64_array & result ) const
{
  real64 const * input_ptrs[4];

  string_array const & inputVarNames = this->getReference<string_array>( dataRepository::keys::inputVarNames );
  string_array const & inputVarTypes = this->getReference<string_array>( dataRepository::keys::inputVarTypes );

  localIndex const numVars = inputVarNames.size();
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
      rtTypes::ApplyTypeLambda1( rtTypes::typeID(typeIndex), [&]( auto type ) -> void
        {
          using varType = decltype(type);
          dataRepository::ViewWrapper<varType> const & view = dynamic_cast< dataRepository::ViewWrapper<varType> const & >(vwb);

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
        input[c] = input_ptrs[a][index+b];
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
