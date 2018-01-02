/*
 * SolverBase.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#ifndef SOLVERBASE_HPP_
#define SOLVERBASE_HPP_



#include <string>
#include <limits>

#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"
#include "../dataRepository/ManagedGroup.hpp"
#include "common/DataTypes.hpp"
#include "mesh/MeshBody.hpp"


namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const courant = "courant";
string const maxDt   = "maxDt";
}
}

class SolverBase : public dataRepository::ManagedGroup
{
public:

  explicit SolverBase( std::string const & name,
                       ManagedGroup * const parent );

  virtual ~SolverBase();

  static string CatalogName() { return "SolverBase"; }

  SolverBase() = default;
  SolverBase( SolverBase const & ) = default;
  SolverBase( SolverBase &&) = default;
  SolverBase& operator=( SolverBase const & ) = default;
  SolverBase& operator=( SolverBase&& ) = default;


//  virtual void Registration( dataRepository::WrapperCollection& domain );


  virtual void TimeStep( real64 const & time_n,
                         real64 const & dt,
                         int const cycleNumber,
                         dataRepository::ManagedGroup * domain );

  // virtual void TimeStep( real64 const & time_n,
  //                        real64 const & dt,
  //                        int const cycleNumber,
  //                        dataRepository::ManagedGroup * domain ) = 0;

  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  using CatalogInterface = cxx_utilities::CatalogInterface< SolverBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

private:

};



} /* namespace ANST */


#endif /* SOLVERBASE_HPP_ */
