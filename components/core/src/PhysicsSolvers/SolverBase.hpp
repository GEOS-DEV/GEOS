/*
 * SolverBase.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#ifndef SOLVERBASE_HPP_
#define SOLVERBASE_HPP_



#include "codingUtilities/ObjectCatalog.hpp"

#include "dataRepository/WrapperCollection.hpp"
#include <string>
#include "common/DataTypes.hpp"


namespace geosx
{

class SolverBase : public dataRepository::WrapperCollection
{
public:

  explicit SolverBase( std::string const & name,
                       WrapperCollection * const parent );

  virtual ~SolverBase();

  SolverBase() = default;
  SolverBase( SolverBase const & ) = default;
  SolverBase( SolverBase &&) = default;
  SolverBase& operator=( SolverBase const & ) = default;
  SolverBase& operator=( SolverBase&& ) = default;


  virtual void Registration( dataRepository::WrapperCollection& domain ) = 0;

  virtual void TimeStep( real64 const & time_n,
                         real64 const & dt,
                         int const cycleNumber,
                         dataRepository::WrapperCollection& domain ) = 0;


  using CatalogInterface = objectcatalog::CatalogInterface< SolverBase, std::string const &, WrapperCollection * const >;
  static CatalogInterface::CatalogType& GetCatalog();

private:

};



} /* namespace ANST */


#endif /* SOLVERBASE_HPP_ */
