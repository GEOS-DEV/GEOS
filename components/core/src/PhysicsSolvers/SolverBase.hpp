/*
 * SolverBase.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#ifndef SOLVERBASE_HPP_
#define SOLVERBASE_HPP_



#include "codingUtilities/ObjectCatalog.hpp"

#include "dataRepository/intrinsic/WrapperCollection.hpp"
#include <string>
#include "common/DataTypes.hpp"


namespace geosx
{

class SolverBase : public dataRepository::WrapperCollection
{
public:

  SolverBase( std::string const & name,
              WrapperCollection * const parent );

  virtual ~SolverBase();

  virtual void Registration( dataRepository::WrapperCollection& domain ) = 0;

  virtual void TimeStep( real64 const & time_n,
                         real64 const & dt,
                         int const cycleNumber,
                         dataRepository::WrapperCollection& domain ) = 0;


  using CatalogInterface = objectcatalog::CatalogInterface< SolverBase, std::string, WrapperCollection * const >;
  static CatalogInterface::CatalogType& GetCatalogue()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

private:

  SolverBase() = delete;
  SolverBase( SolverBase const & ) = delete;
  SolverBase( SolverBase &&) = delete;
  SolverBase& operator=( SolverBase const & ) = delete;
  SolverBase& operator=( SolverBase&& ) = delete;
};



} /* namespace ANST */


#endif /* SOLVERBASE_HPP_ */
