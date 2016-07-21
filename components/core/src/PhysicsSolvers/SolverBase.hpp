/*
 * SolverBase.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#ifndef SOLVERBASE_HPP_
#define SOLVERBASE_HPP_


#include "codingUtilities/ObjectCatalogue.hpp"

#include "dataRepository/DataTypes.hpp"

#include "dataRepository/intrinsic/WrapperCollection.hpp"
#include <string>

namespace geosx
{

class SolverBase : public dataRepository::WrapperCollection
{
public:

  SolverBase( std::string const & name,
              WrapperCollection * const parent );

  virtual ~SolverBase();

  virtual void Registration( dataRepository::WrapperCollection& domain ) = 0;

  virtual void TimeStep( const real64& time_n,
                         const real64& dt,
                         const int cycleNumber,
                         dataRepository::WrapperCollection& domain ) = 0;

//  CATALOGUE( SolverBase, VA_LIST( std::string const & name ), VA_LIST( name ) )

  using CatalogInterface = objectcatalogue::CatalogInterface< SolverBase, std::string, WrapperCollection * const >;
  static CatalogInterface::CatalogType& GetCatalogue()
  {
    static CatalogInterface::CatalogType * const catalog = new CatalogInterface::CatalogType();
    return *catalog;
//    static CatalogInterface::CatalogueType catalogue;
//    return catalogue;
  }

private:

  SolverBase() = delete;
  SolverBase(const SolverBase&) = delete;
  SolverBase(const SolverBase&&) = delete;
  SolverBase& operator=(const SolverBase&) = delete;
  SolverBase& operator=(const SolverBase&&) = delete;
};





} /* namespace ANST */


#endif /* SOLVERBASE_HPP_ */
