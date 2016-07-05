/*
 * SolverBase.hpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#ifndef SOLVERBASE_HPP_
#define SOLVERBASE_HPP_

#include "dataRepository/DataTypes.hpp"
#include "codingUtilities/ObjectCatalogue.hpp"
#include <string>

namespace geosx
{
class DataObjectManager;

class SolverBase
{
public:
//  typedef ObjectCatalogueEntryBase< geosx::SolverBase, std::string > SolverFactory;

  SolverBase( std::string const & name );

  virtual ~SolverBase();


  virtual void RegisterDataObjects( DataObjectManager& domain ) = 0;

  virtual void TimeStep( const real64& time_n,
                         const real64& dt,
                         const int cycleNumber,
                         DataObjectManager& domain ) = 0;






//  class CatalogueEntryBase
//  {
//  public:
//    typedef std::unordered_map<std::string, CatalogueEntryBase*> CatalogueType;
//
//    CatalogueEntryBase()
//    {
//    }
//    virtual std::unique_ptr<SolverBase> Allocate( std::string const & name ) = 0;
//    virtual ~CatalogueEntryBase()
//    {
//    }
//
//    static CatalogueType& GetCatalogue()
//    {
//      static CatalogueType * const catalogue = new CatalogueType();
//      return *catalogue;
//    }
//
//    static std::unique_ptr<SolverBase> Factory( const std::string& objectTypeName, std::string const & name )
//    {
//      std::cout << "Creating solver of type: " << objectTypeName << std::endl;
//      CatalogueEntryBase* const entry = GetCatalogue().at( objectTypeName );
//      return entry->Allocate( name );
//    }
//
//  };
//
//  template<typename TYPE>
//  class CatalogueEntry : public CatalogueEntryBase
//  {
//  public:
//    CatalogueEntry() :
//        CatalogueEntryBase()
//    {
//      std::string name = TYPE::CatalogueName();
//      ( CatalogueEntryBase::GetCatalogue() )[name] = this;
//      std::cout << "Registered Solver: " << name << std::endl;
//    }
//
//    ~CatalogueEntry() final
//    {
//    }
//
//    virtual std::unique_ptr<SolverBase> Allocate( std::string const & name ) final
//    {
//      return std::unique_ptr<SolverBase>( new TYPE( name ) );
//    }
//    /// Compiler directive to simplify autoregistration
//
//  };
//
//  typedef CatalogueEntryBase SolverFactory;
//  #define REGISTER_SOLVER( SOLVER ) namespace{ SolverBase::CatalogueEntry<SOLVER> reg_; }


  CATALOGUE( SolverBase, VA_LIST( std::string const & name ), VA_LIST( name ) )

private:
  std::string m_name;




  SolverBase() = delete;
  SolverBase(const SolverBase&) = delete;
  SolverBase(const SolverBase&&) = delete;
  SolverBase& operator=(const SolverBase&) = delete;
  SolverBase& operator=(const SolverBase&&) = delete;
};





} /* namespace ANST */


#endif /* SOLVERBASE_HPP_ */
