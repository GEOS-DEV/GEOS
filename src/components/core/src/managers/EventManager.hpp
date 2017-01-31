/*
 * EventManager.hpp
 *
 *  Created on: Oct 5, 2016
 *      Author: sherman
 */

#ifndef SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace pugi
{
class xml_node;
}

namespace geosx
{

class EventManager : public dataRepository::ManagedGroup
{
public:
  EventManager( std::string const & name,
                ManagedGroup * const parent );

  virtual ~EventManager();

  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;

  void ReadXML( pugi::xml_node const & problemNode ) override;

  void CheckEventTiming();

};

class SolverApplication : public dataRepository::ManagedGroup
{
public:
  SolverApplication( std::string const & name,
                     ManagedGroup * const parent );

  virtual ~SolverApplication();

  static string CatalogName() { return "SolverApplication"; }

  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_ */
