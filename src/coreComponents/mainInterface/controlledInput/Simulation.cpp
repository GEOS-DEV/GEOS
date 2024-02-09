#include "Simulation.hpp"

#include "Misc.hpp"

namespace geos::input
{

using namespace pugi;

void Simulation::setBegin( string const & begin )
{
  m_begin = begin;
}

void Simulation::setEnd( string const & end )
{
  m_end = end;
}

void Simulation::setSolver( std::vector< std::shared_ptr< solvers::Solver > > const & solvers )
{
  m_solvers = solvers;
}

void Simulation::setNumericalStrategy( std::shared_ptr< numericalStrategies::NumericalStrategies > ns )
{
  m_ns = ns;
}

void Simulation::fillProblemXmlNode( xml_node & problemNode, std::vector< string > const & defaultDomains ) const
{
  xml_node eventsNode = problemNode.select_node( "Events" ).node();
  fillEventsXmlNode( eventsNode );

  for( auto const & solver: m_solvers )
  {
    solver->fillProblemXmlNode( problemNode, defaultDomains );
  }

  m_ns->fillProblemXmlNode( problemNode );
}

void Simulation::fillEventsXmlNode( xml_node & eventsNode ) const
{
  eventsNode.append_attribute( "minTime" ) = convertTime( m_begin );
  eventsNode.append_attribute( "maxTime" ) = convertTime( m_end );
}

void operator>>( const YAML::Node & node,
                 Simulation & simulation )
{
  simulation.setBegin( node["begin"].as< string >() );
  simulation.setEnd( node["end"].as< string >() );

  std::vector< std::shared_ptr< solvers::Solver > > solvers;
  // TODO check solver or solvers
  node["solver"] >> solvers;
  simulation.setSolver( solvers );

  std::shared_ptr< numericalStrategies::NumericalStrategies > ns;
  node["numerical_strategies"] >> ns;
  simulation.setNumericalStrategy( ns );
}

}
