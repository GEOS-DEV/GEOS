#include "NumericalStrategies.hpp"

#include "common/DataTypes.hpp"

namespace geos::input::numericalStrategies
{

using namespace pugi;

class NumericalStrategiesImpl : public NumericalStrategies
{
public:
  explicit NumericalStrategiesImpl( string const & type )
    : m_type( type )
  { }

  void fillProblemXmlNode( pugi::xml_node & problemNode ) const override
  {
    xml_node firstSolver = problemNode.select_node( "Solvers" ).node().first_child();
    // TODO assert there's one unique child in this precise case.
    firstSolver.append_attribute( "timeIntegrationOption" ) = m_timeIntegrationOption.c_str();
    firstSolver.append_child( "LinearSolverParameters" ).append_attribute( "directParallel" ).set_value( m_type == "serial_direct" ? 0 : 1 );
  }

private:
  string m_type;
  string m_timeIntegrationOption = "SteadyState";
};


void operator>>( const YAML::Node & node,
                 std::shared_ptr< NumericalStrategies > & numericalStrategy )
{
  string const linear_solve_type = node["linear_solver"]["type"].as< string >();

  numericalStrategy = std::make_shared< NumericalStrategiesImpl >( linear_solve_type );
}

} // end of namespace
