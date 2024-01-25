#include "Events.hpp"

#include "Solvers.hpp"

#include "Misc.hpp"

#include "common/DataTypes.hpp"


namespace geos::input::solvers
{

using namespace pugi;

struct Variable
{
public:
  string functional_space;  // TOD Use an enum here?
  int fe_order;
};

struct NamedVariable: public Variable
{
public:
  string name;
};

class Laplace : public Solver
{
public:
  explicit Laplace( NamedVariable const & variable,
                    std::vector< string > const & on ) // TODO move the `on` upper in the class hierarchy
    : m_on( on ),
      m_variable( variable )
  { }

  void fillProblemXmlNode( xml_node & problemNode ) const override
  {
    xml_node solvers = problemNode.select_node( "Solvers" ).node();
    xml_node laplaceFem = solvers.append_child( "LaplaceFEM" );
    string const solverName = "__laplace";
    laplaceFem.append_attribute( "name" ) = solverName.c_str();
    string const feSpaceName = "FE" + std::to_string( m_variable.fe_order );
    laplaceFem.append_attribute( "discretization" ) = feSpaceName.c_str();
    laplaceFem.append_attribute( "fieldName" ) = m_variable.name.c_str();
    laplaceFem.append_attribute( "targetRegions" ) = createGeosArray( m_on ).c_str();

    xml_node events = problemNode.select_node( "Events" ).node();
    events::PeriodicEvent pe = events::PeriodicEvent( "/Solvers/" + solverName );
    pe.fillEventsXmlNode( events );

    xml_node constitutive = problemNode.select_node( "Constitutive" ).node();
    constitutive.append_child( "NullModel" ).append_attribute( "name" ).set_value( "nullModel" );

    xml_node numericalMethods = problemNode.select_node( "NumericalMethods" ).node();
    xml_node feSpace = numericalMethods.append_child( "FiniteElements" ).append_child( "FiniteElementSpace" );
    feSpace.append_attribute( "name" ) = feSpaceName.c_str();
    feSpace.append_attribute( "order" ) = m_variable.fe_order;
  }

private:
  std::vector< string > m_on;
  NamedVariable m_variable;
};

void operator>>( const YAML::Node & node,
                 Variable & variable )
{
  variable.functional_space = node["functional_space"].as< string >();
  variable.fe_order = node["fe_order"].as< int >();
}

void operator>>( const YAML::Node & node,
                 NamedVariable & variable )
{
  node >> static_cast< Variable & >(variable);
  variable.name = node["name"].as< string >();
}

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Solver > > & solvers )
{

  for( auto const & kv: node )
  {
    string const solverType = kv.first.as< string >();
    auto const & subNode = kv.second;
    if( solverType == "laplace" )
    {
      NamedVariable nv;
      std::vector< string > on = subNode["on"].as< std::vector< string >>( std::vector< string >{ "Domain" } );
      subNode["formulation"]["variable"] >> nv;

      solvers.push_back( std::make_shared< Laplace >( nv, on ) );
      std::cout << nv.name << std::endl;
      std::cout << nv.functional_space << std::endl;
    }
    else
    {
      GEOS_WARNING( "Discarded solver \"" << solverType << "\"" );
    }
  }
}

} // end of namespace
