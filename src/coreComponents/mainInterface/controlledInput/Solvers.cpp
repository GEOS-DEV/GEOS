#include "Solvers.hpp"

namespace geos::input::solvers
{

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
  explicit Laplace( NamedVariable const & variable )
    : m_variable( variable )
  { }

private:
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
      subNode["formulation"]["variable"] >> nv;
      solvers.push_back( std::make_shared< Laplace >( nv ) );
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
