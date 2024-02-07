#include "BoundaryConditions.hpp"

#include "common/DataTypes.hpp"


namespace geos::input::bc
{

using namespace pugi;

class Uniform: public BoundaryConditions
{
public:
  Uniform( string const & location,
           string const & field,
           real64 value )
    : m_location( location ),
      m_field( field ),
      m_value( value )
  { }

  void fillProblemXmlNode( xml_node & problemNode ) const override
  {
    xml_node functions = problemNode.select_node( "FieldSpecifications" ).node();
    xml_node fs = functions.append_child( "FieldSpecification" );
    fs.append_attribute( "fieldName" ) = m_field.c_str();
    fs.append_attribute( "objectPath" ) = "nodeManager";  // TODO
    fs.append_attribute( "scale" ) = std::to_string( m_value ).c_str();
    fs.append_attribute( "setNames" ) = ( "{ " + m_location + " }" ).c_str();
  }

private:
  string m_location;
  string m_field;
  real64 m_value;
};

class DirichletFct : public BoundaryConditions
{
public:
  DirichletFct( string const & location,
                string const & field,
                string const & functionName )
    : m_location( location ),
      m_field( field ),
      m_functionName( functionName )
  { }

  void fillProblemXmlNode( xml_node & problemNode ) const override
  {
    xml_node functions = problemNode.select_node( "FieldSpecifications" ).node();
    xml_node fs = functions.append_child( "FieldSpecification" );
    fs.append_attribute( "fieldName" ) = m_field.c_str();
    fs.append_attribute( "objectPath" ) = "nodeManager";  // TODO
    fs.append_attribute( "functionName" ) = m_functionName.c_str();
    fs.append_attribute( "scale" ) = "1.0";
    fs.append_attribute( "setNames" ) = ( "{ " + m_location + " }" ).c_str();
  }

private:
  string m_location;
  string m_field;
  string m_functionName;
};

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< BoundaryConditions > > & bcs )
{
  GEOS_ASSERT( node.IsSequence() );

  for( std::size_t i = 0; i < node.size(); i++ )
  {
    auto bc = node[i];
    for( auto const & kv: bc )
    {
      string const location = kv.first.as< string >();
      auto const & fieldSpecifications = kv.second;
      GEOS_ASSERT( fieldSpecifications.IsSequence() );
      for( std::size_t j = 0; j < fieldSpecifications.size(); j++ )
      {
        auto fs = fieldSpecifications[j];
        string const field = fs["field"].as< string >();
        if( auto val = fs["value"] )
        {
          auto uniform = std::make_shared< Uniform >( location, field, val.as< real64 >() );
          bcs.push_back(uniform);
        }
        else if( auto fct = fs["function"] )
        {
          auto dirichletFct = std::make_shared< DirichletFct >( location, field, fct.as< string >() );
          bcs.push_back(dirichletFct);
        }
        else
        {
          GEOS_WARNING( "Discarded bc \"" << fs << "\"" );
        }
      }
    }
  }
}

}