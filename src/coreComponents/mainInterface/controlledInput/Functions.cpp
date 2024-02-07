#include "Functions.hpp"
#include "Misc.hpp"

#include "common/DataTypes.hpp"

#include <vector>

namespace geos::input::functions
{

using namespace pugi;

class TableFunction : public Function
{
public:
  TableFunction( string const & name,
                 std::vector< real64 > const & xs,
                 std::vector< real64 > const & ys )
    : m_name(name),
      m_xs( xs ),
      m_ys( ys )
  { }

  void fillProblemXmlNode( xml_node & problemNode ) const override
  {
    xml_node functions = problemNode.select_node( "Functions" ).node();
    xml_node table = functions.append_child( "TableFunction" );
    table.append_attribute( "name" ) = m_name.c_str();
    table.append_attribute( "inputVarNames" ) = "{ time }";
    table.append_attribute( "coordinates" ) = createGeosArray( m_xs ).c_str();
    table.append_attribute( "values" ) = createGeosArray( m_ys ).c_str();
  }
private:
  string m_name;
  std::vector< real64 > m_xs;
  std::vector< real64 > m_ys;
};


void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Function > > & functions )
{
  GEOS_ASSERT( node.IsSequence() );

  for( std::size_t i = 0; i < node.size(); i++ )
  {
    auto function = node[i];
    for( auto const & kv: function )
    {
      string const type = kv.first.as< string >();
      auto const & subNode = kv.second;
      string const name = subNode["name"].as< string >();
      if( type == "table" )
      {
        auto f = std::make_shared< TableFunction >( name,
                                                    subNode["x"].as< std::vector< real64 > >(),
                                                    subNode["y"].as< std::vector< real64 > >() );
        functions.push_back( f );
      }
      else
      {
        GEOS_WARNING( "Discarded function \"" << type << "\"" );
      }
    }
  }

}

}
