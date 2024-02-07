#include "Geometry.hpp"

#include "Misc.hpp"

#include "common/DataTypes.hpp"

#include <memory>

namespace geos::input::geometry
{

using namespace pugi;

class Box : public Geometry
{
public:
  Box( string const & name,
       std::vector< real64 > const & min,
       std::vector< real64 > const & max )
    : m_name( name ),
      m_min( min ),
      m_max( max )
  { }

  void fillProblemXmlNode( pugi::xml_node & problemNode ) const override
  {
    xml_node geometry = problemNode.select_node( "Geometry" ).node();
    xml_node box = geometry.append_child( "Box" );
    box.append_attribute( "name" ) = m_name.c_str();
    box.append_attribute( "xMin" ) = createGeosArray( m_min ).c_str();
    box.append_attribute( "xMax" ) = createGeosArray( m_max ).c_str();
  }

private:
  string m_name;
  std::vector< real64 > m_min;
  std::vector< real64 > m_max;
};


void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Geometry > > & geometries )
{
  GEOS_ASSERT( node.IsSequence() );

  for( std::size_t i = 0; i < node.size(); i++ )
  {
    for( auto const & kv: node[i] )
    {
      string const type = kv.first.as< string >();
      auto const & subNode = kv.second;
      string const name = subNode["name"].as< string >();
      if( type == "box" )
      {
        auto box = std::make_shared< Box >( name,
                                            subNode["min"].as< std::vector< real64 > >(),
                                            subNode["max"].as< std::vector< real64 > >() );
        geometries.push_back( box );
      }
      else
      {
        GEOS_WARNING( "Discarded geometry \"" << type << "\"" );
      }
    }
  }

}

}
