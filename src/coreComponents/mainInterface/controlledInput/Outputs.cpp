#include "Outputs.hpp"

namespace geos::input::outputs
{

using namespace pugi;

class VtkOutput : public Output
{
public:
  VtkOutput( int counter,
             string const & every,
             std::vector< string > const & at,
             std::vector< string > const & fields,
             bool writeGhostCells,
             string fileRoot )
    : Output( "__vtk-" + std::to_string( counter ), every, at ),
      m_fields( fields ),
      m_writeGhostCells( writeGhostCells ),
      m_fileRoot( std::move( fileRoot ) )
  { }

  void fillOutputsXmlNode( xml_node & outputsNode ) const override
  {
    xml_node vtkOutput = outputsNode.append_child( "VTK" );
    vtkOutput.append_attribute( "name" ) = m_name.c_str();
    vtkOutput.append_attribute( "writeGhostCells" ) = m_writeGhostCells ? "1" : "0";
    vtkOutput.append_attribute( "fieldNames" ) = createGeosArray( m_fields ).c_str();
    if( !m_fileRoot.empty() )
    {
      vtkOutput.append_attribute( "plotFileRoot" ) = m_fileRoot.c_str();
    }
  }

private:
  std::vector< string > m_fields;
  bool m_writeGhostCells;
  string m_fileRoot;
};


class RestartOutput : public Output
{
public:
  RestartOutput( int counter,
                 string const & every,
                 std::vector< string > const & at )
    : Output( "__restart-" + std::to_string( counter ), every, at )
  { }

  void fillOutputsXmlNode( xml_node & outputsNode ) const override
  {
    xml_node vtkOutput = outputsNode.append_child( "Restart" );
    vtkOutput.append_attribute( "name" ) = m_name.c_str();
  }
};

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Output > > & outputs )
{
  GEOS_ASSERT( node.IsSequence() );

  for( std::size_t i = 0; i < node.size(); i++ )
  {
    auto output = node[i];
    for( auto const & kv: output )
    {
      string const outputType = kv.first.as< string >();
      auto const & subNode = kv.second;
      if( outputType == "vtk" )
      {
        auto sp = std::make_shared< VtkOutput >( i,
                                                 subNode["every"].as< string >( "" ),
                                                 subNode["at"].as< std::vector< string > >( std::vector< string >() ),
                                                 subNode["fields"].as< std::vector< string > >( std::vector< string >() ),
                                                 subNode["write_ghost_cells"].as< bool >( false ),
                                                 subNode["file_root"].as< string >( "vtkOutput" ) );
        outputs.push_back( sp );
      }
      else if( outputType == "restart" )
      {
        auto sp = std::make_shared< RestartOutput >( i,
                                                     subNode["every"].as< string >( "" ),
                                                     subNode["at"].as< std::vector< string > >( std::vector< string >() ) );
        outputs.push_back( sp );
      }
      else
      {
        GEOS_WARNING( "Discarded output \"" << outputType << "\"" );
      }
    }
  }
}

} // end of namespace