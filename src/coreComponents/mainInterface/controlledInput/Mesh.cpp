#include "Mesh.hpp"
#include "Misc.hpp"

#include "common/DataTypes.hpp"

namespace geos::input::meshes
{

using namespace pugi;

class InternalMesh : public Mesh // TODO make VtkMesh inherit from Mesh
{
public:
  InternalMesh( string const & elementType,
                std::vector< string > const & xRange,
                std::vector< string > const & yRange,
                std::vector< string > const & zRange,
                std::vector< string > const & nx,
                std::vector< string > const & ny,
                std::vector< string > const & nz )
    : m_elementType( elementType ),
      m_xRange( xRange ),
      m_yRange( yRange ),
      m_zRange( zRange ),
      m_nx( nx ),
      m_ny( ny ),
      m_nz( nz )
  { }

  void fillMeshXmlNode( xml_node & meshNode ) const override
  {
    xml_node internal = meshNode.append_child( "InternalMesh" );
    internal.append_attribute( "name" ) = "__internal_mesh";
    internal.append_attribute( "elementTypes" ) = ( "{ " + convertYamlElementTypeToGeosElementType( m_elementType ) + " }" ).c_str();
    internal.append_attribute( "xCoords" ) = createGeosArray( m_xRange ).c_str();
    internal.append_attribute( "yCoords" ) = createGeosArray( m_yRange ).c_str();
    internal.append_attribute( "zCoords" ) = createGeosArray( m_zRange ).c_str();
    internal.append_attribute( "nx" ) = createGeosArray( m_nx ).c_str();
    internal.append_attribute( "ny" ) = createGeosArray( m_ny ).c_str();
    internal.append_attribute( "nz" ) = createGeosArray( m_nz ).c_str();
    internal.append_attribute( "cellBlockNames" ) = "{ cb1 }";  // TODO Improve the cell block mgmt!
  }

private:
  string m_elementType;
  std::vector< string > m_xRange;
  std::vector< string > m_yRange;
  std::vector< string > m_zRange;
  std::vector< string > m_nx;
  std::vector< string > m_ny;
  std::vector< string > m_nz;
};

void operator>>( const YAML::Node & node,
                 std::shared_ptr< Mesh > & mesh )
{
  if( node["internal"] )
  {
    const YAML::Node & internal = node["internal"];
    mesh = std::make_shared< InternalMesh >(
      internal["element_type"].as< string >(),
      internal["x_range"].as< std::vector< string > >(),
      internal["y_range"].as< std::vector< string > >(),
      internal["z_range"].as< std::vector< string > >(),
      internal["nx"].as< std::vector< string > >(),
      internal["ny"].as< std::vector< string > >(),
      internal["nz"].as< std::vector< string > >()
    );
  }
}


} // end of namespaces