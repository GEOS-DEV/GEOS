#include "Mesh.hpp"
#include "Misc.hpp"

#include <map>
#include <vector>

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
                std::vector< string > const & nz,
                std::map< string, std::vector< std::vector< int > > > const & regionToIndices )
    : m_elementType( elementType ),
      m_xRange( xRange ),
      m_yRange( yRange ),
      m_zRange( zRange ),
      m_nx( nx ),
      m_ny( ny ),
      m_nz( nz ),
      m_regionToIndices( regionToIndices )
  { }

  std::vector< string > getDomains() const override
  {
    std::vector< string > result;
    for( auto kv: m_regionToIndices )
    {
      result.emplace_back( kv.first );
    }
    return result.empty() ? std::vector< string >{ "__domain" } : result;
  }

  void fillProblemXmlNode( pugi::xml_node & problemNode ) const override
  {
    xml_node meshNode = problemNode.select_node( "Mesh" ).node();
    xml_node internal = meshNode.append_child( "InternalMesh" );
    internal.append_attribute( "name" ) = "__internal_mesh";
    internal.append_attribute( "elementTypes" ) = ( "{ " + convertYamlElementTypeToGeosElementType( m_elementType ) + " }" ).c_str();
    internal.append_attribute( "xCoords" ) = createGeosArray( m_xRange ).c_str();
    internal.append_attribute( "yCoords" ) = createGeosArray( m_yRange ).c_str();
    internal.append_attribute( "zCoords" ) = createGeosArray( m_zRange ).c_str();
    internal.append_attribute( "nx" ) = createGeosArray( m_nx ).c_str();
    internal.append_attribute( "ny" ) = createGeosArray( m_ny ).c_str();
    internal.append_attribute( "nz" ) = createGeosArray( m_nz ).c_str();

    if( m_regionToIndices.empty() )
    {
      fillProblemXmlNodeWithoutRegions( problemNode, internal );
    }
    else
    {
      fillProblemXmlNodeWithRegions( problemNode, internal );
    }
  }

private:

  void fillProblemXmlNodeWithoutRegions( pugi::xml_node & problemNode,
                                         pugi::xml_node & internal ) const
  {
    string const cbs = "{ __cb___single }";
    internal.append_attribute( "cellBlockNames" ) = cbs.c_str();
    xml_node er = problemNode.select_node( "ElementRegions" ).node();
    xml_node cer = er.append_child( "CellElementRegion" );
    cer.append_attribute( "name" ) = getDomains().front().c_str();
    cer.append_attribute( "cellBlocks" ) = cbs.c_str();
    cer.append_attribute( "materialList" ) = "{ nullModel }";
  }

  void fillProblemXmlNodeWithRegions( pugi::xml_node & problemNode,
                                      pugi::xml_node & internal ) const
  {
    std::vector< string > cellBlocks( m_nx.size() * m_ny.size() * m_nz.size(), "__cb___inactive" );
    std::map< string, string > regionToCellBlock;
    for( auto const & [regionName, indices]: m_regionToIndices )  // TODO manage empty case
    {
      for (std::vector< int > const & is: indices)
      {
        int const i = is[0], j = is[1], k = is[2];
        int const offset = k + m_nz.size() * j + m_nz.size() * m_ny.size() * i;
        cellBlocks[offset] = "__cb_" + regionName;
      }
      regionToCellBlock[regionName] = "__cb_" + regionName;
    }

    internal.append_attribute( "cellBlockNames" ) = createGeosArray( cellBlocks ).c_str();
    xml_node er = problemNode.select_node( "ElementRegions" ).node();
    for( auto const & [regionName, cellBlockName]: regionToCellBlock )
    {
      xml_node cer = er.append_child( "CellElementRegion" );
      cer.append_attribute( "name" ) = regionName.c_str();
      cer.append_attribute( "cellBlocks" ) = ( "{ " + cellBlockName + " }" ).c_str();
      cer.append_attribute( "materialList" ) = "{ nullModel }";
    }
    if( std::find( cellBlocks.cbegin(), cellBlocks.cend(), "__cb___inactive" ) != cellBlocks.cend() )
    {
      xml_node cer = er.append_child( "CellElementRegion" );
      cer.append_attribute( "name" ) = "__inactive_cer";
      cer.append_attribute( "cellBlocks" ) = "{ __cb___inactive }";
      cer.append_attribute( "materialList" ) = "{  }";
    }
  }

  string m_elementType;
  std::vector< string > m_xRange;
  std::vector< string > m_yRange;
  std::vector< string > m_zRange;
  std::vector< string > m_nx;
  std::vector< string > m_ny;
  std::vector< string > m_nz;
  std::map< string, std::vector< std::vector< int > > > const m_regionToIndices;
};

void operator>>( const YAML::Node & node,
                 std::shared_ptr< Mesh > & mesh )
{
  if( YAML::Node const & internal = node["internal"] )
  {
    std::map< string, std::vector< std::vector< int > > > regionToIndices;
    if( YAML::Node const & regions = internal["regions"] )
    {
      for( auto const & kv: regions["mapping"] )
      {
        string const regionName = kv.first.as< string >();
        std::vector< std::vector< int > > const indices = kv.second.as< std::vector< std::vector< int > > >();
        regionToIndices[regionName] = indices;
      }
    }

    mesh = std::make_shared< InternalMesh >(
      internal["element_type"].as< string >(),
      internal["x_range"].as< std::vector< string > >(),
      internal["y_range"].as< std::vector< string > >(),
      internal["z_range"].as< std::vector< string > >(),
      internal["nx"].as< std::vector< string > >(),
      internal["ny"].as< std::vector< string > >(),
      internal["nz"].as< std::vector< string > >(),
      regionToIndices
    );
  }
}


} // end of namespaces