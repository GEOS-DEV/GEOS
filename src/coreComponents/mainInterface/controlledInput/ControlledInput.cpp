#include "ControlledInput.hpp"

#include "BoundaryConditions.hpp"
#include "Events.hpp"
#include "Functions.hpp"
#include "Geometry.hpp"
#include "Mesh.hpp"
#include "Outputs.hpp"
#include "Simulation.hpp"

#include "dataRepository/xmlWrapper.hpp"

#include <yaml-cpp/yaml.h>

#include <vector>

namespace geos::input
{

using namespace pugi;

class Deck
{
public:

  void setSimulation( Simulation const & simulation )
  {
    m_simulation = simulation;
  }

  void setOutputs( std::vector< std::shared_ptr< outputs::Output > > const & outputs )
  {
    m_outputs = outputs;
  }

  void setMesh( std::shared_ptr< meshes::Mesh > const & mesh )
  {
    m_mesh = mesh;
  }

  void setFunctions( std::vector< std::shared_ptr< functions::Function > > const & functions )
  {
    m_functions = functions;
  }

  void setBoundaryConditions( std::vector< std::shared_ptr< bc::BoundaryConditions > > const & bcs )
  {
    m_bcs = bcs;
  }

  void setGeometries( std::vector< std::shared_ptr< geometry::Geometry > > const & geometries )
  {
    m_geometries = geometries;
  }

  void fillProblemXmlNode( xml_node & problemNode ) const
  {
    problemNode.append_child( "Constitutive" );
    xml_node xmlEvents = problemNode.append_child( "Events" );
    problemNode.append_child( "Mesh" );
    problemNode.append_child( "NumericalMethods" );
    xml_node xmlOutputs = problemNode.append_child( "Outputs" );
    problemNode.append_child( "Solvers" );
    xml_node fieldSpecifications = problemNode.append_child( "FieldSpecifications" );
    problemNode.append_child( "ElementRegions" );

    m_simulation.fillProblemXmlNode( problemNode, m_mesh->getDomains() );

    for( std::shared_ptr< outputs::Output > output: m_outputs )
    {
      output->fillOutputsXmlNode( xmlOutputs );
      for( std::shared_ptr< events::Event > event: output->getEvents() )
      {
        event->fillEventsXmlNode( xmlEvents );
      }
    }

    m_mesh->fillProblemXmlNode( problemNode );

    // Add name to all the events.
    int iEvent = 0;
    for( auto it = xmlEvents.children().begin(); it != xmlEvents.children().end(); ++it, ++iEvent )
    {
      it->append_attribute( "name" ) = ( "__event-" + std::to_string( iEvent ) ).c_str();
    }

    if( !m_functions.empty() )
    {
      problemNode.append_child( "Functions" );
      for( std::shared_ptr< functions::Function > function: m_functions )
      {
        function->fillProblemXmlNode( problemNode );
      }
    }

    for( std::shared_ptr< bc::BoundaryConditions > bc: m_bcs )
    {
      bc->fillProblemXmlNode( problemNode );
    }

    // Add name to all the events.
    int iBC = 0;
    for( auto it = fieldSpecifications.children().begin(); it != fieldSpecifications.children().end(); ++it, ++iBC )
    {
      it->append_attribute( "name" ) = ( "__fs-" + std::to_string( iBC ) ).c_str();
    }

    if( !m_geometries.empty() )
    {
      problemNode.append_child( "Geometry" );
      for( std::shared_ptr< geometry::Geometry > geometry: m_geometries )
      {
        geometry->fillProblemXmlNode( problemNode );
      }
    }
  }

private:
  Simulation m_simulation;
  std::vector< std::shared_ptr< functions::Function > > m_functions;
  std::vector< std::shared_ptr< bc::BoundaryConditions > > m_bcs;
  std::vector< std::shared_ptr< outputs::Output > > m_outputs;
  std::shared_ptr< meshes::Mesh > m_mesh;
  std::vector< std::shared_ptr< geometry::Geometry > > m_geometries;
};


void operator>>( const YAML::Node & node,
                 Deck & deck )
{
  Simulation simulation;
  node["simulation"] >> simulation;
  deck.setSimulation( simulation );

  std::vector< std::shared_ptr< outputs::Output > > outputs;
  node["outputs"] >> outputs;
  deck.setOutputs( outputs );

  std::shared_ptr< meshes::Mesh > mesh;
  auto yamlMesh = node["mesh"];
  yamlMesh >> mesh;
  deck.setMesh( mesh );

  std::vector< std::shared_ptr< functions::Function > > functions;
  node["functions"] >> functions;
  deck.setFunctions( functions );

  std::vector< std::shared_ptr< bc::BoundaryConditions > > bcs;
  node["boundary_conditions"] >> bcs;
  deck.setBoundaryConditions( bcs );

  std::vector< std::shared_ptr< geometry::Geometry > > geometries;
  if( auto yamlSurfaces = yamlMesh["surfaces"] )
  {
    yamlSurfaces >> geometries;
    deck.setGeometries( geometries );
  }
}

void fillWithMissingXmlInfo( xml_node & problem )
{
  problem.select_node("Constitutive").node().append_child( "NullModel" ).append_attribute( "name" ).set_value( "nullModel" ); // TODO null model hard coded
}

void checkInputFileVersion( YAML::Node const & input )
{
  string const version = input["geos"]["input_file_version"].as< string >();
  GEOS_ERROR_IF( version != "0.0.1", "Not supported input file version " << version );
}

void convert( string const & stableInputFileName,
              xmlWrapper::xmlDocument & doc )
{
  xml_node problem = doc.appendChild( "Problem" );

  YAML::Node const input = YAML::LoadFile( stableInputFileName );
  checkInputFileVersion( input );
  Deck deck;
  input >> deck;

  deck.fillProblemXmlNode( problem );

  fillWithMissingXmlInfo( problem );

  doc.getPugiDocument().save( std::cout, "    ", pugi::format_indent | pugi::format_indent_attributes );
}

}
