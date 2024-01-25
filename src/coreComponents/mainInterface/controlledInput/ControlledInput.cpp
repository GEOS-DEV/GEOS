#include "ControlledInput.hpp"

#include "Events.hpp"
#include "Mesh.hpp"
#include "Outputs.hpp"
#include "Solvers.hpp"
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

  void fillProblemXmlNode( xml_node & problemNode ) const
  {
    problemNode.append_child( "Constitutive" );
    xml_node xmlEvents = problemNode.append_child( "Events" );
    xml_node xmlMesh = problemNode.append_child( "Mesh" );
    problemNode.append_child( "NumericalMethods" );
    xml_node xmlOutputs = problemNode.append_child( "Outputs" );
    problemNode.append_child( "Solvers" );

    m_simulation.fillProblemXmlNode( problemNode );

    for( std::shared_ptr< outputs::Output > output: m_outputs )
    {
      output->fillOutputsXmlNode( xmlOutputs );
      for( std::shared_ptr< events::Event > event: output->getEvents() )
      {
        event->fillEventsXmlNode( xmlEvents );
      }
    }

    // Create and populate the mesh node
    m_mesh->fillMeshXmlNode( xmlMesh );

    // Add name to all the events.
    int iEvent = 0;
    for( auto it = xmlEvents.children().begin(); it != xmlEvents.children().end(); ++it, ++iEvent )
    {
      it->append_attribute( "name" ) = ( "__event-" + std::to_string( iEvent ) ).c_str();
    }
  }

private:
  Simulation m_simulation;
  std::vector< std::shared_ptr< outputs::Output > > m_outputs;
  std::shared_ptr< meshes::Mesh > m_mesh;
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
  node["mesh"] >> mesh;
  deck.setMesh( mesh );
}

void fillWithMissingXmlInfo( xml_node & problem )
{
  xml_node laplaceFem = problem.select_node( "Solvers/LaplaceFEM" ).node();
  laplaceFem.append_attribute( "timeIntegrationOption" ) = "SteadyState";
  laplaceFem.append_child( "LinearSolverParameters" ).append_attribute( "directParallel" ).set_value( 0 );

  xml_node fieldSpecifications = problem.append_child( "FieldSpecifications" );
  {
    xml_node fs = fieldSpecifications.append_child( "FieldSpecification" );
    fs.append_attribute( "name" ) = "sourceTerm";
    fs.append_attribute( "fieldName" ) = "temperature";
    fs.append_attribute( "objectPath" ) = "nodeManager";
    fs.append_attribute( "functionName" ) = "DirichletTimeFunction";
    fs.append_attribute( "scale" ) = "1.0";
    fs.append_attribute( "setNames" ) = "{ source }";
  }
  {
    xml_node fs = fieldSpecifications.append_child( "FieldSpecification" );
    fs.append_attribute( "name" ) = "sinkTerm";
    fs.append_attribute( "fieldName" ) = "temperature";
    fs.append_attribute( "objectPath" ) = "nodeManager";
    fs.append_attribute( "scale" ) = "0.0";
    fs.append_attribute( "setNames" ) = "{ sink }";
  }

  xml_node geometry = problem.append_child( "Geometry" );
  {
    xml_node box = geometry.append_child( "Box" );
    box.append_attribute( "name" ) = "source";
    box.append_attribute( "xMin" ) = "{ -0.01, -0.01, -0.01 }";
    box.append_attribute( "xMax" ) = "{ +0.01, +1.01, +1.01 }";
  }
  {
    xml_node box = geometry.append_child( "Box" );
    box.append_attribute( "name" ) = "sink";
    box.append_attribute( "xMin" ) = "{ +0.99, -0.01, -0.01 }";
    box.append_attribute( "xMax" ) = "{ +1.01, +1.01, +1.01 }";
  }

  xml_node tableFunctions = problem.append_child( "Functions" ).append_child( "TableFunction" );
  tableFunctions.append_attribute( "name" ) = "DirichletTimeFunction";
  tableFunctions.append_attribute( "inputVarNames" ) = "{ time }";
  tableFunctions.append_attribute( "coordinates" ) = "{ 0.0, 1.0, 2.0 }";
  tableFunctions.append_attribute( "values" ) = "{ 0.0, 3.e2, 4.e3 }";

  xml_node cesr = problem.append_child( "ElementRegions" ).append_child( "CellElementRegion" );
  cesr.append_attribute( "name" ) = "Domain";
  cesr.append_attribute( "cellBlocks" ) = "{ cb1 }";
  cesr.append_attribute( "materialList" ) = "{ nullModel }";
}

void convert( string const & stableInputFileName,
              xmlWrapper::xmlDocument & doc )
{
//  xml_document & pugiDoc = doc.getPugiDocument();
//  pugiDoc.select_node( "/Problem/Events" );
  xml_node problem = doc.appendChild( "Problem" );

  YAML::Node const input = YAML::LoadFile( stableInputFileName );
  Deck deck;
  input >> deck;

  deck.fillProblemXmlNode( problem );

  fillWithMissingXmlInfo( problem );

  doc.getPugiDocument().save( std::cout, "    ", pugi::format_indent | pugi::format_indent_attributes );
}

}
