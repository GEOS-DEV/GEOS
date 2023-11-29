#include "ControlledInput.hpp"

#include "dataRepository/xmlWrapper.hpp"

#include "codingUtilities/StringUtilities.hpp"

#include <yaml-cpp/yaml.h>

#include <map>
#include <vector>

namespace geos::input {

using namespace pugi;

real64 convertTime(string const & time)
{
  std::vector< string > tokens = stringutilities::tokenize< std::vector >( time, " " );
  std::size_t const numTokens = tokens.size();
  GEOS_ASSERT( numTokens == 1 || numTokens == 2 );

  real64 const t = std::stod( tokens.front() ); // TODO check cast.

  real64 scaling = 1.;
  if( numTokens == 2 )
  {
    GEOS_ASSERT_EQ( tokens.size(), 2 );
    int constexpr sec = 1;
    int constexpr min = 60;
    int constexpr hour = 60 * min;
    int constexpr day = 24 * hour;
    int constexpr week = 7 * day;
    real64 constexpr year = 365.25 * day;
    real64 constexpr month = year / 12.;

    std::map< string, real64 > const conv{  // TODO To real64 because of the 365.25
      { "s",       sec },
      { "sec",     sec },
      { "second",  sec },
      { "second(s)",  sec },
      { "seconds", sec },
      { "min",     min },
      { "minute",  min },
      { "minute(s)",  min },
      { "minutes", min },
      { "h",       hour },
      { "hour",    hour },
      { "hour(s)",    hour },
      { "hours",   hour },
      { "d",       day },
      { "day",     day },
      { "day(s)",     day },
      { "days",    day },
      { "w",       week },
      { "week",    week },
      { "week(s)",    week },
      { "weeks",   week },
      { "m",       month },
      { "month",   month },
      { "month(s)",   month },
      { "months",  month },
      { "y",       year },
      { "year",    year },
      { "year(s)",    year },
      { "years",   year },
    };
    scaling = conv.at( tokens.back() );// TODO check if value is found
  }
  return t * scaling;
}

string convertYamlElementTypeToGeosElementType( string const yamlElementType )
{
  std::map< string, string > m{
    { "tetrahedra", "C3D4" },
    { "pyramids", "C3D5" },
    { "wedges", "C3D6" },
    { "hexahedra", "C3D8" },
    { "pentagonal_prism", "PentagonalPrism" },
    { "hexagonal_prism", "HexagonalPrism" },
    { "heptagonal_prism", "HeptagonalPrism" },
    { "octagonal_prism", "OctagonalPrism" },
    { "nonagonal_prism", "NonagonalPrism" },
    { "decagonal_prism", "DecagonalPrism" },
    { "hendecagonal_prism", "HendecagonalPrism" },
    { "polyhedron", "Polyhedron" },
  };

  auto const geosIt = m.find( yamlElementType );
  if( geosIt == m.cend() )
  {
    GEOS_ERROR( "Could not find element type " << yamlElementType << " in supported element list." );
  }
  return geosIt->second;
}

template< class T >
string createGeosArray( T const & t )
{
  return "{ " + stringutilities::join( t, ", " ) + " }";
}

class Simulation
{
public:

  void setBegin( string const & begin )
  {
    m_begin = begin;
  }

  void setEnd( string const & end )
  {
    m_end = end;
  }

  void fillEventsXmlNode( xml_node & eventsNode ) const
  {
    eventsNode.append_attribute( "minTime" ) = convertTime( m_begin );
    eventsNode.append_attribute( "maxTime" ) = convertTime( m_end );
  }

private:
  string m_begin;
  string m_end;
};

class Event
{
public:
  explicit Event( string const & target )
    : m_target( target )
  { }

  virtual ~Event() = default;

  virtual void fillEventsXmlNode( xml_node & eventsNode ) const = 0;
protected:
  string m_target;
};

class PeriodicEvent : public Event
{
public:
  PeriodicEvent( string const & target, string const & every )
    : Event( target ),
      m_every( every )
  { }

  void fillEventsXmlNode( xml_node & eventsNode ) const override
  {
    xml_node periodicEvent = eventsNode.append_child( "PeriodicEvent" );
    periodicEvent.append_attribute( "timeFrequency" ) = convertTime( m_every );
    periodicEvent.append_attribute( "target" ) = m_target.c_str();
  }

private:
  string m_every;
};

class SoloEvent : public Event
{
public:
  SoloEvent( string const & target, string const & at )
    : Event( target ),
      m_at( at )
  { }

  void fillEventsXmlNode( xml_node & eventsNode ) const override
  {
    xml_node soloEvent = eventsNode.append_child( "SoloEvent" );
    soloEvent.append_attribute( "targetTime" ) = convertTime( m_at );
    soloEvent.append_attribute( "target" ) = m_target.c_str();
  }

private:
  string m_at;
};

class Output
{
public:
  Output( string name,
          string const & every,
          std::vector< string > const & at )
    : m_name( name ),
      m_every( every ),
      m_at( at )
  { }

  virtual ~Output() = default;

  virtual void fillOutputsXmlNode( xml_node & outputsNode ) const = 0;

  std::vector< std::shared_ptr< Event > > getEvents() const
  {
    std::vector< std::shared_ptr< Event > > result;

    if( !m_every.empty() )
    {
      result.emplace_back( std::make_shared< PeriodicEvent >( "/Outputs/" + m_name, m_every ) );
    }

    for( string const & at: m_at )
    {
      if( !at.empty() )
      {
        result.emplace_back( std::make_shared< SoloEvent >( "/Outputs/" + m_name, at ) );
      }
    }

    return result;
  };

protected:
  string m_name;

private:
  string m_every;
  std::vector< string > m_at;
};

class VtkOutput: public Output
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


class RestartOutput: public Output
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

class Mesh
{
public:
  virtual ~Mesh() = default;
  virtual void fillMeshXmlNode( xml_node & meshNode ) const = 0;
};

class InternalMesh: public Mesh // TODO make VtkMesh inherit from Mesh
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

class Deck
{
public:

  void setSimulation( Simulation const & simulation )
  {
    m_simulation = simulation;
  }

  void setOutputs( std::vector< std::shared_ptr< Output > > const & outputs )
  {
    m_outputs = outputs;
  }

  void setMesh( std::shared_ptr< Mesh > const & mesh )
  {
    m_mesh = mesh;
  }

  void fillProblemXmlNode( xml_node & problemNode ) const
  {
    xml_node xmlOutputs = problemNode.append_child( "Outputs" );
    xml_node xmlEvents = problemNode.append_child( "Events" );
    xml_node xmlMesh = problemNode.append_child( "Mesh" );

    m_simulation.fillEventsXmlNode( xmlEvents );

    for( std::shared_ptr< Output > output: m_outputs )
    {
      output->fillOutputsXmlNode( xmlOutputs );
      for( std::shared_ptr< Event > event: output->getEvents() )
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
  std::vector< std::shared_ptr< Output > > m_outputs;
  std::shared_ptr< Mesh > m_mesh;
};

void operator>>( const YAML::Node & node,
                 Simulation & simulation )
{
  simulation.setBegin( node["begin"].as< string >() );
  simulation.setEnd( node["end"].as< string >() );
}

void operator>>( const YAML::Node & node,
                 std::vector< std::shared_ptr< Output > > & outputs )
{
  GEOS_ASSERT( node.IsSequence() );

  for( std::size_t i = 0; i < node.size(); i++ )
  {
    auto output = node[i];
    for (auto const & kv: output)
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

void operator>>( const YAML::Node & node,
                 Deck & deck )
{
  Simulation simulation;
  node["simulation"] >> simulation;
  deck.setSimulation( simulation );

  std::vector< std::shared_ptr< Output > > outputs;
  node["outputs"] >> outputs;
  deck.setOutputs( outputs );

  std::shared_ptr< Mesh > mesh;
  node["mesh"] >> mesh;
  deck.setMesh( mesh );
}

void fillWithMissingXmlInfo( xml_node & problem )
{
  xml_node laplaceFem = problem.append_child( "Solvers" ).append_child( "LaplaceFEM" );
  laplaceFem.append_attribute( "name" ) = "laplace";
  laplaceFem.append_attribute( "discretization" ) = "FE1";
  laplaceFem.append_attribute( "timeIntegrationOption" ) = "SteadyState";
  laplaceFem.append_attribute( "fieldName" ) = "Temperature";
  laplaceFem.append_attribute( "targetRegions" ) = "{ Domain }";
  laplaceFem.append_child( "LinearSolverParameters" ).append_attribute( "directParallel" ).set_value( 0 );

  xpath_node events = problem.select_node( "Events" );
  {
    xml_node pe = events.node().append_child( "PeriodicEvent" );
    pe.append_attribute( "name" ) = "solverApplications";
    pe.append_attribute( "forceDt" ) = "1.0";
    pe.append_attribute( "target" ) = "/Solvers/laplace";
  }

  problem.append_child( "Constitutive" ).append_child( "NullModel" ).append_attribute( "name" ).set_value( "nullModel" );

  xml_node fieldSpecifications = problem.append_child( "FieldSpecifications" );
  {
    xml_node fs = fieldSpecifications.append_child( "FieldSpecification" );
    fs.append_attribute( "name" ) = "sourceTerm";
    fs.append_attribute( "fieldName" ) = "Temperature";
    fs.append_attribute( "objectPath" ) = "nodeManager";
    fs.append_attribute( "functionName" ) = "DirichletTimeFunction";
    fs.append_attribute( "scale" ) = "1.0";
    fs.append_attribute( "setNames" ) = "{ source }";
  }
  {
    xml_node fs = fieldSpecifications.append_child( "FieldSpecification" );
    fs.append_attribute( "name" ) = "sinkTerm";
    fs.append_attribute( "fieldName" ) = "Temperature";
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

  xml_node fes = problem.append_child( "NumericalMethods" ).append_child( "FiniteElements" ).append_child( "FiniteElementSpace" );
  fes.append_attribute( "name" ) = "FE1";
  fes.append_attribute( "order" ) = 1;
}

void convert( string const & stableInputFileName, xmlWrapper::xmlDocument & doc )
{
//  xml_document & pugiDoc = doc.getPugiDocument();
//  pugiDoc.select_node( "/Problem/Events" );
  xml_node problem = doc.appendChild( "Problem" );

  YAML::Node const input = YAML::LoadFile( stableInputFileName );
  Deck deck;
  input >> deck;

  deck.fillProblemXmlNode( problem );

  fillWithMissingXmlInfo(problem );

  doc.getPugiDocument().save( std::cout, "    ", pugi::format_indent | pugi::format_indent_attributes );
}

}
