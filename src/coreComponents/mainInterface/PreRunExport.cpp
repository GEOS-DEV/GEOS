/* SPDX-License-Identifier: LGPL-2.1-only */
#include "PreRunExport.hpp"

#include "common/GeosxConfig.hpp"
#include "common/initializeEnvironment.hpp"
#include "common/MpiWrapper.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "functions/FunctionManager.hpp"
#ifdef GEOS_PRERUN_USE_VTK
#include "fileIO/vtk/VTKPolyDataWriterInterface.hpp"
#endif
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/version.hpp"

#include <cerrno>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

namespace geos
{
namespace
{
namespace fs = std::filesystem;
using dataRepository::Group;
using dataRepository::WrapperBase;
using xmlWrapper::xmlNode;

string quote( string const & value )
{
  std::ostringstream out;
  out << '"';
  for( unsigned char c : value )
  {
    switch( c )
    {
      case '"': out << "\\\""; break;
      case '\\': out << "\\\\"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if( c < 32 ) out << "\\u" << std::hex << std::setw( 4 ) << std::setfill( '0' ) << int( c ) << std::dec;
        else out << c;
    }
  }
  return out.str() + '"';
}

void writeFile( fs::path const & path, string const & text )
{
  std::ofstream file( path, std::ios::binary | std::ios::trunc );
  file << text;
  file.close();
  if( !file ) throw std::runtime_error( "Cannot write export file: " + path.string() );
}

string readFile( fs::path const & path )
{
  std::ifstream file( path, std::ios::binary );
  if( !file ) throw std::runtime_error( "Cannot read export file: " + path.string() );
  return { std::istreambuf_iterator< char >( file ), std::istreambuf_iterator< char >() };
}

// A deterministic schema identity, not a cryptographic authenticity claim.
string fingerprint( string const & bytes )
{
  uint64_t hash = UINT64_C( 14695981039346656037 );
  for( unsigned char c : bytes ) { hash ^= c; hash *= UINT64_C( 1099511628211 ); }
  std::ostringstream out;
  out << std::hex << std::setw( 16 ) << std::setfill( '0' ) << hash;
  return out.str();
}

// link() makes publication exclusive and atomic, including a racing existing destination.
void publishFile( fs::path const & staged, fs::path const & final )
{
  if( ::link( staged.c_str(), final.c_str() ) != 0 )
    throw std::runtime_error( "Cannot publish export (destination must be new): " + final.string() );
}

struct PropertyMetadata
{
  string description;
  string type;
  string inputFlag;
  string pattern;
  string limits;
  int limitsMode = 0;
};
using Metadata = std::map< string, PropertyMetadata >;

void collectMetadata( Group const & group, Metadata & metadata )
{
  for( auto const & entry : group.wrappers() )
  {
    WrapperBase const & wrapper = *entry.second;
    if( wrapper.getInputFlag() <= dataRepository::InputFlags::FALSE ) continue;
    string const key = group.getName() + "Type/" + wrapper.getName();
    PropertyMetadata value{ wrapper.getDescription(), wrapper.getRTTypeName(),
                            dataRepository::InputFlagToString( wrapper.getInputFlag() ),
                            wrapper.getTypeRegex().m_regexStr, wrapper.getLimitsString(),
                            static_cast< int >( wrapper.getLimitsMode() ) };
    auto const previous = metadata.find( key );
    // SchemaConstruction merges identically named types. Refuse to silently invent
    // one description/range if two registered wrappers disagree about that type.
    if( previous != metadata.end() &&
        ( previous->second.type != value.type || previous->second.limits != value.limits ||
          previous->second.pattern != value.pattern || previous->second.inputFlag != value.inputFlag ) )
      throw std::runtime_error( "Conflicting input metadata for schema property: " + key );
    metadata.emplace( key, std::move( value ) );
  }
  group.forSubGroups( [&]( Group const & child ) { collectMetadata( child, metadata ); } );
}

string numericBounds( string const & interval )
{
  static std::regex const pattern( R"(^([\[(])\s*([^,]+),\s*([^\])]+)\s*([\])])$)" );
  std::smatch match;
  if( !std::regex_match( interval, match, pattern ) ) return "";
  std::ostringstream out;
  auto add = [&]( string const & text, char const * key, bool exclusive )
  {
    try
    {
      size_t used = 0;
      double const value = std::stod( text, &used );
      if( text.find_first_not_of( " \t", used ) != string::npos || !std::isfinite( value ) ) return;
      out << ',' << quote( key ) << ':' << std::setprecision( 17 ) << value;
      out << ',' << quote( string( "exclusive" ) + ( string( key ) == "minimum" ? "Minimum" : "Maximum" ) )
          << ':' << ( exclusive ? "true" : "false" );
    }
    catch( std::exception const & ) {} // The original interval remains available verbatim.
  };
  add( match[2], "minimum", match[1] == "(" );
  add( match[3], "maximum", match[4] == ")" );
  return out.str();
}

string choices( string pattern )
{
  // Only a complete literal alternation is a choice list. Arbitrary regex syntax
  // is preserved as pattern, never mistaken for an enum.
  if( pattern.size() >= 2 && pattern.front() == '(' && pattern.back() == ')' )
    pattern = pattern.substr( 1, pattern.size() - 2 );
  if( !std::regex_match( pattern, std::regex( "[A-Za-z0-9_-]+(\\|[A-Za-z0-9_-]+)+" ) ) ) return "[]";
  std::ostringstream out;
  out << '[';
  std::istringstream in( pattern );
  string item;
  bool first = true;
  while( std::getline( in, item, '|' ) )
  {
    if( !first ) out << ',';
    first = false;
    out << quote( item );
  }
  return out.str() + ']';
}

void appendElement( std::ostream & out, xmlNode schema, xmlNode element,
                    string const & parent, Metadata const & metadata,
                    std::set< string > ancestors, bool & first )
{
  string const name = element.attribute( "name" ).value();
  string const typeName = element.attribute( "type" ).value();
  string const path = parent + '/' + name;
  xmlNode const type = schema.find_child_by_attribute( "xsd:complexType", "name", typeName.c_str() );
  if( !type ) throw std::runtime_error( "Missing schema type: " + typeName );
  if( !first ) out << ",\n";
  first = false;
  out << "{\"type\":" << quote( name ) << ",\"schemaType\":" << quote( typeName )
      << ",\"path\":" << quote( path ) << ",\"group\":" << quote( parent )
      << ",\"description\":\"\",\"properties\":[";
  bool firstProperty = true;
  for( xmlNode const & property : type.children( "xsd:attribute" ) )
  {
    string const propertyName = property.attribute( "name" ).value();
    auto const found = metadata.find( typeName + '/' + propertyName );
    PropertyMetadata value;
    value.type = property.attribute( "type" ).value();
    bool const required = string( property.attribute( "use" ).value() ) == "required";
    value.inputFlag = required ? "REQUIRED" : "OPTIONAL";
    if( found != metadata.end() ) value = found->second;
    else
    {
      xmlNode const comment = property.previous_sibling();
      string const prefix = propertyName + " => ";
      string const description = comment.value();
      if( comment.type() == xmlWrapper::xmlNodeType::node_comment && description.rfind( prefix, 0 ) == 0 )
        value.description = description.substr( prefix.size() );
    }
    if( value.pattern.empty() )
    {
      auto simple = schema.find_child_by_attribute( "xsd:simpleType", "name", property.attribute( "type" ).value() );
      value.pattern = simple.child( "xsd:restriction" ).child( "xsd:pattern" ).attribute( "value" ).value();
    }
    if( !firstProperty ) out << ',';
    firstProperty = false;
    out << "{\"name\":" << quote( propertyName ) << ",\"path\":" << quote( path + "/@" + propertyName )
        << ",\"type\":" << quote( value.type ) << ",\"schemaType\":" << quote( property.attribute( "type" ).value() )
        << ",\"description\":" << quote( value.description )
        << ",\"required\":" << ( required ? "true" : "false" )
        << ",\"inputFlag\":" << quote( value.inputFlag )
        << ",\"units\":null,\"unitsStatus\":\"unknown\",\"mutability\":\"unknown\""
        << ",\"pattern\":" << quote( value.pattern ) << ",\"choices\":" << choices( value.pattern )
        << ",\"limits\":" << quote( value.limits ) << ",\"limitsMode\":" << value.limitsMode
        << numericBounds( value.limits );
    if( property.attribute( "default" ) ) out << ",\"default\":" << quote( property.attribute( "default" ).value() );
    out << '}';
  }
  xmlNode const choice = type.child( "xsd:choice" );
  out << "],\"childChoice\":{\"minOccurs\":" << choice.attribute( "minOccurs" ).as_int( 1 )
      << ",\"maxOccurs\":" << quote( choice.attribute( "maxOccurs" ).as_string( "1" ) )
      << "},\"children\":[";
  bool firstChild = true;
  for( xmlNode const & child : type.child( "xsd:choice" ).children( "xsd:element" ) )
  {
    if( !firstChild ) out << ',';
    firstChild = false;
    out << "{\"name\":" << quote( child.attribute( "name" ).value() )
        << ",\"type\":" << quote( child.attribute( "type" ).value() )
        << ",\"minOccurs\":" << child.attribute( "minOccurs" ).as_int( 1 )
        << ",\"maxOccurs\":" << quote( child.attribute( "maxOccurs" ).as_string( "1" ) ) << '}';
  }
  bool const recursive = !ancestors.insert( typeName ).second;
  out << "],\"recursive\":" << ( recursive ? "true" : "false" ) << '}';
  if( recursive ) return; // Children describe the recursive schema edge without infinite expansion.
  for( xmlNode const & child : type.child( "xsd:choice" ).children( "xsd:element" ) )
    appendElement( out, schema, child, path, metadata, ancestors, first );
}

void writeCatalog( ProblemManager & problem, fs::path const & schemaPath, fs::path const & output )
{
  Metadata metadata;
  collectMetadata( problem, metadata );
  collectMetadata( problem.getFunctionManager(), metadata );
  collectMetadata( problem.getFieldSpecificationManager(), metadata );
  pugi::xml_document doc;
  if( !doc.load_file( schemaPath.c_str(), pugi::parse_default | pugi::parse_comments ) )
    throw std::runtime_error( "Cannot read generated schema" );
  xmlNode schema = doc.child( "xsd:schema" );
  std::ostringstream out;
  out << "{\"formatVersion\":1,\"schemaVersion\":1,\"geosVersion\":" << quote( getVersion() )
      << ",\"scope\":\"global\",\"mpiSize\":" << MpiWrapper::commSize() << ",\"schemaIdentity\":{\"algorithm\":\"fnv1a64\",\"digest\":"
      << quote( fingerprint( readFile( schemaPath ) ) )
      << "},\"unitMetadata\":\"unavailable-in-wrapper-registry\",\"elements\":[\n";
  bool first = true;
  appendElement( out, schema, schema.find_child_by_attribute( "xsd:element", "name", "Problem" ),
                 "", metadata, {}, first );
  out << "\n]}\n";
  writeFile( output, out.str() );
}

#ifdef GEOS_PRERUN_USE_VTK
void rewriteMeshReferences( xmlNode node, fs::path const & sourceRoot, string const & prefix )
{
  if( node.attribute( "file" ) )
  {
    fs::path const relative( node.attribute( "file" ).value() );
    if( relative.is_absolute() || relative.string().find( ".." ) != string::npos ||
        !fs::is_regular_file( sourceRoot / relative ) || fs::file_size( sourceRoot / relative ) == 0 )
      throw std::runtime_error( "Mesh export has a missing or unsafe sidecar: " + relative.string() );
    node.attribute( "file" ).set_value( ( prefix + '/' + relative.generic_string() ).c_str() );
  }
  for( xmlNode child : node.children() ) rewriteMeshReferences( child, sourceRoot, prefix );
}
#endif
}

bool handleCapabilitiesCommand( int argc, char * argv[], int & exitCode )
{
  bool requested = false;
  for( int i = 1; i < argc; ++i ) if( string( argv[i] ) == "--capabilities" ) requested = true;
  if( !requested ) return false;
  bool valid = argc == 2 || ( argc == 3 && string( argv[2] ) == "--format=json" ) ||
               ( argc == 4 && string( argv[2] ) == "--format" && string( argv[3] ) == "json" );
  valid = valid && string( argv[1] ) == "--capabilities";
  exitCode = valid ? 0 : 2;
  if( !valid )
  {
    std::cout << "{\"error\":{\"code\":\"invalid-capabilities-arguments\",\"message\":\"Use --capabilities --format=json as a standalone command\"}}\n";
    return true;
  }
#ifdef GEOS_PRERUN_USE_VTK
  char const * meshAvailable = "true";
#else
  char const * meshAvailable = "false";
#endif
  std::cout << "{\"formatVersion\":1,\"version\":" << quote( getVersion() )
            << ",\"schemaVersion\":1,\"input-catalog\":true,\"export-mesh\":" << meshAvailable << ","
            << "\"inputCatalogFormatVersion\":1,\"meshExportFormatVersion\":1,"
            << "\"meshExportFormat\":\"vtm\",\"inputCatalogScope\":\"global\","
            << "\"arbitrary-plane-phase-initialization\":false,\"generated-set-phase-initialization\":false,"
            << "\"well-neighborhood-refinement\":false}\n";
  return true;
}

void reportPreRunExportError( int argc, char * argv[], char const * message )
{
  for( int i = 1; i < argc; ++i )
  {
    string const argument( argv[i] );
    if( argument == "--input-catalog" || argument == "--export-mesh" ||
        argument.rfind( "--input-catalog=", 0 ) == 0 || argument.rfind( "--export-mesh=", 0 ) == 0 )
    {
      std::cout << "{\"error\":{\"code\":\"pre-run-export-failed\",\"message\":" << quote( message ) << "}}" << std::endl;
      return;
    }
  }
}

void runPreRunExport( std::unique_ptr< CommandLineOptions > options )
{
  bool const catalog = !options->inputCatalog.empty();
  fs::path const target = fs::absolute( catalog ? options->inputCatalog : options->exportMesh );
  fs::path const sidecars = target.string() + ".data";
  if( target.filename().empty() || !fs::is_directory( target.parent_path() ) || fs::exists( target ) ||
      ( !catalog && ( target.extension() != ".vtm" || fs::exists( sidecars ) ) ) )
    throw std::runtime_error( "Export requires a new output file in an existing directory (mesh suffix: .vtm)" );
  int owner = static_cast< int >( ::getpid() );
  MpiWrapper::bcast( &owner, 1, 0 );
  fs::path const stage = target.parent_path() / ( "." + target.filename().string() + ".partial-" + std::to_string( owner ) );
  if( MpiWrapper::commRank() == 0 && !fs::create_directory( stage ) )
    throw std::runtime_error( "Cannot create exclusive export staging directory" );
  MpiWrapper::barrier( MPI_COMM_GEOS );
  // Input preprocessing may write a combined XML file. Confine it, solver setup
  // diagnostics, and all generated intermediates to the unpublished staging tree.
  options->outputDirectory = stage.string();
  if( catalog )
  {
    if( !options->inputFileNames.empty() )
    {
      auto validationOptions = std::make_unique< CommandLineOptions >( *options );
      GeosxState validation( std::move( validationOptions ) );
      validation.initializeDataRepository();
      validation.applyInitialConditions();
    }
    options->inputFileNames.clear();
    fs::path const schema = stage / ( "schema-" + std::to_string( MpiWrapper::commRank() ) + ".xsd" );
    options->schemaName = schema.string();
    GeosxState state( std::move( options ) );
    state.initializeDataRepository();
    if( MpiWrapper::commRank() == 0 )
    {
      writeCatalog( state.getProblemManager(), schema, stage / "catalog.json" );
      publishFile( stage / "catalog.json", target );
    }
  }
  else
  {
#ifdef GEOS_PRERUN_USE_VTK
    GeosxState state( std::move( options ) );
    state.initializeDataRepository();
    vtk::VTKPolyDataWriterInterface writer( "mesh" );
    writer.setOutputLocation( stage.string(), "mesh" );
    writer.writeMesh( state.getProblemManager().getDomainPartition() );
    MpiWrapper::barrier( MPI_COMM_GEOS );
    if( MpiWrapper::commRank() == 0 )
    {
      xmlWrapper::xmlDocument doc;
      if( !doc.loadFile( ( stage / "mesh/000000.vtm" ).string() ) )
        throw std::runtime_error( "Mesh export did not write a VTM manifest" );
      rewriteMeshReferences( doc.getChild( "VTKFile" ), stage / "mesh", sidecars.filename().string() );
      if( !doc.saveFile( ( stage / "mesh.vtm" ).string() ) ) throw std::runtime_error( "Cannot write mesh manifest" );
      writeFile( stage / "mesh/metadata.json",
                 "{\"formatVersion\":1,\"kind\":\"geos-pre-run-mesh\",\"geosVersion\":" + quote( getVersion() ) +
                 ",\"coordinateUnits\":null,\"unitsStatus\":\"unknown\",\"ordering\":\"native-vtk-writer\",\"coordinatePrecision\":\"float64\","
                 "\"ghostCells\":true,\"timeLoopEntered\":false,\"initialConditionsApplied\":false}\n" );
      fs::remove( stage / "mesh/000000.vtm" );
      // Publish all completed sidecars first; the authoritative VTM is the commit marker.
      fs::rename( stage / "mesh", sidecars );
      publishFile( stage / "mesh.vtm", target );
    }
#else
    throw std::runtime_error( "--export-mesh is unavailable: this executable was built without VTK" );
#endif
  }
  MpiWrapper::barrier( MPI_COMM_GEOS );
  if( MpiWrapper::commRank() == 0 ) fs::remove_all( stage );
}
}
