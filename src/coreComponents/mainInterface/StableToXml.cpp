/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "StableToXml.hpp"

#include <yaml-cpp/yaml.h>

#include <ostream>

namespace geosx
{
namespace api
{

// TODO Duplicated
template< typename IT, typename S = char >
string join( IT first, IT last, S const & delim = S() )
{
  if( first == last )
  {
    return {};
  }
  std::ostringstream oss;
  oss << *first;
  while( ++first != last )
  {
    oss << delim << *first;
  }
  return oss.str();
}

// TODO Duplicated
template< typename CONTAINER, typename S = char >
string join( CONTAINER const & cont, S const & delim = S() )
{
  return join( std::begin( cont ), std::end( cont ), delim );
}

// TODO Duplicated
std::vector< string > tokenize( string const & str,
                                string const & delimiters )
{
  if( str.empty() )
  {
    return {};
  }

  auto const isNonSpace = []( char const c ){ return !isspace( c ); };
  bool const usesNonWhitespaceDelimiters = std::any_of( delimiters.begin(), delimiters.end(), isNonSpace );

  // When only whitespace delimiters, skip multiple adjacent delimiters; otherwise don't and keep empty tokens
  std::vector<string> tokens;
  size_t lastPos = usesNonWhitespaceDelimiters ? 0 : str.find_first_not_of( delimiters, 0 );
  size_t newPos;
  while( ( newPos = str.find_first_of( delimiters, lastPos ) ) != string::npos )
  {
    tokens.emplace_back( str.substr( lastPos, newPos - lastPos ) );
    lastPos = usesNonWhitespaceDelimiters ? newPos + 1 : str.find_first_not_of( delimiters, newPos );
  }
  if( lastPos != string::npos )
  {
    tokens.emplace_back( str.substr( lastPos ) );
  }

  return tokens;
}



using namespace pugi;

string serialize( YAML::Node const & value )
{
  if( value.IsSequence() )
  {
    auto tmp = value.as< std::vector< string > >();
    return "{ " + join( tmp, ", " ) + " }";
  }
  else if( value.as< string >() == "true" )
  { return "1"; }
  else if( value.as< string >() == "false" )
  { return "0"; }
  else
  {
    return value.as< string >();
  }
}

string get( YAML::Node const & input,
            string key,
            string defaultValue )
{
  return input[key] ? serialize( input[key] ) : defaultValue;
}

void BuildMeshNode( YAML::Node const & input, pugi::xml_node & problem )
{
  int counter = 0;

  xml_node xmlMesh = problem.append_child( "Mesh" );

  std::map< string, string > const conv = { { "generated", "InternalMesh" },
                                            { "vtk",       "VTKMesh" },
                                            { "pamela",    "PAMELAMesh" } };

  // Let's loop on all the meshes to find the GEOSX mesh.
  for( YAML::Node const & mesh: input["mesh"] )
  {
    if( conv.count( mesh["type"].as< string >() ) == 0 ) continue;

    string const geosxMesh = conv.at( mesh["type"].as< string >() );

    if( geosxMesh == "InternalMesh" )
    {
      std::map< string, string > api2geosx = { { "name",         "name" },
                                               { "element_type", "elementTypes" },
                                               { "x_range",      "xCoords" },
                                               { "y_range",      "yCoords" },
                                               { "z_range",      "zCoords" },
                                               { "nx",           "nx" },
                                               { "ny",           "ny" },
                                               { "nz",           "nz" },
                                               { "cell_blocks",  "cellBlockNames" } };
      xml_node result = xmlMesh.append_child( geosxMesh.c_str() );

      result.append_attribute( api2geosx["name"].c_str() ) = get( mesh, "name", "gen-mesh-" + std::to_string( counter++ ) ).c_str();
      for( string k: { "element_type", "x_range", "y_range", "z_range", "nx", "ny", "nz", "cell_blocks" } )
      {
        result.append_attribute( api2geosx[k].c_str() ) = serialize( mesh[k] ).c_str();
      }
    }
    else
    {
      GEOSX_ERROR("not supported");
    }
  }
}

void BuildGeometryNode( YAML::Node const & input, pugi::xml_node & problem )
{
  int counter = 0;

  xml_node xmlGeometry = problem.append_child( "Geometry" );

  std::map< string, string > const conv = { { "box", "Box" } };

  // Let's loop on all the meshes to find the GEOSX mesh.
  for( YAML::Node const & mesh: input["mesh"] )
  {
    if( conv.count( mesh["type"].as< string >() ) == 0 ) continue;

    string const geosxGeo = conv.at( mesh["type"].as< string >() );

    if( geosxGeo == "Box" )
    {
      std::map< string, string > api2geosx = { { "name",         "name" },
                                               { "point_min",      "xMin" },
                                               { "point_max",      "xMax" } };
      xml_node result = xmlGeometry.append_child( geosxGeo.c_str() );

      result.append_attribute( api2geosx["name"].c_str() ) = get( mesh, "name", "gen-box-" + std::to_string( counter++ ) ).c_str();
      for( string k: { "point_min", "point_max" } )
      {
        result.append_attribute( api2geosx[k].c_str() ) = serialize( mesh[k] ).c_str();
      }
    }
    else
    {
      GEOSX_ERROR("not supported");
    }
  }
}

void BuildOutputNode( YAML::Node const & input, pugi::xml_node & problem )
{
  int counter = 0;

  xml_node xmlOutputs = problem.append_child( "Outputs" );
  xml_node xmlEvents = problem.append_child( "Events" );
  xmlEvents.append_attribute( "maxTime" ) = input["simulation"]["max_time"].as< string >().c_str();

  std::map< string, string > const conv = { { "vtk",     "VTK" },
                                            { "silo",    "Silo" },
                                            { "restart", "Restart" } };

  // Let's loop on all the meshes to find the GEOSX mesh.
  for( YAML::Node const & output: input["outputs"] )
  {
    string const geosxOutput = conv.at( output["type"].as< string >() );

    std::map< string, string > api2geosx = { { "name", "name" },
                                             { "exact_timestep", "targetExactTimestep" },
                                             { "period", "timeFrequency" } };

    xml_node xmlOutput = xmlOutputs.append_child( geosxOutput.c_str() );
    string const outputName = get( output, "name", "gen-output-" + std::to_string( counter++ ) );
    xmlOutput.append_attribute( api2geosx["name"].c_str() ) = outputName.c_str();
    xml_node xmlEvent = xmlEvents.append_child( "PeriodicEvent" );
    xmlEvent.append_attribute( "name" ) = ( "gen-event-" + std::to_string( counter++ ) ).c_str();
    xmlEvent.append_attribute( api2geosx["exact_timestep"].c_str() ) = get( output, "exact_timestep", "1" ).c_str();
    for( string k: { "period" } )
    {
      xmlEvent.append_attribute( api2geosx[k].c_str() ) = serialize( output[k] ).c_str();
    }
    xmlEvent.append_attribute( "target" ) = ( "/Outputs/" + outputName ).c_str();

    if( geosxOutput == "Silo" )
    {
      string const parallelThreads = get( output, "parallel_threads", {} );
      if( not parallelThreads.empty() )
      {
        xmlOutput.append_attribute( "parallelThreads" ) = parallelThreads.c_str();
      }
    }
  }
}

void BuildFunctionsNode( YAML::Node const & input, pugi::xml_node & problem )
{
  int counter = 0;

  xml_node xmlGeometry = problem.append_child( "Functions" );

  std::map< string, string > const conv = { { "table", "TableFunction" } };

  // Let's loop on all the meshes to find the GEOSX mesh.
  for( YAML::Node const & function: input["functions"] )
  {
//    if( conv.count( function["type"].as< string >() ) == 0 ) continue;

    string const geosxFct = conv.at( function["type"].as< string >() );

    if( geosxFct == "TableFunction" )
    {
      std::map< string, string > api2geosx = { { "name",        "name" },
                                               { "variables",   "inputVarNames" },
                                               { "coordinates", "coordinates" },
                                               { "values",      "values" } };
      xml_node result = xmlGeometry.append_child( geosxFct.c_str() );

      result.append_attribute( api2geosx["name"].c_str() ) = get( function, "name", "gen-fct-" + std::to_string( counter++ ) ).c_str();
      for( auto kv: api2geosx )
      {
        result.append_attribute( kv.second.c_str() ) = serialize( function[kv.first] ).c_str();
      }
    }
    else
    {
      GEOSX_ERROR("not supported");
    }
  }
}


void BuildFieldsNode( YAML::Node const & input,
                      pugi::xml_node & problem )
{
  int counter = 0;
  std::map< string, string > const fields2geosx = { { "temperature",          "Temperature" },
                                                    { "global_comp_fraction", "globalCompFraction" } };

  std::map< string, string > const loc2mgr = { { "nodes",    "nodeManager" },
                                               { "elements", "ElementRegions" } };

  xml_node xmlGeometry = problem.append_child( "FieldSpecifications" );

//  std::map< string, string > const conv = { { "table", "FieldSpecification" } };

  // Let's loop on all the meshes to find the GEOSX mesh.
  for( YAML::Node const & field: input["sources_and_sinks"] )
  {
    std::map< string, string > const api2geosx = { { "name",     "name" },
                                                   { "field",    "fieldName" },
                                                   { "function", "functionName" },
                                                   { "value",    "scale" } };
    xml_node xmlField = xmlGeometry.append_child( "FieldSpecification" );

    xmlField.append_attribute( api2geosx.at( "name" ).c_str() ) = get( field, "name", "gen-field-" + std::to_string( counter++ ) ).c_str();

    std::vector< string > const location = tokenize( field["location"].as< string >(), "/" );
    xmlField.append_attribute( "setNames" ) = ( "{ " + location.back() + " }" ).c_str();
    xmlField.append_attribute( "objectPath" ) = loc2mgr.at( location.front() ).c_str();

    xmlField.append_attribute( api2geosx.at( "value" ).c_str() ) = get( field,  "value", "1." ).c_str();
    xmlField.append_attribute( api2geosx.at( "field" ).c_str() ) = fields2geosx.at( field["field"].as< string >() ).c_str();



//    field.count( function["type"].as< string >() ) == 0
    for( string k: { "function" } )
    {
      if( field[k] )
      {
        xmlField.append_attribute( api2geosx.at( k ).c_str() ) = field[k].as< string >().c_str();
      }
    }
  }
}


struct ElementRegionData
{
  std::vector< string > regionNames;
  std::vector< string > constitutiveNames;
};

// Returns the regions names
ElementRegionData BuildRegionsNode( YAML::Node const & input, pugi::xml_node & problem )
{
  int counter = 0;
  ElementRegionData erd;

  xml_node xmlRegions = problem.append_child( "ElementRegions" );

  std::map< string, string > const conv = { { "cells", "CellElementRegion" } };

  // Let's loop on all the meshes to find the GEOSX mesh.
  for( YAML::Node const & region: input["regions"] )
  {
//    if( conv.count( function["type"].as< string >() ) == 0 ) continue;

    string const geosxRegion = conv.at( region["type"].as< string >() );

    if( geosxRegion == "CellElementRegion" )
    {
      std::map< string, string > api2geosx = { { "name",        "name" },
                                               { "variables",   "inputVarNames" },
                                               { "coordinates", "coordinates" },
                                               { "values",      "values" } };
      xml_node result = xmlRegions.append_child( geosxRegion.c_str() );

      string const & regionName = get( region, "name", "gen_reg_" + std::to_string( counter++ ) );
      erd.regionNames.push_back( regionName );
      result.append_attribute( api2geosx["name"].c_str() ) = regionName.c_str();
      result.append_attribute( "cellBlocks" ) = serialize( region["cell_blocks"] ).c_str();
      string const & materialName = get( region, "material_list", "gen-mat-" + std::to_string( counter++ ) );
      erd.constitutiveNames.push_back( materialName );
      result.append_attribute( "materialList" ) = ( "{ " + materialName + " }" ).c_str();
    }
    else
    {
      GEOSX_ERROR("not supported");
    }
  }
  return erd;
}

void BuildConstitutiveNode( YAML::Node const & input, ElementRegionData const & erd, pugi::xml_node & problem )
{
  int counter = 0;

  xml_node xmlConstitutive = problem.append_child( "Constitutive" );

  for (string const & constitutiveName: erd.constitutiveNames)
  {
    xml_node result = xmlConstitutive.append_child( "NullModel" );
    result.append_attribute( "name" ) = constitutiveName.c_str();
  }
}

void BuildSolverNode( YAML::Node const & input, ElementRegionData const & erd, pugi::xml_node & problem )
{
  int counter = 0;

  auto solvers = input["simulation"]["solvers"];
  auto laplace = solvers[0];

  xml_node xmlSolvers = problem.append_child( "Solvers" );
  xml_node xmlSolver = xmlSolvers.append_child( "LaplaceFEM" );
  string const & solverName = get( laplace, "name", "gen_solver_" + std::to_string( counter++ ) );
  xmlSolver.append_attribute( "name") = solverName.c_str();
  xmlSolver.append_attribute("timeIntegrationOption") = "SteadyState";
  xmlSolver.append_attribute("fieldName") = "Temperature";
  xmlSolver.append_attribute("discretization") = "FE1";
  xmlSolver.append_attribute( "targetRegions" ) = ( "{ " + join( erd.regionNames, ", " ) + " }" ).c_str(); // Only be default!
  xml_node lsp = xmlSolver.append_child("LinearSolverParameters");
  lsp.append_attribute("directParallel") = "0";

  auto nm = problem.append_child( "NumericalMethods" );
  auto fe = nm.append_child( "FiniteElements" );
  auto fes = fe.append_child( "FiniteElementSpace" );
  fes.append_attribute("name") = "FE1";
  fes.append_attribute("order") = "1";

  auto xmlEvents = problem.child("Events");
  auto solvEvent = xmlEvents.append_child("PeriodicEvent");
  solvEvent.append_attribute("name") = "gen-sol-2134";
  solvEvent.append_attribute("forceDt") = "1.0";
  solvEvent.append_attribute( "target" ) = ( "/Solvers/" + solverName ).c_str();
}


void Convert( string const & stableInputFileName, xmlWrapper::xmlDocument & doc )
{
  xml_node problem = doc.append_child( "Problem" );
//  xml_node problem = doc.append_child("Problem");
  YAML::Node const input = YAML::LoadFile( stableInputFileName );

  string const apiVersion = input["version"].as<string>();
  GEOSX_LOG(apiVersion);

//  xml_node const meshNode = BuildMeshNode( input );
//  xml_node xmlMesh = problem.append_child("Mesh");
  BuildMeshNode( input, problem );
  BuildGeometryNode( input, problem );
  BuildOutputNode( input, problem );
  BuildFunctionsNode( input, problem );
  BuildFieldsNode( input, problem );
  ElementRegionData erd = BuildRegionsNode( input, problem );
  BuildConstitutiveNode( input, erd, problem );
  BuildSolverNode( input, erd, problem );

  doc.save( std::cout );
}

}
}
