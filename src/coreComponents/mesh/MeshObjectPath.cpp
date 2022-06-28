/*
 * MeshPath.cpp
 *
 *  Created on: Jun 21, 2022
 *      Author: settgast
 */

#include "MeshObjectPath.hpp"
#include "MeshBody.hpp"

namespace geosx
{

MeshObjectPath::MeshObjectPath( string const path,
                                dataRepository::Group const & meshBodies ):
    m_objectType( getObjectType(path) ),
    m_pathPermutations()
{
  processPath( path, meshBodies );
}

MeshObjectPath::~MeshObjectPath()
{
  // TODO Auto-generated destructor stub
}

string MeshObjectPath::getObjectType( string const & path )
{
  string objectType;
  if( path.find( MeshLevel::groupStructKeys::nodeManagerString() ) != std::string::npos )
  {
    objectType = MeshLevel::groupStructKeys::nodeManagerString();
  }
  else if( path.find( MeshLevel::groupStructKeys::edgeManagerString() ) != std::string::npos )
  {
    objectType = MeshLevel::groupStructKeys::edgeManagerString();
  }
  else if( path.find( MeshLevel::groupStructKeys::faceManagerString() ) != std::string::npos )
  {
    objectType = MeshLevel::groupStructKeys::faceManagerString();
  }
  else if( path.find( MeshLevel::groupStructKeys::elemManagerString() ) != std::string::npos )
  {
    objectType = MeshLevel::groupStructKeys::elemManagerString();
  }
  return objectType;
}

void MeshObjectPath::fillPathTokens( string const & path,
                                     std::vector<string> & pathTokens,
                                     dataRepository::Group const & meshBodies )
{
  // find where the object specification is in the path
  auto findObjectIndex = [&]() -> int
  {
    for( size_t a=0; a<pathTokens.size(); ++a )
    {
      if( pathTokens[a] == m_objectType )
      {
        return a;
      }
    }
    return -1;
  };

  int objectIndex = findObjectIndex();

  GEOSX_ERROR_IF( objectIndex==-1,
                  GEOSX_FMT( "Path ({}) does not contain a valid object type. "
                             "Must contain one of ({},{},{},{})",
                             path,
                             MeshLevel::groupStructKeys::nodeManagerString(),
                             MeshLevel::groupStructKeys::edgeManagerString(),
                             MeshLevel::groupStructKeys::faceManagerString(),
                             MeshLevel::groupStructKeys::elemManagerString() ) );

  // No MeshBody or MeshLevels were specified. add all of them
  if( objectIndex==0 )
  {
    pathTokens.insert( pathTokens.begin(), "{:}" );
    pathTokens.insert( pathTokens.begin(), "{:}" );
  }
  // MeshBody OR MeshLevel specified. Check which one, and add all of the other.
  else if( objectIndex==1 )
  {
    string const unidentifiedName = pathTokens[0];
    // if the MeshBody is specified, add all MeshLevels
    if( meshBodies.hasGroup( unidentifiedName ) )
    {
      pathTokens.insert( pathTokens.begin(), unidentifiedName );
      pathTokens.insert( pathTokens.begin()+1, "{:}" );
    }
    // It wasn't the MeshBody that was specified, it was the MeshLevel?? Check and add.
    else
    {
      pathTokens.insert( pathTokens.begin(), "{:}" );

      string existingMeshBodyAndLevel;
      bool levelNameFound = false;
      meshBodies.forSubGroups<MeshBody>( [&]( MeshBody const & meshBody )
      {
        existingMeshBodyAndLevel += meshBody.getName() + ": ";
        meshBody.forMeshLevels( [&](MeshLevel const & meshLevel)
        {
          existingMeshBodyAndLevel += meshLevel.getName() + ", ";
          levelNameFound = ( unidentifiedName==meshLevel.getName() ) ? true : levelNameFound;
        });
        existingMeshBodyAndLevel += "/n";
      } );

      GEOSX_ERROR_IF( !levelNameFound,
                      GEOSX_FMT( "Path ({}) specifies an invalid MeshBody or MeshLevel. ",
                                 "existing MeshBodies: MeshLevels /n",
                                 path,
                                 existingMeshBodyAndLevel) );
      pathTokens.insert( pathTokens.begin()+1, unidentifiedName );
    }
  }

  // checks to make sure we have properly inserted thus far
  objectIndex = findObjectIndex();
  size_t targetTokenLength = pathTokens.size();

  GEOSX_ERROR_IF_NE_MSG( objectIndex, 2,
                         "Filling of MeshBody and/or MeshLevel in path has failed. Object Index should be 2" );

  GEOSX_ERROR_IF( targetTokenLength < 2,
                  "Filling of MeshBody and/or MeshLevel in path has failed. targetTokenLength should be greater than 2" );

  // now we need to fill in any missing region/subregion specifications.

  if( m_objectType == MeshLevel::groupStructKeys::elemManagerString() )
  {
    // there are no regions specified
    if( targetTokenLength == 3 )
    {
      pathTokens.push_back( "{:}" );
      pathTokens.push_back( "{:}" );
    }
    // there are no subregions specified
    else if( targetTokenLength == 4)
    {
      pathTokens.push_back( "{:}" );
    }
  }


}


void MeshObjectPath::ExpandPathTokens( std::vector<string> & pathTokens,
                                       dataRepository::Group const & meshBodies )
{
  for( size_t a=0; a<pathTokens.size(); ++a )
  {
    pathTokens[a].erase( std::remove( pathTokens[a].begin(), pathTokens[a].end(), '{'), pathTokens[a].end());
    pathTokens[a].erase( std::remove( pathTokens[a].begin(), pathTokens[a].end(), '}'), pathTokens[a].end());
    pathTokens[a].erase( std::remove( pathTokens[a].begin(), pathTokens[a].end(), ' '), pathTokens[a].end());
  }

  // **** MeshBody *****
  if( pathTokens[0] == ":" )
  {
    meshBodies.forSubGroups<MeshBody>( [&]( MeshBody const & meshBody )
    {
      m_pathPermutations[ meshBody.getName() ];
    } );

  }

  printPermutations();


  // **** MeshLevel *****

  for( auto & meshBodyPair : m_pathPermutations )
  {
    string const & meshBodyName = meshBodyPair.first;
    auto & meshBodyNode = meshBodyPair.second;
    MeshBody const & meshBody = meshBodies.getGroup<MeshBody>( meshBodyName );

    if( pathTokens[1] == ":" )
    {
      meshBody.forMeshLevels([&](MeshLevel const & meshLevel)
      {
        auto & meshLevelNode = meshBodyNode[ meshLevel.getName() ];
        expandElementRegions( pathTokens, meshLevel, meshLevelNode );
      });
    }
    else
    {
      for( auto & meshLevelName : stringutilities::tokenize< std::vector<string> >( pathTokens[1], "," ) )
      {
        GEOSX_ERROR_IF( !(meshBody.getMeshLevels().hasGroup<MeshLevel>(meshLevelName)),
                        "Inconsistent MeshLevel specifications in path. "
                        "Specified MeshLevel(s) must exist in all specified MeshBody(s)."
                        "The list of specified MeshLevel are: \n"<<
                        pathTokens[1]<<std::endl<<
                        "while the MeshBody \""<<meshBody.getName()<<
                        "\" does not contain MeshLevel \""<<meshLevelName<<"\".");

        auto & meshLevelNode = meshBodyNode[ meshLevelName ];
        expandElementRegions( pathTokens,
                              meshBody.getMeshLevel(meshLevelName),
                              meshLevelNode );
      }
    }
  }
  printPermutations();
}

void MeshObjectPath::expandElementRegions( std::vector<string> const & pathTokens,
                                           MeshLevel const & meshLevel,
                                           std::map< string, std::vector< string > > & meshLevelNode )
{
  // ***** ElementRegions *****
  if( m_objectType == MeshLevel::groupStructKeys::elemManagerString() )
  {
    // expand element regions
    ElementRegionManager const & elemRegions = meshLevel.getElemManager();
    if( pathTokens[3] == ":" )
    {
      elemRegions.forElementRegions( [&]( ElementRegionBase const & elemRegion )
      {
        auto & elemRegionNode = meshLevelNode[ elemRegion.getName() ];
        expandElementSubRegions( pathTokens, elemRegion, elemRegionNode );
      });
    }
    else
    {
      for( auto & regionName : stringutilities::tokenize< std::vector<string> >( pathTokens[3], "," ) )
      {
        GEOSX_ERROR_IF( !(elemRegions.hasRegion(regionName)),
                        "Inconsistent ElementRegion specifications in path. "
                        "Specified ElementRegion(s) must exist in all specified MeshLevel(s)."
                        "The list of specified ElementRegion(s) are: \n"<<
                        pathTokens[3]<<std::endl<<
                        "while the MeshLevel \""<<meshLevel.getName()<<
                        "\" does not contain ElementRegion \""<<regionName<<"\".");
        auto & elemRegionNode =meshLevelNode[ regionName ];
        expandElementSubRegions( pathTokens,
                                 elemRegions.getRegion(regionName),
                                 elemRegionNode );
      }
    }
  }
}

void MeshObjectPath::expandElementSubRegions( std::vector<string> const & pathTokens,
                                              ElementRegionBase const & elemRegion,
                                              std::vector< string > & elemRegionNode )
{
  if( pathTokens[4] == ":" )
  {
    elemRegion.forElementSubRegions( [&]( ElementSubRegionBase const & subRegion )
    {
      elemRegionNode.push_back( subRegion.getName() );
    });
  }
  else
  {
    for( auto & subRegionName : stringutilities::tokenize< std::vector<string> >( pathTokens[4], "," ) )
    {
      GEOSX_ERROR_IF( !(elemRegion.hasSubRegion(subRegionName)),
                      "Inconsistent ElementSubRegion specifications in path. "
                      "Specified ElementSubRegion(s) must exist in all specified ElementRegion(s)."
                      "The list of specified ElementSubRegion(s) are: \n"<<
                      pathTokens[4]<<std::endl<<
                      "while the ElementRegion \""<<elemRegion.getName()<<
                      "\" does not contain ElementSubRegion \""<<subRegionName<<"\".");
      elemRegionNode.push_back( subRegionName  );
    }
  }
}


void MeshObjectPath::processPath( string const objectPath,
                                  dataRepository::Group const & meshBodies )
{
  std::vector<string> pathTokens = stringutilities::tokenize< std::vector<string> >( objectPath, "/" );
  fillPathTokens( objectPath, pathTokens, meshBodies );
  ExpandPathTokens( pathTokens, meshBodies );


//  return targetTokens;
}

void MeshObjectPath::printPermutations() const
{
  for( auto const & meshBodyPair : m_pathPermutations )
  {
    std::cout<<meshBodyPair.first<<": "<<std::endl;
    for( auto const & meshLevelPair : meshBodyPair.second )
    {
      std::cout<<"  "<<meshLevelPair.first<<": "<<std::endl;
      for( auto const & elemRegionPair : meshLevelPair.second )
      {
        std::cout<<"    "<<elemRegionPair.first<<": "<<std::endl;
      }
    }
  }
  std::cout<<std::endl<<std::endl;
}


} /* namespace geosx */
