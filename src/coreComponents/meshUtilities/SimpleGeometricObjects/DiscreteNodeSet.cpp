/*
 * DiscreteNodeSet.cpp
 *
 *  Created on: Feb 24, 2020
 *      Author: ron
 */

#include "DiscreteNodeSet.hpp"

namespace geosx
{
using namespace dataRepository;

DiscreteNodeSet::DiscreteNodeSet( const std::string& name, Group * const parent ):
		  SimpleGeometricObjectBase( name, parent ),
          m_nodeSetFile()
{
	 registerWrapper( viewKeyStruct::nodeSetFileString, &m_nodeSetFile, false )->
	    setInputFlag(InputFlags::REQUIRED)->
	    setDescription("node set file name of the DiscreteNodeSet");
}

DiscreteNodeSet::~DiscreteNodeSet()
{
}

template< typename T >
void DiscreteNodeSet::parse_file( array1d<T> & target, string const & filename, char delimiter )
{
  std::ifstream inputStream(filename.c_str());
  std::string lineString;
  T value;

  GEOSX_ERROR_IF( !inputStream, "Could not read input file: " << filename );

  // Read the file
  // TODO: Update this to handle large parallel jobs
  while (!inputStream.eof())
  {
    std::getline(inputStream, lineString);
    std::istringstream ss( lineString );

    while(ss.peek() == delimiter || ss.peek() == ' ')
    {
      ss.ignore();
    }
    while( ss>>value )
    {
      target.push_back( value );
      while(ss.peek() == delimiter || ss.peek() == ' ')
      {
        ss.ignore();
      }
    }
  }

  inputStream.close();
}

bool DiscreteNodeSet::IsCoordInObject( const R1Tensor& coord ) const
{
	R1Tensor coord0(coord);
	bool rval = false;
	return rval;
}

void DiscreteNodeSet::PostProcessInput()
{
	parse_file( m_nodeIndexes, m_nodeSetFile, ',' );
	m_nodeNum = m_nodeIndexes.size();
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, DiscreteNodeSet, std::string const &, Group * const )

} /* namespace geosx */
