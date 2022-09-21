#include "FaceBlock.hpp"

namespace geosx
{

ArrayOfArrays< localIndex > myConvert( std::vector< std::vector< localIndex > > const & vv )
{
  ArrayOfArrays< localIndex > res;

  for( std::size_t i = 0; i < vv.size(); ++i )
  {
    auto const & vvv = vv[i];
    res.appendArray( vvv.begin(), vvv.end() );
  }

  return res;
}

array2d< localIndex > myConvert2d( std::vector< std::vector< localIndex > > const & vv )
{
  for( auto const & vvv: vv )
  {
    GEOSX_ERROR_IF_NE_MSG( vvv.size(), 2, "my conversion function sucks." );
  }
  array2d< localIndex > res( vv.size(), 2 );

  for( std::size_t i = 0; i < vv.size(); ++i )
  {
    for( std::size_t j = 0; j < vv[i].size(); ++j )
    {
      res[i][j] = vv[i][j];
    }
  }

  return res;
}

array1d< localIndex > myConvert( std::vector< localIndex > const & v )
{
  array1d< localIndex > result;
  for( std::size_t i = 0; i < v.size(); ++i )
  {
    result.emplace_back( v[i] );
  }
  return result;
}


using namespace dataRepository;

//FaceBlock::FaceBlock( string const & name,
//                      Group * const parent )
//  :
//  FaceBlockABC( name, parent )
//{ }
//
//localIndex FaceBlock::num2dElements() const
//{
//  return 10;
//}
//
//ArrayOfArrays< localIndex > FaceBlock::get2dElemToNodes() const
//{
//  std::vector< std::vector< localIndex > > result;
//  for( std::size_t i = 0; i < 10; ++i )
//  {
//    std::vector< localIndex > tmp( 8, -1 );
//    tmp[0] = 12 * i + 5;
//    tmp[1] = 12 * ( i + 1 ) + 5;
//    tmp[3] = 12 * ( i + 1 ) + 5 + 132;
//    tmp[2] = 12 * i + 5 + 132;
//    tmp[4] = 12 * i + 5 + 1;
//    tmp[5] = 12 * ( i + 1 ) + 5 + 1;
//    tmp[7] = 12 * ( i + 1 ) + 5 + 132 + 1;
//    tmp[6] = 12 * i + 5 + 132 + 1;
//    result.push_back( tmp );
//  }
//  return myConvert( result );
//}
//
//ArrayOfArrays< localIndex > FaceBlock::get2dElemToEdges() const
//{
//  std::map< localIndex, std::vector< localIndex > > const f2ed{ // face to edges, WRONG!
//    { 15,  { 15,  50,  372, 16 } },
//    { 47,  { 49,  84,  394, 50 } },
//    { 79,  { 83,  118, 416, 84 } },
//    { 111, { 117, 152, 438, 118 } },
//    { 143, { 151, 186, 460, 152 } },
//    { 175, { 185, 220, 482, 186 } },
//    { 207, { 219, 254, 504, 220 } },
//    { 239, { 253, 288, 526, 254 } },
//    { 271, { 287, 322, 548, 288 } },
//    { 303, { 321, 350, 570, 322 } }
//  };
//
////  std::vector< localIndex > const leftCommonFaces{ 205, 209, 213, 217, 221, 225, 229, 233, 237, 241 };
//  std::vector< localIndex > const leftCommonFaces{ 15, 47, 79, 111, 143, 175, 207, 239, 271, 303 };
//
//  localIndex const numNodesPerElement = 4;
//  ArrayOfArrays< localIndex > result;
//  result.resize( this->num2dElements(), numNodesPerElement );
//
//  for( localIndex i = 0; i < this->num2dElements(); ++i )
//  {
//    result.resizeArray( i, numNodesPerElement );
//    for( std::size_t j = 0; j < numNodesPerElement; ++j )
//    {
//      result[i][j] = f2ed.at( leftCommonFaces[i] )[j];
//    }
//  }
//
//  return result;
//}
//
//ToCellRelation< array2d< localIndex > > FaceBlock::get2dElemToElems() const
//{
//  ToCellRelation< array2d< localIndex > > result;
//
//  std::vector< std::vector< localIndex > > toElements;
//  for( localIndex i = 0; i < 10; ++i )
//  {
//    toElements.push_back( { 4 + i * 10, 5 + i * 10 } );
//  }
//
//  std::vector< std::vector< localIndex > > toCellBlock;
//  for( localIndex i = 0; i < 10; ++i )
//  {
//    toCellBlock.push_back( { 0, 0 } );
//  }
//
//  return { myConvert2d( toCellBlock ), myConvert2d( toElements ) };
//}
//
//array2d< localIndex > FaceBlock::get2dElemToFaces() const
//{
//  std::vector< std::vector< localIndex > > result;
//  for( localIndex i = 0; i < 10; ++i )
//  {
//    result.push_back( std::vector< localIndex >{ 15 + 32 * i, 15 + 3 + 32 * i } );
//  }
//
//  return myConvert2d( result );
//}
//
//array1d< localIndex > FaceBlock::get2dFaceToEdge() const
//{
//  std::vector< localIndex > result;
//  for( int i = 0; i < num2dElements(); ++i )
//  {
//    result.push_back( 15 + i * 34 );
//    result.push_back( 16 + i * 34 );
//    result.push_back( 372 + i * 22 );
//  }
//  result.push_back( 350 );
//
//  return myConvert( result );
//}
//
//localIndex FaceBlock::num2dFaces() const
//{
//  return 31;
//}
//
//ArrayOfArrays< localIndex > FaceBlock::get2dFaceTo2dElems() const
//{
//  ArrayOfArrays< localIndex > result;
//  result.resize( num2dFaces() );
//  for( int i = 0; i < num2dFaces(); ++i )
//  {
////    if( i > 0 and ( i % 3 ) == 0 and i < 30 )
//    if( ( ( i - 1 ) % 3 ) == 0 and i < 30 and i > 1 )
//    {
//      result.resizeArray( i, 2 );
//      result[i][0] = localIndex( i / 3 ) - 1;
//      result[i][1] = localIndex( i / 3 );
//    }
//    else
//    {
//      result.resizeArray( i, 1 );
//      result[i][0] = localIndex( i / 3 );
//    }
//  }
//  result[30][0] = 9;
//
//  return result;
//}

localIndex FaceBlock::num2dElements() const
{
  return m_num2dElements;
}

localIndex FaceBlock::num2dFaces() const
{
  return m_num2dFaces;
}

ArrayOfArrays< localIndex > FaceBlock::get2dElemToNodes() const
{
  return myConvert( m_2dElemToNodes );
}

ArrayOfArrays< localIndex > FaceBlock::get2dElemToEdges() const
{
  return myConvert( m_2dElemToEdges );
}

array2d< localIndex > FaceBlock::get2dElemToFaces() const
{
  return myConvert2d( m_2dElemToFaces );
}

ToCellRelation< array2d< localIndex > > FaceBlock::get2dElemToElems() const
{
  array2d< localIndex > converted = myConvert2d( m_2dElemToElems );
  array2d< localIndex > blocks( converted );
  blocks.setValues< serialPolicy >( 0 );

  return { blocks, converted };
}

array1d< localIndex > FaceBlock::get2dFaceToEdge() const
{
  return myConvert( m_2dFaceToEdge );
}

ArrayOfArrays< localIndex > FaceBlock::get2dFaceTo2dElems() const
{
  return myConvert( m_2dFaceTo2dElems );
}

void FaceBlock::setNum2DElements( localIndex num2DElements )
{
  m_num2dElements = num2DElements;
}

void FaceBlock::setNum2DFaces( localIndex num2DFaces )
{
  m_num2dFaces = num2DFaces;
}

void FaceBlock::set2dElemToNodes(std::vector< std::vector< localIndex > > _2dElemToNodes)
{
  m_2dElemToNodes = _2dElemToNodes;
}

void FaceBlock::set2dElemToEdges(std::vector< std::vector< localIndex > > _2dElemToEdges)
{
  m_2dElemToEdges = _2dElemToEdges;
}

void FaceBlock::set2dElemToFaces(std::vector< std::vector< localIndex > > _2dElemToFaces)
{
  m_2dElemToFaces = _2dElemToFaces;
}

void FaceBlock::set2dFaceTo2dElems( std::vector< std::vector< localIndex > > _2dFaceTo2dElems)
{
  m_2dFaceTo2dElems = _2dFaceTo2dElems;
}

void FaceBlock::set2dFaceToEdge( std::vector< localIndex > _2dFaceToEdge)
{
  m_2dFaceToEdge = _2dFaceToEdge;
}

void FaceBlock::set2dElemToElems(std::vector< std::vector< localIndex > > _2dElemToElems)
{
  m_2dElemToElems = _2dElemToElems;
}

}