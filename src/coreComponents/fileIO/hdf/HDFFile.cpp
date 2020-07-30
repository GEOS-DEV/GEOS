#include "HDFFile.hpp"

#include "managers/Outputs/TimeHistoryOutput.hpp"

#include "mpiCommunications/MpiWrapper.hpp"

#include <hdf5.h>

namespace geosx
{

/**
 * @brief Get the HDF data type for the specified type
 * @tparam T the type to get info for
 * @return The HDF-specific data type info for T.
 */
template< typename T >
inline hid_t GetHDFDataType();
template<>
inline hid_t GetHDFDataType< char >() { return H5T_NATIVE_CHAR; }
template<>
inline hid_t GetHDFDataType< signed char >() { return H5T_NATIVE_CHAR; }
template<>
inline hid_t GetHDFDataType< float >() { return H5T_NATIVE_FLOAT; }
template<>
inline hid_t GetHDFDataType< double >() { return H5T_NATIVE_DOUBLE; }
template<>
inline hid_t GetHDFDataType< int >() { return H5T_NATIVE_INT; }
template<>
inline hid_t GetHDFDataType< long >() { return H5T_NATIVE_LONG; }
template<>
inline hid_t GetHDFDataType< long long >() { return H5T_NATIVE_LLONG; }

/**
 * @brief Get the HDF data type from the type_index of the type.
 * @param type the std::type_index(typeid(T)) of the type T
 * @return The HDF data type.
 */
inline hid_t GetHDFDataType( std::type_index const & type )
{
  if( type == std::type_index( typeid(char)) )
  {
    return GetHDFDataType< char >();
  }
  else if( type == std::type_index( typeid(signed char)) )
  {
    return GetHDFDataType< signed char >();
  }
  else if( type == std::type_index( typeid(real32)) )
  {
    return GetHDFDataType< real32 >();
  }
  else if( type == std::type_index( typeid(real64)) )
  {
    return GetHDFDataType< real64 >();
  }
  else if( type == std::type_index( typeid(integer)) )
  {
    return GetHDFDataType< integer >();
  }
  else if( type == std::type_index( typeid(localIndex)) )
  {
    return GetHDFDataType< localIndex >();
  }
  else if( type == std::type_index( typeid(globalIndex)) )
  {
    return GetHDFDataType< globalIndex >();
  }
  else
  {
    return GetHDFDataType< char >();
  }
}

/**
 * @brief Get an HDF array data type.
 * @param type The std::type_index(typeid(T)) of the type T inside the array
 * @param rank The rank of the array.
 * @param dims The extent of each dimension of the array.
 */
inline hid_t GetHDFArrayDataType( std::type_index const & type, hsize_t const rank, hsize_t const * dims )
{
  return H5Tarray_create( GetHDFDataType( type ), rank, dims );
}

HDFFile::HDFFile( string const & fnm, bool deleteExisting, bool parallelAccess, MPI_Comm comm ):
  m_filename( ),
  m_fileId( 0 ),
  m_faplId( 0 ),
  m_mpioFapl( parallelAccess ),
  m_comm( comm )
{
  int rnk = MpiWrapper::Comm_rank( comm );
  if( m_mpioFapl )
  {
    m_faplId = H5Pcreate( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio( m_faplId, m_comm, MPI_INFO_NULL );
    m_filename = fnm + ".hdf5";
  }
  else
  {
    m_faplId = H5P_DEFAULT;
    m_filename = fnm + "." + std::to_string( rnk ) + ".hdf5";
  }
  // check if file already exists
  htri_t exists = 0;
  H5E_BEGIN_TRY {
    exists = H5Fis_hdf5( m_filename.c_str() );
  } H5E_END_TRY
  // if there is an non-hdf file with the same name,
  // and we're either not using parallel access or we're rank 0
  if( exists == 0 && ( !m_mpioFapl || rnk == 0 ) )
  {
    remove( m_filename.c_str() );
  }
  if( exists > 0 && !deleteExisting )
  {
    m_fileId = H5Fopen( m_filename.c_str(), H5F_ACC_RDWR, m_faplId );
  }
  else if( exists >= 0 && deleteExisting )
  {
    m_fileId = H5Fcreate( m_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, m_faplId );
  }
  else if( exists < 0 )
  {
    m_fileId = H5Fcreate( m_filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, m_faplId );
  }
}

HDFFile::~HDFFile()
{
  if( m_mpioFapl )
  {
    H5Pclose( m_faplId );
  }
  H5Fclose( m_fileId );
}

HDFHistIO::HDFHistIO( string const & filename,
                      localIndex rank,
                      std::vector< localIndex > const & dims,
                      string const & name,
                      std::type_index typeId,
                      localIndex writeHead,
                      localIndex initAlloc,
                      MPI_Comm comm ):
  BufferedHistoryIO(),
  m_filename( filename ),
  m_overallocMultiple( initAlloc ),
  m_globalIdxOffset( 0 ),
  m_globalIdxCount( 0 ),
  m_writeLimit( initAlloc ),
  m_writeHead( writeHead ),
  m_hdfType( GetHDFDataType( typeId )),
  m_typeSize( H5Tget_size( m_hdfType )),
  m_typeCount( 1 ),
  m_rank( LvArray::integerConversion< hsize_t >( rank )),
  m_dims( rank ),
  m_name( name ),
  m_comm( comm ),
  m_subcomm( MPI_COMM_NULL )
{
  for( hsize_t dd = 0; dd < m_rank; ++dd )
  {
    m_dims[dd] = LvArray::integerConversion< hsize_t >( dims[dd] );
    m_typeCount *= m_dims[dd];
  }
  m_dataBuffer.resize( initAlloc * m_typeSize * m_typeCount );
}

void HDFHistIO::init( bool exists_okay )
{
  globalIndex localIdxCount = LvArray::integerConversion< globalIndex >( m_dims[0] );

  int size = MpiWrapper::Comm_size( m_comm );
  int rank = MpiWrapper::Comm_rank( m_comm );

  int color = MPI_UNDEFINED;
  if( localIdxCount > 0 )
  {
    color = 0;
  }

  std::vector< globalIndex > counts( size );
  MpiWrapper::Allgather( &localIdxCount, 1, &counts[0], 1, m_comm );

  hsize_t minIdxCount = std::numeric_limits< hsize_t >::max( );
  int key = 0;
  for( int ii = 0; ii < size; ii++ )
  {
    if( counts[ii] != 0 && ii < rank )
    {
      key++;
    }
    if( ii == rank )
    {
      m_globalIdxOffset = m_globalIdxCount;
    }
    m_globalIdxCount += counts[ii];
    if( counts[ii] > 0 &&  LvArray::integerConversion< hsize_t >( counts[ii] ) < minIdxCount )
    {
      minIdxCount = counts[ii];
    }
  }

  m_subcomm = MpiWrapper::Comm_split( m_comm, color, key );
  // create a dataset in the file if needed, don't erase file
  if( m_subcomm != MPI_COMM_NULL )
  {
    std::vector< hsize_t > historyFileDims( m_rank+1 );
    historyFileDims[0] = LvArray::integerConversion< hsize_t >( m_writeLimit );

    std::vector< hsize_t > dimChunks( m_rank+1 );
    dimChunks[0] = 1;

    for( hsize_t dd = 1; dd < m_rank+1; ++dd )
    {
      // a process with chunk size 0 is considered incorrect by hdf5, hence the subcomm
      dimChunks[dd] = m_dims[dd-1];
      historyFileDims[dd] = m_dims[dd-1];
    }
    dimChunks[1] = minIdxCount;
    historyFileDims[1] = LvArray::integerConversion< hsize_t >( m_globalIdxCount );

    HDFFile target( m_filename, false, true, m_subcomm );

    bool inTarget = target.CheckInTarget( m_name );
    if( !inTarget )
    {
      hid_t dcplId = 0;
      std::vector< hsize_t > maxFileDims( historyFileDims );
      // chunking is required to create an extensible dataset
      dcplId = H5Pcreate( H5P_DATASET_CREATE );
      H5Pset_chunk( dcplId, m_rank+1, &dimChunks[0] );
      maxFileDims[0] = H5S_UNLIMITED;
      hid_t space = H5Screate_simple( m_rank+1, &historyFileDims[0], &maxFileDims[0] );
      hid_t dataset = H5Dcreate( target, m_name.c_str(), m_hdfType, space, H5P_DEFAULT, dcplId, H5P_DEFAULT );
      H5Dclose( dataset );
      H5Sclose( space );
    }
    else if( exists_okay )
    {
      // todo:
      // hid_t dataset = H5Dopen(target, m_name.c_str( ), H5P_DEFAULT);
      // check that the extent of the filespace is compatible with the data
    }
    GEOSX_ERROR_IF( inTarget && !exists_okay, "Dataset (" + m_name + ") already exists in output file: " + m_filename );
  }
}

void HDFHistIO::write( )
{
  // don't need to write if nothing is buffered, this should only happen if the output event occurs before the collection event
  if( m_subcomm != MPI_COMM_NULL )
  {
    localIndex maxBuffered = 0;
    MpiWrapper::allReduce( &m_bufferedCount, &maxBuffered, 1, MPI_MAX, m_subcomm );
    if( maxBuffered > 0 )
    {
      resizeFileIfNeeded( maxBuffered );
      HDFFile target( m_filename, false, true, m_subcomm );

      hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
      hid_t filespace = H5Dget_space( dataset );

      std::vector< hsize_t > fileOffset( m_rank+1 );
      fileOffset[0] = LvArray::integerConversion< hsize_t >( m_writeHead );
      fileOffset[1] = LvArray::integerConversion< hsize_t >( m_globalIdxOffset );

      std::vector< hsize_t > bufferedCounts( m_rank+1 );
      bufferedCounts[0] = LvArray::integerConversion< hsize_t >( maxBuffered );
      for( hsize_t dd = 1; dd < m_rank+1; ++dd )
      {
        bufferedCounts[dd] = m_dims[dd-1];
      }
      hid_t memspace = H5Screate_simple( m_rank+1, &bufferedCounts[0], nullptr );

      hid_t fileHyperslab = filespace;
      H5Sselect_hyperslab( fileHyperslab, H5S_SELECT_SET, &fileOffset[0], nullptr, &bufferedCounts[0], nullptr );

      buffer_unit_type * dataBuffer = nullptr;
      // if local rank is writting nothing, m_dataBuffer is never alloc'd so don't try to access it
      if( m_typeCount != 0 )
      {
        dataBuffer = &m_dataBuffer[0];
      }
      H5Dwrite( dataset, m_hdfType, memspace, fileHyperslab, H5P_DEFAULT, dataBuffer );

      H5Sclose( memspace );
      H5Sclose( filespace );
      H5Dclose( dataset );

      m_writeHead += maxBuffered;
      emptyBuffer( );
    }
  }
}

void HDFHistIO::compressInFile( )
{
  if( m_subcomm != MPI_COMM_NULL )
  {
    HDFFile target( m_filename, false, true, m_subcomm );
    std::vector< hsize_t > maxFileDims( m_rank+1 );
    maxFileDims[0] = LvArray::integerConversion< hsize_t >( m_writeHead );
    maxFileDims[1] = LvArray::integerConversion< hsize_t >( m_globalIdxCount );
    for( hsize_t dd = 2; dd < m_rank+1; ++dd )
    {
      maxFileDims[dd] = m_dims[dd-1];
    }
    hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
    H5Dset_extent( dataset, &maxFileDims[0] );
    H5Dclose( dataset );
    m_writeLimit = m_writeHead;
  }
}

inline void HDFHistIO::resizeFileIfNeeded( localIndex buffered_count )
{
  if( m_subcomm != MPI_COMM_NULL )
  {
    HDFFile target( m_filename, false, true, m_subcomm );
    if( m_writeHead + buffered_count > m_writeLimit )
    {
      while( m_writeHead + buffered_count > m_writeLimit )
      {
        m_writeLimit *= m_overallocMultiple;
      }
      std::vector< hsize_t > maxFileDims( m_rank+1 );
      maxFileDims[0] = LvArray::integerConversion< hsize_t >( m_writeLimit );
      maxFileDims[1] = LvArray::integerConversion< hsize_t >( m_globalIdxCount );
      for( hsize_t dd = 2; dd < m_rank+1; ++dd )
      {
        maxFileDims[dd] = m_dims[dd-1];
      }
      hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
      H5Dset_extent( dataset, &maxFileDims[0] );
      H5Dclose( dataset );
    }
  }
}

void HDFHistIO::resizeBuffer( )
{
  size_t osize = m_dataBuffer.size();
  // if needed, resize the buffer
  if( (m_bufferedCount + 1) * (m_typeCount * m_typeSize) > osize )
  {
    m_dataBuffer.resize( osize + ( m_typeCount * m_typeSize ) * m_overallocMultiple );
  }
  // advance the buffer head
  m_bufferHead = (&m_dataBuffer[0]) + m_bufferedCount * ( m_typeCount * m_typeSize );
}



HDFSerialHistIO::HDFSerialHistIO( string const & filename,
                                  localIndex rank,
                                  std::vector< localIndex > const & dims,
                                  string const & name,
                                  std::type_index typeId,
                                  localIndex writeHead,
                                  localIndex initAlloc,
                                  MPI_Comm comm ):
  BufferedHistoryIO(),
  m_filename( filename ),
  m_overallocMultiple( initAlloc ),
  m_writeLimit( initAlloc ),
  m_writeHead( writeHead ),
  m_hdfType( GetHDFDataType( typeId )),
  m_typeSize( H5Tget_size( m_hdfType )),
  m_typeCount( 1 ),
  m_rank( LvArray::integerConversion< hsize_t >( rank )),
  m_dims( rank ),
  m_name( name ),
  m_comm( comm )
{
  for( hsize_t dd = 0; dd < m_rank; ++dd )
  {
    m_dims[dd] = LvArray::integerConversion< hsize_t >( dims[dd] );
    m_typeCount *= m_dims[dd];
  }
  m_dataBuffer.resize( initAlloc * m_typeSize * m_typeCount );
}

void HDFSerialHistIO::init( bool exists_okay )
{
  // create a dataset in the file if needed, don't erase file
  if( m_typeCount > 0 )
  {
    std::vector< hsize_t > historyFileDims( m_rank+1 );
    historyFileDims[0] = LvArray::integerConversion< hsize_t >( m_writeLimit );

    std::vector< hsize_t > dimChunks( m_rank+1 );
    dimChunks[0] = 1;

    for( hsize_t dd = 1; dd < m_rank+1; ++dd )
    {
      // a process with chunk size 0 is considered incorrect by hdf5, hence the subcomm
      dimChunks[dd] = m_dims[dd-1];
      historyFileDims[dd] = m_dims[dd-1];
    }

    // why is this killing things, only in this context, only on lassen?
    HDFFile target( m_filename, false, false, m_comm );
    //
    bool inTarget = target.CheckInTarget( m_name );
    if( !inTarget )
    {
      hid_t dcplId = 0;
      std::vector< hsize_t > maxFileDims( historyFileDims );
      // chunking is required to create an extensible dataset
      dcplId = H5Pcreate( H5P_DATASET_CREATE );
      H5Pset_chunk( dcplId, m_rank+1, &dimChunks[0] );
      maxFileDims[0] = H5S_UNLIMITED;
      hid_t space = H5Screate_simple( m_rank+1, &historyFileDims[0], &maxFileDims[0] );
      hid_t dataset = H5Dcreate( target, m_name.c_str(), m_hdfType, space, H5P_DEFAULT, dcplId, H5P_DEFAULT );
      H5Dclose( dataset );
      H5Sclose( space );
    }
    else if( exists_okay )
    {
      // todo:
      // hid_t dataset = H5Dopen(target, m_name.c_str( ), H5P_DEFAULT);
      // check that the extent of the filespace is compatible with the data
    }
    GEOSX_ERROR_IF( inTarget && !exists_okay, "Dataset (" + m_name + ") already exists in output file: " + m_filename );
  }
}

void HDFSerialHistIO::write( )
{
  // don't need to write if nothing is buffered, this should only happen if the output event occurs before the collection event
  if( m_typeCount > 0 && m_bufferedCount > 0 )
  {
    resizeFileIfNeeded( m_bufferedCount );
    HDFFile target( m_filename, false, false, m_comm );

    hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
    hid_t filespace = H5Dget_space( dataset );

    std::vector< hsize_t > fileOffset( m_rank+1 );
    fileOffset[0] = LvArray::integerConversion< hsize_t >( m_writeHead );
    fileOffset[1] = 0;

    std::vector< hsize_t > bufferedCounts( m_rank+1 );
    bufferedCounts[0] = LvArray::integerConversion< hsize_t >( m_bufferedCount );
    for( hsize_t dd = 1; dd < m_rank+1; ++dd )
    {
      bufferedCounts[dd] = m_dims[dd-1];
    }
    hid_t memspace = H5Screate_simple( m_rank+1, &bufferedCounts[0], nullptr );

    hid_t fileHyperslab = filespace;
    H5Sselect_hyperslab( fileHyperslab, H5S_SELECT_SET, &fileOffset[0], nullptr, &bufferedCounts[0], nullptr );

    buffer_unit_type * dataBuffer = nullptr;
    // if local rank is writting nothing, m_dataBuffer is never alloc'd so don't try to access it
    if( m_typeCount != 0 )
    {
      dataBuffer = &m_dataBuffer[0];
    }
    H5Dwrite( dataset, m_hdfType, memspace, fileHyperslab, H5P_DEFAULT, dataBuffer );

    H5Sclose( memspace );
    H5Sclose( filespace );
    H5Dclose( dataset );

    m_writeHead += m_bufferedCount;
    emptyBuffer( );
  }
}

void HDFSerialHistIO::compressInFile( )
{
  if( m_typeCount > 0 )
  {
    HDFFile target( m_filename, false, false, m_comm );
    std::vector< hsize_t > maxFileDims( m_rank+1 );
    maxFileDims[0] = LvArray::integerConversion< hsize_t >( m_writeHead );
    maxFileDims[1] = m_dims[0];
    for( hsize_t dd = 2; dd < m_rank+1; ++dd )
    {
      maxFileDims[dd] = m_dims[dd-1];
    }
    hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
    H5Dset_extent( dataset, &maxFileDims[0] );
    H5Dclose( dataset );
    m_writeLimit = m_writeHead;
  }
}

inline void HDFSerialHistIO::resizeFileIfNeeded( localIndex buffered_count )
{
  if( m_typeCount > 0 )
  {
    HDFFile target( m_filename, false, false, m_comm );
    if( m_writeHead + buffered_count > m_writeLimit )
    {
      while( m_writeHead + buffered_count > m_writeLimit )
      {
        m_writeLimit *= m_overallocMultiple;
      }
      std::vector< hsize_t > maxFileDims( m_rank+1 );
      maxFileDims[0] = LvArray::integerConversion< hsize_t >( m_writeLimit );
      maxFileDims[1] = m_dims[0];
      for( hsize_t dd = 2; dd < m_rank+1; ++dd )
      {
        maxFileDims[dd] = m_dims[dd-1];
      }
      hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
      H5Dset_extent( dataset, &maxFileDims[0] );
      H5Dclose( dataset );
    }
  }
}

void HDFSerialHistIO::resizeBuffer( )
{
  size_t osize = m_dataBuffer.size();
  // if needed, resize the buffer
  if( (m_bufferedCount + 1) * (m_typeCount * m_typeSize) > osize )
  {
    m_dataBuffer.resize( osize + ( m_typeCount * m_typeSize ) * m_overallocMultiple );
  }
  // advance the buffer head
  m_bufferHead = (&m_dataBuffer[0]) + m_bufferedCount * ( m_typeCount * m_typeSize );
}

}
