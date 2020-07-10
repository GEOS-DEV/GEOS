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
inline hid_t GetHDFDataType< real32 >() { return H5T_NATIVE_FLOAT; }
template<>
inline hid_t GetHDFDataType< real64 >() { return H5T_NATIVE_DOUBLE; }
template<>
inline hid_t GetHDFDataType< integer >() { return H5T_NATIVE_INT; }
template<>
inline hid_t GetHDFDataType< localIndex >() { return H5T_NATIVE_LONG; }
template<>
inline hid_t GetHDFDataType< globalIndex >() { return H5T_NATIVE_LLONG; }

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

HDFFile::HDFFile( string const & fnm, bool delete_existing, bool parallel_access, MPI_Comm comm ):
  m_filename( ),
  m_file_id( 0 ),
  m_fapl_id( 0 ),
  m_mpio_fapl( parallel_access ),
  m_comm( comm )
{
  int rnk = MpiWrapper::Comm_rank( comm );
  if ( m_mpio_fapl )
  {
    m_fapl_id = H5Pcreate( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio( m_fapl_id, m_comm, MPI_INFO_NULL );
    m_filename = fnm;
  }
  else
  {
    m_filename = fnm + "." + std::to_string(rnk) + ".hdf5";
  }
  // check if file already exists
  htri_t exists = 0;
  H5E_BEGIN_TRY {
    exists = H5Fis_hdf5( m_filename.c_str() );
  } H5E_END_TRY
  if( exists > 0 && !delete_existing )
  {
    m_file_id = H5Fopen( m_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
  }
  else if( exists >= 0 && delete_existing )
  {
    m_file_id = H5Fcreate( m_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  }
  else if( exists < 0 )
  {
    m_file_id = H5Fcreate( m_filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT );
  }
  GEOSX_ERROR_IF( exists == 0, string( "Existing file " ) + m_filename + string( " is not an HDF5 file, cannot use for HDF5 output." ) );
}

HDFFile::~HDFFile()
{
  if ( m_mpio_fapl )
  {
    H5Pclose( m_fapl_id );
  }
  H5Fclose( m_file_id );
}

HDFHistIO::HDFHistIO( string const & filename,
                      localIndex rank,
                      const localIndex * dims,
                      string const & name,
                      std::type_index type_id,
                      localIndex write_head,
                      localIndex init_alloc,
                      MPI_Comm comm ):
  BufferedHistoryIO(),
  m_filename( filename ),
  m_overalloc_multiple( init_alloc ),
  m_global_idx_offset( 0 ),
  m_global_idx_count( 0 ),
  m_write_limit( init_alloc ),
  m_write_head( write_head ),
  m_hdf_type( GetHDFDataType( type_id )),
  m_type_size( H5Tget_size( m_hdf_type )),
  m_type_count( 1 ),
  m_rank( LvArray::integerConversion< hsize_t >( rank )),
  m_dims( rank ),
  m_name( name ),
  m_comm( comm ),
  m_subcomm( MPI_COMM_NULL )
{
  for( hsize_t dd = 0; dd < m_rank; ++dd )
  {
    m_dims[dd] = LvArray::integerConversion< hsize_t >( dims[dd] );
    m_type_count *= m_dims[dd];
  }
  m_data_buffer.resize( init_alloc * m_type_size * m_type_count );
}

void HDFHistIO::Init( bool exists_okay )
{
  globalIndex local_idx_count = LvArray::integerConversion< globalIndex >( m_dims[0] );

  int size = MpiWrapper::Comm_size( m_comm );
  int rank = MpiWrapper::Comm_rank( m_comm );

  int color = MPI_UNDEFINED;
  if( local_idx_count > 0 )
  {
    color = 0;
  }

  std::vector< globalIndex > counts( size );
  MpiWrapper::Allgather( &local_idx_count, 1, &counts[0], 1, m_comm );

  hsize_t min_idx_count = std::numeric_limits< hsize_t >::max( );
  int key = 0;
  for( int ii = 0; ii < size; ii++ )
  {
    if( counts[ii] != 0 && ii < rank )
    {
      key++;
    }
    if( ii == rank )
    {
      m_global_idx_offset = m_global_idx_count;
    }
    m_global_idx_count += counts[ii];
    if( counts[ii] > 0 &&  LvArray::integerConversion< hsize_t >( counts[ii] ) < min_idx_count )
    {
      min_idx_count = counts[ii];
    }
  }

  m_subcomm = MpiWrapper::Comm_split( m_comm, color, key );
  // create a dataset in the file if needed, don't erase file
  if( m_subcomm != MPI_COMM_NULL )
  {
    std::vector< hsize_t > history_file_dims( m_rank+1 );
    history_file_dims[0] = LvArray::integerConversion< hsize_t >( m_write_limit );

    std::vector< hsize_t > dim_chunks( m_rank+1 );
    dim_chunks[0] = 1;

    for( hsize_t dd = 1; dd < m_rank+1; ++dd )
    {
      // a process with chunk size 0 is considered incorrect by hdf5, hence the subcomm
      dim_chunks[dd] = m_dims[dd-1];
      history_file_dims[dd] = m_dims[dd-1];
    }
    dim_chunks[1] = min_idx_count;
    history_file_dims[1] = LvArray::integerConversion< hsize_t >( m_global_idx_count );

    // why is this killing things, only in this context, only on lassen?
    HDFFile target( m_filename, false, true, m_comm );
    //
    bool in_target = target.CheckInTarget( m_name );
    if( !in_target )
    {
      hid_t dcpl_id = 0;
      std::vector< hsize_t > max_file_dims( history_file_dims );
      // chunking is required to create an extensible dataset
      dcpl_id = H5Pcreate( H5P_DATASET_CREATE );
      H5Pset_chunk( dcpl_id, m_rank+1, &dim_chunks[0] );
      max_file_dims[0] = H5S_UNLIMITED;
      hid_t space = H5Screate_simple( m_rank+1, &history_file_dims[0], &max_file_dims[0] );
      hid_t dataset = H5Dcreate( target, m_name.c_str(), m_hdf_type, space, H5P_DEFAULT, dcpl_id, H5P_DEFAULT );
      H5Dclose( dataset );
      H5Sclose( space );
    }
    else if( exists_okay )
    {
      // todo:
      // hid_t dataset = H5Dopen(target, m_name.c_str( ), H5P_DEFAULT);
      // check that the extent of the filespace is compatible with the data
    }
    GEOSX_ERROR_IF( in_target && !exists_okay, "Dataset (" + m_name + ") already exists in output file: " + m_filename );
  }
}

void HDFHistIO::Write( )
{
  // don't need to write if nothing is buffered, this should only happen if the output event occurs before the collection event
  if( m_subcomm != MPI_COMM_NULL )
  {
    localIndex max_buffered = 0;
    MpiWrapper::allReduce( &m_buffered_count, &max_buffered, 1, MPI_MAX, m_comm );
    if( max_buffered > 0 )
    {
      ResizeFileIfNeeded( max_buffered );
      HDFFile target( m_filename, false, true, m_comm );

      hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
      hid_t filespace = H5Dget_space( dataset );

      std::vector< hsize_t > file_offset( m_rank+1 );
      file_offset[0] = LvArray::integerConversion< hsize_t >( m_write_head );
      file_offset[1] = LvArray::integerConversion< hsize_t >( m_global_idx_offset );

      std::vector< hsize_t > buffered_counts( m_rank+1 );
      buffered_counts[0] = LvArray::integerConversion< hsize_t >( max_buffered );
      for( hsize_t dd = 1; dd < m_rank+1; ++dd )
      {
        buffered_counts[dd] = m_dims[dd-1];
      }
      hid_t memspace = H5Screate_simple( m_rank+1, &buffered_counts[0], nullptr );

      hid_t file_hyperslab = filespace;
      H5Sselect_hyperslab( file_hyperslab, H5S_SELECT_SET, &file_offset[0], nullptr, &buffered_counts[0], nullptr );

      buffer_unit_type * data_buffer = nullptr;
      // if local rank is writting nothing, m_data_buffer is never alloc'd so don't try to access it
      if( m_type_count != 0 )
      {
        data_buffer = &m_data_buffer[0];
      }
      H5Dwrite( dataset, m_hdf_type, memspace, file_hyperslab, H5P_DEFAULT, data_buffer );

      H5Sclose( memspace );
      H5Sclose( filespace );
      H5Dclose( dataset );

      m_write_head += max_buffered;
      EmptyBuffer( );
    }
  }
}

void HDFHistIO::CompressInFile( )
{
  if( m_subcomm != MPI_COMM_NULL )
  {
    HDFFile target( m_filename, false, true, m_comm );
    std::vector< hsize_t > max_file_dims( m_rank+1 );
    max_file_dims[0] = LvArray::integerConversion< hsize_t >( m_write_head );
    max_file_dims[1] = LvArray::integerConversion< hsize_t >( m_global_idx_count );
    for( hsize_t dd = 2; dd < m_rank+1; ++dd )
    {
      max_file_dims[dd] = m_dims[dd-1];
    }
    hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
    H5Dset_extent( dataset, &max_file_dims[0] );
    H5Dclose( dataset );
    m_write_limit = m_write_head;
  }
}

inline void HDFHistIO::ResizeFileIfNeeded( localIndex buffered_count )
{
  if( m_subcomm != MPI_COMM_NULL )
  {
    HDFFile target( m_filename, false, true, m_comm );
    if( m_write_head + buffered_count > m_write_limit )
    {
      while( m_write_head + buffered_count > m_write_limit )
      {
        m_write_limit *= m_overalloc_multiple;
      }
      std::vector< hsize_t > max_file_dims( m_rank+1 );
      max_file_dims[0] = LvArray::integerConversion< hsize_t >( m_write_limit );
      max_file_dims[1] = LvArray::integerConversion< hsize_t >( m_global_idx_count );
      for( hsize_t dd = 2; dd < m_rank+1; ++dd )
      {
        max_file_dims[dd] = m_dims[dd-1];
      }
      hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
      H5Dset_extent( dataset, &max_file_dims[0] );
      H5Dclose( dataset );
    }
  }
}

globalIndex HDFHistIO::GetRankOffset( )
{
  return m_global_idx_offset;
}

void HDFHistIO::resizeBuffer( )
{
  size_t osize = m_data_buffer.size();
  // if needed, resize the buffer
  if( (m_buffered_count + 1) * (m_type_count * m_type_size) > osize )
  {
    m_data_buffer.resize( osize + ( m_type_count * m_type_size ) * m_overalloc_multiple );
  }
  // advance the buffer head
  m_buffer_head = (&m_data_buffer[0]) + m_buffered_count * ( m_type_count * m_type_size );
}




HDFSerialHistIO::HDFSerialHistIO( string const & filename,
                      localIndex rank,
                      const localIndex * dims,
                      string const & name,
                      std::type_index type_id,
                      localIndex write_head,
                      localIndex init_alloc,
                      MPI_Comm comm ):
  BufferedHistoryIO(),
  m_filename( filename ),
  m_overalloc_multiple( init_alloc ),
  m_write_limit( init_alloc ),
  m_write_head( write_head ),
  m_hdf_type( GetHDFDataType( type_id )),
  m_type_size( H5Tget_size( m_hdf_type )),
  m_type_count( 1 ),
  m_rank( LvArray::integerConversion< hsize_t >( rank )),
  m_dims( rank ),
  m_name( name ),
  m_comm( comm )
{
  for( hsize_t dd = 0; dd < m_rank; ++dd )
  {
    m_dims[dd] = LvArray::integerConversion< hsize_t >( dims[dd] );
    m_type_count *= m_dims[dd];
  }
  m_data_buffer.resize( init_alloc * m_type_size * m_type_count );
}

void HDFSerialHistIO::Init( bool exists_okay )
{
  // create a dataset in the file if needed, don't erase file
  if( m_type_count > 0 )
  {
    std::vector< hsize_t > history_file_dims( m_rank+1 );
    history_file_dims[0] = LvArray::integerConversion< hsize_t >( m_write_limit );

    std::vector< hsize_t > dim_chunks( m_rank+1 );
    dim_chunks[0] = 1;

    for( hsize_t dd = 1; dd < m_rank+1; ++dd )
    {
      // a process with chunk size 0 is considered incorrect by hdf5, hence the subcomm
      dim_chunks[dd] = m_dims[dd-1];
      history_file_dims[dd] = m_dims[dd-1];
    }

    // why is this killing things, only in this context, only on lassen?
    HDFFile target( m_filename, false, false, m_comm );
    //
    bool in_target = target.CheckInTarget( m_name );
    if( !in_target )
    {
      hid_t dcpl_id = 0;
      std::vector< hsize_t > max_file_dims( history_file_dims );
      // chunking is required to create an extensible dataset
      dcpl_id = H5Pcreate( H5P_DATASET_CREATE );
      H5Pset_chunk( dcpl_id, m_rank+1, &dim_chunks[0] );
      max_file_dims[0] = H5S_UNLIMITED;
      hid_t space = H5Screate_simple( m_rank+1, &history_file_dims[0], &max_file_dims[0] );
      hid_t dataset = H5Dcreate( target, m_name.c_str(), m_hdf_type, space, H5P_DEFAULT, dcpl_id, H5P_DEFAULT );
      H5Dclose( dataset );
      H5Sclose( space );
    }
    else if( exists_okay )
    {
      // todo:
      // hid_t dataset = H5Dopen(target, m_name.c_str( ), H5P_DEFAULT);
      // check that the extent of the filespace is compatible with the data
    }
    GEOSX_ERROR_IF( in_target && !exists_okay, "Dataset (" + m_name + ") already exists in output file: " + m_filename );
  }
}

void HDFSerialHistIO::Write( )
{
  // don't need to write if nothing is buffered, this should only happen if the output event occurs before the collection event
  if( m_type_count > 0 && m_buffered_count > 0)
  {
    ResizeFileIfNeeded( m_buffered_count );
    HDFFile target( m_filename, false, false, m_comm );

    hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
    hid_t filespace = H5Dget_space( dataset );

    std::vector< hsize_t > file_offset( m_rank+1 );
    file_offset[0] = LvArray::integerConversion< hsize_t >( m_write_head );
    file_offset[1] = 0;

    std::vector< hsize_t > buffered_counts( m_rank+1 );
    buffered_counts[0] = LvArray::integerConversion< hsize_t >( m_buffered_count );
    for( hsize_t dd = 1; dd < m_rank+1; ++dd )
    {
      buffered_counts[dd] = m_dims[dd-1];
    }
    hid_t memspace = H5Screate_simple( m_rank+1, &buffered_counts[0], nullptr );

    hid_t file_hyperslab = filespace;
    H5Sselect_hyperslab( file_hyperslab, H5S_SELECT_SET, &file_offset[0], nullptr, &buffered_counts[0], nullptr );

    buffer_unit_type * data_buffer = nullptr;
    // if local rank is writting nothing, m_data_buffer is never alloc'd so don't try to access it
    if( m_type_count != 0 )
    {
      data_buffer = &m_data_buffer[0];
    }
    H5Dwrite( dataset, m_hdf_type, memspace, file_hyperslab, H5P_DEFAULT, data_buffer );

    H5Sclose( memspace );
    H5Sclose( filespace );
    H5Dclose( dataset );

    m_write_head += m_buffered_count;
    EmptyBuffer( );
  }
}

void HDFSerialHistIO::CompressInFile( )
{
  if( m_type_count > 0 )
  {
    HDFFile target( m_filename, false, false, m_comm );
    std::vector< hsize_t > max_file_dims( m_rank+1 );
    max_file_dims[0] = LvArray::integerConversion< hsize_t >( m_write_head );
    max_file_dims[1] = m_dims[0];
    for( hsize_t dd = 2; dd < m_rank+1; ++dd )
    {
      max_file_dims[dd] = m_dims[dd-1];
    }
    hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
    H5Dset_extent( dataset, &max_file_dims[0] );
    H5Dclose( dataset );
    m_write_limit = m_write_head;
  }
}

inline void HDFSerialHistIO::ResizeFileIfNeeded( localIndex buffered_count )
{
  if( m_type_count > 0 )
  {
    HDFFile target( m_filename, false, false, m_comm );
    if( m_write_head + buffered_count > m_write_limit )
    {
      while( m_write_head + buffered_count > m_write_limit )
      {
        m_write_limit *= m_overalloc_multiple;
      }
      std::vector< hsize_t > max_file_dims( m_rank+1 );
      max_file_dims[0] = LvArray::integerConversion< hsize_t >( m_write_limit );
      max_file_dims[1] = m_dims[0];
      for( hsize_t dd = 2; dd < m_rank+1; ++dd )
      {
        max_file_dims[dd] = m_dims[dd-1];
      }
      hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
      H5Dset_extent( dataset, &max_file_dims[0] );
      H5Dclose( dataset );
    }
  }
}

globalIndex HDFSerialHistIO::GetRankOffset( )
{
  return 0;
}

void HDFSerialHistIO::resizeBuffer( )
{
  size_t osize = m_data_buffer.size();
  // if needed, resize the buffer
  if( (m_buffered_count + 1) * (m_type_count * m_type_size) > osize )
  {
    m_data_buffer.resize( osize + ( m_type_count * m_type_size ) * m_overalloc_multiple );
  }
  // advance the buffer head
  m_buffer_head = (&m_data_buffer[0]) + m_buffered_count * ( m_type_count * m_type_size );
}

}
