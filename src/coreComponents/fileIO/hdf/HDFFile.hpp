#ifndef GEOSX_HDFFILE_HPP_
#define GEOSX_HDFFILE_HPP_

#include "LvArray/src/Array.hpp"
#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"

#include "managers/TimeHistory/HistoryDataSpec.hpp"
#include "managers/Outputs/TimeHistoryOutput.hpp"

#include "mpiCommunications/MpiWrapper.hpp"

#include <hdf5.h>
#include <hdf5_hl.h>

#include <string>

namespace geosx
{

template < typename T >
inline hid_t GetHDFDataType();
template <>
inline hid_t GetHDFDataType<char>() { return H5T_NATIVE_CHAR; }
template <>
inline hid_t GetHDFDataType<signed char>() { return H5T_NATIVE_CHAR; }
template <>
inline hid_t GetHDFDataType<real32>() { return H5T_NATIVE_FLOAT; }
template <>
inline hid_t GetHDFDataType<real64>() { return H5T_NATIVE_DOUBLE; }
template <>
inline hid_t GetHDFDataType<integer>() { return H5T_NATIVE_INT; }
template <>
inline hid_t GetHDFDataType<localIndex>() { return H5T_NATIVE_LONG; }
template <>
inline hid_t GetHDFDataType<globalIndex>() { return H5T_NATIVE_LLONG; }

inline hid_t GetHDFDataType(std::type_index const & type)
{
  if ( type == std::type_index(typeid(char)) )
  {
    return GetHDFDataType<char>();
  }
  else if ( type == std::type_index(typeid(signed char)) )
  {
    return GetHDFDataType<signed char>();
  }
  else if ( type == std::type_index(typeid(real32)) )
  {
    return GetHDFDataType<real32>();
  }
  else if ( type == std::type_index(typeid(real64)) )
  {
    return GetHDFDataType<real64>();
  }
  else if ( type == std::type_index(typeid(integer)) )
  {
    return GetHDFDataType<integer>();
  }
  else if ( type == std::type_index(typeid(localIndex)) )
  {
    return GetHDFDataType<localIndex>();
  }
  else if ( type == std::type_index(typeid(globalIndex)) )
  {
    return GetHDFDataType<globalIndex>();
  }
  else
  {
    return GetHDFDataType<char>();
  }
}

inline hid_t GetHDFArrayDataType(std::type_index const & type, hsize_t const rank, hsize_t const * dims)
{
  return H5Tarray_create(GetHDFDataType(type),rank,dims);
}

class HDFTarget
{
public:
  virtual operator hid_t() { return 0; }
  virtual bool CheckInTarget( const string & name )
  {
    htri_t exists = 0;
    H5E_BEGIN_TRY {
      exists = H5Gget_objinfo(this->operator hid_t(), name.c_str(), 0, NULL);
    } H5E_END_TRY
        return (exists == 0);
  }
};

class HDFFile : public HDFTarget
{
public:
  HDFFile(string const & fnm, bool delete_existing = false, MPI_Comm comm = MPI_COMM_GEOSX) :
    m_filename(fnm),
    m_file_id(0),
    m_fapl_id(0),
    m_dxpl_id(0),
    m_comm(comm)
  {
    // check if file already exists
    htri_t exists = 0;
    H5E_BEGIN_TRY {
      exists = H5Fis_hdf5(m_filename.c_str() );
    } H5E_END_TRY
    m_fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(m_fapl_id, m_comm, MPI_INFO_NULL);
    if( exists > 0 && !delete_existing )
    {
      m_file_id = H5Fopen(m_filename.c_str(), H5F_ACC_RDWR, m_fapl_id);
    }
    else if ( exists > 0 && delete_existing )
    {
      m_file_id = H5Fcreate(m_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, m_fapl_id);
    }
    else if ( exists < 0 )
    {
      m_file_id = H5Fcreate(m_filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, m_fapl_id);
    }
    GEOSX_ERROR_IF( exists == 0, string("Existing file ") + m_filename + string(" is not an HDF5 file, cannot use for HDF5 output.") );

    // use indepent parallel io to access the file
    m_dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(m_dxpl_id,H5FD_MPIO_INDEPENDENT);
  }
  ~HDFFile()
  {
    H5Pclose(m_dxpl_id);
    H5Pclose(m_fapl_id);
    H5Fclose(m_file_id);
  }
  virtual operator hid_t() final { return m_file_id; }
private:
  string m_filename;
  hid_t m_file_id;
  hid_t m_fapl_id;
  hid_t m_dxpl_id;
  MPI_Comm m_comm;
};

class HDFHistIO : public BufferedHistoryIO
{
  public:
  HDFHistIO( string const & filename,
             localIndex rank,
             const localIndex * dims,
             string const & name,
             std::type_index type_id,
             localIndex write_head = 0,
             localIndex init_alloc = 4,
             MPI_Comm comm = MPI_COMM_GEOSX) :
    BufferedHistoryIO(),
    m_filename(filename),
    m_overalloc_multiple(init_alloc),
    m_global_idx_offset(0),
    m_global_idx_count(0),
    m_write_limit(init_alloc),
    m_write_head(write_head),
    m_hdf_type(GetHDFDataType(type_id)),
    m_type_size(H5Tget_size(m_hdf_type)),
    m_type_count(1),
    m_rank(LvArray::integerConversion<hsize_t>(rank)),
    m_dims(rank),
    m_name(name),
    m_comm(comm)
  {
    for(hsize_t dd = 0; dd < m_rank; ++dd)
    {
      m_dims[dd] = LvArray::integerConversion<hsize_t>( dims[dd] );
      m_type_count *= m_dims[dd];
    }
  }

  HDFHistIO( string const & filename, const HistoryMetadata & spec, localIndex write_head = 0, localIndex init_alloc = 4, MPI_Comm comm = MPI_COMM_GEOSX ) :
    HDFHistIO( filename, spec.getRank(), spec.getDims(), spec.getName(), spec.getType(), write_head, init_alloc, comm )
  { }

  virtual ~HDFHistIO() { }

  virtual void Init( bool exists_okay, bool once ) override
  {
    array1d<hsize_t> history_file_dims(m_rank+1);
    history_file_dims[0] = LvArray::integerConversion<hsize_t>(m_write_limit);

    array1d<hsize_t> dim_chunks(m_rank+1);
    dim_chunks[0] = 1;

    for(hsize_t dd = 1; dd < m_rank+1; ++dd)
    {
      // a process with chunk size 0 is considered incorrect by hdf5
      dim_chunks[dd] = m_dims[dd-1];
      history_file_dims[dd] = m_dims[dd-1];
    }

    globalIndex local_idx_count = LvArray::integerConversion<globalIndex>(m_dims[0]);

    MpiWrapper::allReduce(&local_idx_count,&m_global_idx_count,1,MPI_SUM,m_comm);
    MpiWrapper::exscan(&local_idx_count,&m_global_idx_offset,1,MPI_SUM,m_comm);

    history_file_dims[1] = LvArray::integerConversion<hsize_t>(m_global_idx_count);

    // create a dataset in the file if needed, don't erase file
    HDFFile target( m_filename, false, m_comm );
    bool in_target = target.CheckInTarget( m_name );
    if( !in_target )
    {
      hid_t dcpl_id = 0;
      array1d<hsize_t> max_file_dims(history_file_dims);
      if ( ! once )
      {
        // chunking is required to create an extensible dataset
        dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(dcpl_id, m_rank+1, &dim_chunks[0]);
        max_file_dims[0] = H5S_UNLIMITED;
      }
      else
      {
        // chunking can't be used to create a dataset with zero chunk size on one processor
        //  if the data is only being written once ( this makes no sense but seems to be the case )
        dcpl_id = H5P_DEFAULT;
      }

      hid_t space = H5Screate_simple(m_rank+1,&history_file_dims[0],&max_file_dims[0]);
      hid_t dataset = H5Dcreate(target,m_name.c_str(),m_hdf_type,space,H5P_DEFAULT,dcpl_id,H5P_DEFAULT);

      H5Dclose(dataset);
      H5Sclose(space);
    }
    else if ( exists_okay )
    {
      // todo:
      // hid_t dataset = H5Dopen(target, m_name.c_str( ), H5P_DEFAULT);
      // check that the extent of the filespace is compatible with the data
    }
    GEOSX_ERROR_IF( in_target && !exists_okay, "Dataset (" + m_name + ") already exists in output file: " + m_filename );
  }

  virtual void Write( ) override
  {
    // don't need to write if nothing is buffered, this should only happen if the output event occurs before the collection event
    localIndex max_buffered = 0;
    MpiWrapper::allReduce(&m_buffered_count,&max_buffered,1,MPI_MAX,m_comm);
    if( max_buffered > 0 )
    {
      ResizeFileIfNeeded( max_buffered );
      HDFFile target( m_filename, false,  m_comm );

      hid_t dataset = H5Dopen(target, m_name.c_str(), H5P_DEFAULT);
      hid_t filespace = H5Dget_space(dataset);

      array1d<hsize_t> file_offset(m_rank+1);
      forValuesInSlice(file_offset.toSlice(),[]( hsize_t & val ) { val = 0; });
      file_offset[0] = LvArray::integerConversion<hsize_t>(m_write_head);
      file_offset[1] = LvArray::integerConversion<hsize_t>(m_global_idx_offset);

      array1d<hsize_t> buffered_counts(m_rank+1);
      buffered_counts[0] = LvArray::integerConversion<hsize_t>(max_buffered);
      for( hsize_t dd = 1; dd < m_rank+1; ++dd )
      {
        buffered_counts[dd] = m_dims[dd-1];
      }
      hid_t memspace = H5Screate_simple(m_rank+1,&buffered_counts[0],nullptr);

      hid_t file_hyperslab = filespace;
      H5Sselect_hyperslab(file_hyperslab, H5S_SELECT_SET, &file_offset[0], nullptr, &buffered_counts[0], nullptr);

      buffer_unit_type * data_buffer = nullptr;
      // if local rank is writting nothing, m_data_buffer is never alloc'd so don't try to access it
      if ( m_type_count != 0 )
      {
        data_buffer = &m_data_buffer[0];
      }
      H5Dwrite(dataset,m_hdf_type,memspace,file_hyperslab,H5P_DEFAULT,data_buffer);

      H5Sclose(memspace);
      H5Sclose(filespace);
      H5Dclose(dataset);

      m_write_head += max_buffered;
      EmptyBuffer( );
    }
  }

  virtual void CompressInFile( ) override
  {
    HDFFile target( m_filename, false, m_comm );
    array1d<hsize_t> max_file_dims(m_rank+1);
    max_file_dims[0] = LvArray::integerConversion<hsize_t>(m_write_head);
    max_file_dims[1] = LvArray::integerConversion<hsize_t>(m_global_idx_count);
    for( hsize_t dd = 2; dd < m_rank+1; ++dd )
    {
      max_file_dims[dd] = m_dims[dd-1];
    }
    hid_t dataset = H5Dopen(target, m_name.c_str(), H5P_DEFAULT);
    H5Dset_extent(dataset,&max_file_dims[0]);
    H5Dclose(dataset);
    m_write_limit = m_write_head;
  }

  inline void ResizeFileIfNeeded( localIndex buffered_count )
  {
    HDFFile target( m_filename, false, m_comm );
    if( m_write_head + buffered_count > m_write_limit )
    {
      while( m_write_head + buffered_count > m_write_limit )
      {
        m_write_limit *= m_overalloc_multiple;
      }
      array1d<hsize_t> max_file_dims(m_rank+1);
      max_file_dims[0] = LvArray::integerConversion<hsize_t>(m_write_limit);
      max_file_dims[1] = LvArray::integerConversion<hsize_t>(m_global_idx_count);
      for( hsize_t dd = 2; dd < m_rank+1; ++dd )
      {
        max_file_dims[dd] = m_dims[dd-1];
      }

      hid_t dataset = H5Dopen(target, m_name.c_str(), H5P_DEFAULT);

      H5Dset_extent(dataset,&max_file_dims[0]);
      H5Dclose(dataset);
    }
  }

  virtual globalIndex GetRankOffset( ) override
  {
    return m_global_idx_offset;
  }

  protected:
  virtual void resizeBuffer( ) override
  {
    size_t osize = m_data_buffer.size();
    m_data_buffer.resize(osize + m_type_count * m_type_size);
  }

  private:
    // file io params
    string m_filename;
    const localIndex m_overalloc_multiple;
    globalIndex m_global_idx_offset;
    globalIndex m_global_idx_count;
    localIndex m_write_limit;
    localIndex m_write_head;

    // history metadata
    hsize_t m_hdf_type;
    size_t m_type_size;
    hsize_t m_type_count; // prod(dims[0:n])
    hsize_t m_rank;
    array1d<hsize_t> m_dims;
    string m_name;

    MPI_Comm m_comm;
};

}

#endif
