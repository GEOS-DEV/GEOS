#ifndef GEOSX_HDFFILE_HPP_
#define GEOSX_HDFFILE_HPP_

#include "cxx-utilities/src/Array.hpp"
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
  else if ( type == std::type_index(typeid(localIndex)) )
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
  HDFFile(string const & fnm) :
    filename(fnm),
    file_id(0),
    fapl_id(0),
    dxpl_id(0)
  {
    // check if file already exists
    htri_t exists = 0;
    H5E_BEGIN_TRY {
      exists = H5Fis_hdf5(filename.c_str() );
    } H5E_END_TRY
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_GEOSX, MPI_INFO_NULL);
    if( exists > 0 )
    {
      // file create access property
      file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl_id);
    }
    else if ( exists < 0 )
    {
      // this will fail if the file exists already
      file_id = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, fapl_id);
    }
    GEOSX_ERROR_IF( exists == 0, string("Existing file ") + fnm + string(" is not and HDF5 final, cannot use for HDF5 output.") );

    // use indepent parallel io to access the file
    dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id,H5FD_MPIO_INDEPENDENT);
  }
  ~HDFFile()
  {
    H5Pclose(dxpl_id);
    H5Pclose(fapl_id);
    H5Fclose(file_id);
  }
  virtual operator hid_t() final { return file_id; }
private:
  string filename;
  hid_t file_id;
  hid_t fapl_id;
  hid_t dxpl_id;
};

class HDFHistIO : public BufferedHistoryIO
{
  public:
  HDFHistIO( string const & filename,
             localIndex rank,
             const localIndex * dims,
             string const & name,
             std::type_index type_id,
             localIndex init_alloc = 4 ) :
    BufferedHistoryIO(),
    m_filename(filename),
    m_overalloc_multiple(init_alloc),
    m_global_idx_offset(0),
    m_global_idx_count(0),
    m_write_limit(init_alloc),
    m_write_head(0),
    m_hdf_type(GetHDFDataType(type_id)),
    m_type_size(H5Tget_size(m_hdf_type)),
    m_type_count(1),
    m_rank(integer_conversion<hsize_t>(rank)),
    m_dims(rank),
    m_name(name)
  {
    for(hsize_t dd = 0; dd < m_rank; ++dd)
    {
      m_dims[dd] = integer_conversion<hsize_t>( dims[dd] );
      m_type_count *= m_dims[dd];
    }
  }

  HDFHistIO( string const & filename, const HistoryMetadata & spec ) :
    HDFHistIO( filename, spec.getRank(), spec.getDims(), spec.getName(), spec.getType() )
  { }

  virtual void Init( bool exists_okay ) override
  {
    array1d<hsize_t> history_file_dims(m_rank+1);
    history_file_dims[0] = integer_conversion<hsize_t>(m_write_limit);

    array1d<hsize_t> dim_chunks(m_rank+1);
    dim_chunks[0] = 1;

    // should verify that all dims other than 0 and 1 are identical across the comm
    for(hsize_t dd = 1; dd < m_rank+1; ++dd)
    {
      dim_chunks[dd] = history_file_dims[dd] = m_dims[dd-1];
    }

    globalIndex local_idx_count = integer_conversion<globalIndex>(m_dims[0]);

    MpiWrapper::allReduce(&local_idx_count,&m_global_idx_count,1,MPI_SUM,MPI_COMM_GEOSX);
    MpiWrapper::exscan(&local_idx_count,&m_global_idx_offset,1,MPI_SUM,MPI_COMM_GEOSX);

    history_file_dims[1] = integer_conversion<hsize_t>(m_global_idx_count);

    // create a dataset in the file if needed
    HDFFile target( m_filename );
    bool in_target = target.CheckInTarget( m_name );
    if( !in_target )
    {
      hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_chunk(dcpl_id, m_rank+1, &dim_chunks[0]);

      array1d<hsize_t> max_file_dims(history_file_dims);
      max_file_dims[0] = H5S_UNLIMITED;

      hid_t space = H5Screate_simple(m_rank+1,&history_file_dims[0],&max_file_dims[0]);
      hid_t dataset = H5Dcreate(target,m_name.c_str(),m_hdf_type,space,H5P_DEFAULT,dcpl_id,H5P_DEFAULT);

      H5Dclose(dataset);
      H5Sclose(space);
    }
    else if ( exists_okay )
    {
      // hid_t dataset = H5Dopen(target, m_name.c_str( ), H5P_DEFAULT);
      // check that the extent of the filespace is compatible with the data
    }
    GEOSX_ERROR_IF( in_target && !exists_okay, "Dataset (" + m_name + ") already exists in output file: " + m_filename );
  }

  virtual void Write( ) override
  {
    ResizeFileIfNeeded( );
    HDFFile target( m_filename );

    hid_t dataset = H5Dopen(target, m_name.c_str(), H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dataset);

    array1d<hsize_t> file_offset(m_rank+1);
    forValuesInSlice(file_offset.toSlice(),[]( hsize_t & val ) { val = 0; });
    file_offset[0] = integer_conversion<hsize_t>(m_write_head);
    file_offset[1] = integer_conversion<hsize_t>(m_global_idx_offset);

    array1d<hsize_t> buffered_counts(m_rank+1);
    buffered_counts[0] = integer_conversion<hsize_t>(m_buffered_count);
    for( hsize_t dd = 1; dd < m_rank+1; ++dd )
    {
      buffered_counts[dd] = m_dims[dd-1];
    }
    hid_t memspace = H5Screate_simple(m_rank+1,&buffered_counts[0],nullptr);

    hid_t file_hyperslab = filespace;
    H5Sselect_hyperslab(file_hyperslab, H5S_SELECT_SET, &file_offset[0], nullptr, &buffered_counts[0], nullptr);

    H5Dwrite(dataset,m_hdf_type,memspace,file_hyperslab,H5P_DEFAULT,&m_data_buffer[0]);

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dataset);

    m_write_head += m_buffered_count;
  }

  inline void ResizeFileIfNeeded( )
  {
    HDFFile target( m_filename );
    if( m_write_head + m_buffered_count > m_write_limit )
    {
      while( m_write_head + m_buffered_count > m_write_limit )
      {
        m_write_limit *= m_overalloc_multiple;
      }
      array1d<hsize_t> max_file_dims(m_rank+1);
      max_file_dims[0] = integer_conversion<hsize_t>(m_write_limit);
      max_file_dims[1] = integer_conversion<hsize_t>(m_global_idx_count);
      for( hsize_t dd = 2; dd < m_rank+1; ++dd )
      {
        max_file_dims[dd] = m_dims[dd-1];
      }

      hid_t dataset = H5Dopen(target, m_name.c_str(), H5P_DEFAULT);

      H5Dset_extent(dataset,&max_file_dims[0]);
      H5Dclose(dataset);
    }
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
    hsize_t m_type_count; // prod(dims[0..n])
    hsize_t m_rank;
    array1d<hsize_t> m_dims;
    string m_name;
};


// // this and the spec likely also need to be added into the data repo to some extent
// class HDFTableIO : public BufferedHistoryIO
// {
// public:
//   HDFTableIO( ) : BufferedHistoryIO( ) { }

//   virtual void Init( string const & target_name, DataSpec * spec, bool exists_okay ) override
//   {
//     HDFFile target( target_name );
//     localIndex data_count = spec->getDiscreteDataCount();

//     std::vector<const char*> data_name_ptrs(data_count);
//     std::vector<hid_t> hdf_data_types(data_count);

//     std::type_index const * data_types = spec->getDataTypes();
//     string const * data_names = spec->getDataNames();
//     size_t const * data_sizes = spec->getDataSizes();
//     size_t const * type_sizes = spec->getTypeSizes();
//     size_t const * data_counts = spec->getDataSubcounts();

//     for( localIndex col = 0; col < data_count; ++col )
//     {
//       hsize_t dims[2] = { integer_conversion<hsize_t>(data_counts[col]), integer_conversion<hsize_t>(data_sizes[col] / ( data_counts[col] * type_sizes[col] )) };
//       hdf_data_types[col] =  dims[1] == 1 ? GetHDFArrayDataType(data_types[col],1,&dims[0]) : GetHDFArrayDataType(data_types[col],2,&dims[0]);
//       data_name_ptrs[col] = data_names[col].c_str();
//     }

//     bool in_target = CheckInTarget( target, spec );

//     if ( !in_target )
//     {
//       H5TBmake_table(spec->getTitle().c_str(),
//                      target,
//                      spec->getID().c_str(),
//                      spec->getDiscreteDataCount(),
//                      0,
//                      spec->getTotalDataSize(),
//                      &data_name_ptrs[0],
//                      spec->getDataOffsets(),
//                      &hdf_data_types[0],
//                      40,
//                      nullptr,
//                      0,
//                      nullptr);
//     }
//     else if ( in_target && !exists_okay )
//     {
//       GEOSX_ERROR( "HDFTableIO: A table with the same hdf_id already exists in the write target!");
//     }
//     else
//     {
//       GEOSX_ERROR_IF( ! CheckCompatible( target, spec ), "HDFTableIO: A table with the same hdf_id already exists in the write target, but is not compatible with the specification.");
//     }
//   }

//   virtual void Write( string const & target_name, DataSpec const * spec ) override
//   {
//     HDFFile target( target_name );
//     // MPI::Reduce(m_need_file_realloc)
//     // if (m_need_file_realloc)
//     // m_target_row_limit *= 2;
//     // H5TBreserve(m_target_row_limit)
//     // if ( do_verify ) hdf_tbl->Verify( target );
//     H5TBappend_records(target,spec->getID().c_str(),m_buffered_count,spec->getTotalDataSize(),spec->getDataOffsets(),spec->getDataSizes(),&m_data_buffer[0]);
//     EmptyBuffer();
//   }

//   virtual void ClearAfter( string const & target_name, DataSpec const * spec, localIndex last_good ) override
//   {
//     HDFFile target( target_name );
//     hsize_t num_cols = 0;
//     hsize_t num_rows = 0;
//     char const * hdf_id = spec->getID().c_str();
//     H5TBget_table_info(target,hdf_id,&num_cols,&num_rows);
//     H5TBdelete_record(target,hdf_id,last_good, num_rows - last_good);
//   }

//   inline void Verify( HDFTarget & target, DataSpec const * spec) const
//   {
//     GEOSX_ERROR_IF( ! (CheckInTarget(target, spec) && CheckCompatible(target, spec)), "HDFTable: Compatible table not found in the write target. Make sure to CreateInTarget().");
//   }

//   inline bool CheckInTarget( HDFTarget & target, DataSpec const * spec ) const
//   {
//     htri_t exists = 0;
//     H5E_BEGIN_TRY {
//       exists = H5Gget_objinfo(target, spec->getID().c_str(), 0, NULL);
//     } H5E_END_TRY
//     return ( exists  == 0 );
//   }

//   inline bool CheckCompatible( HDFTarget & target, DataSpec const * spec ) const
//   {
//     hsize_t o_col_count = 0;
//     char const * hdf_id = spec->getID().c_str();
//     H5TBget_table_info(target,hdf_id,&o_col_count,NULL);
//     if ( integer_conversion<localIndex>(o_col_count) != spec->getDiscreteDataCount() ) return false;
//     std::vector<size_t> o_col_sizes(o_col_count);
//     H5TBget_field_info(target,hdf_id,NULL,&o_col_sizes[0],NULL,NULL);
//     size_t const * col_sizes = spec->getDataSizes();
//     for( size_t col = 0; col < o_col_count; ++col)
//     {
//       if( o_col_sizes[col] != col_sizes[col] ) return false;
//     }
//     return true;
//   }
// };



}

#endif