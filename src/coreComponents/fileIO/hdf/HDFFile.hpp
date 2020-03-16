#ifndef GEOSX_HDFFILE_HPP_
#define GEOSX_HDFFILE_HPP_

#include "cxx-utilities/src/Array.hpp"
#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"

#include "managers/TimeHistory/HistoryDataSpec.hpp"
#include "managers/Outputs/TimeHistoryOutput.hpp"

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
};

class HDFFile : public HDFTarget
{
public:
  HDFFile(string const & fnm) :
    filename(fnm),
    file_id(0)
  {
    // check if file already exists
    htri_t exists = 0;
    H5E_BEGIN_TRY {
      exists = H5Fis_hdf5(filename.c_str() );
    } H5E_END_TRY
    if( exists > 0 )
    {
      file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else if ( exists < 0 )
    {
      // this will fail if the file exists already
      file_id = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    }
    GEOSX_ERROR_IF( exists == 0, string("Existing file ") + fnm + string(" is not and HDF5 final, cannot use for HDF5 output.") );
  }
  ~HDFFile()
  {
    H5Fclose(file_id);
  }
  virtual operator hid_t() final { return file_id; }
private:
  string filename;
  hid_t file_id;
};

// this and the spec likely also need to be added into the data repo to some extent
class HDFTableIO : public BufferedHistoryIO
{
public:
  HDFTableIO( ) : BufferedHistoryIO( ) { }

  virtual void Init( string const &  m_write_target_name, DataSpec * spec, bool exists_okay ) override
  {
    HDFFile target( m_write_target_name );
    localIndex data_count = spec->getDiscreteDataCount();

    std::vector<const char*> data_name_ptrs(data_count);
    std::vector<hid_t> hdf_data_types(data_count);

    std::type_index const * data_types = spec->getDataTypes();
    string const * data_names = spec->getDataNames();
    size_t const * data_sizes = spec->getDataSizes();
    size_t const * type_sizes = spec->getTypeSizes();
    size_t const * data_counts = spec->getDataSubcount();

    for( localIndex col = 0; col < data_count; ++col )
    {
      hsize_t dims[2] = { integer_conversion<hsize_t>(data_counts[col]), integer_conversion<hsize_t>(data_sizes[col] / ( data_counts[col] * type_sizes[col] )) };
      hdf_data_types[col] =  dims[1] == 1 ? GetHDFArrayDataType(data_types[col],1,&dims[0]) : GetHDFArrayDataType(data_types[col],2,&dims[0]);
      data_name_ptrs[col] = data_names[col].c_str();
    }

    bool in_target = CheckInTarget( target, spec );
    if ( !in_target )
    {
      H5TBmake_table(spec->getTitle().c_str(),
                     target,
                     spec->getID().c_str(),
                     spec->getDiscreteDataCount(),
                     0,
                     spec->getTotalDataSize(),
                     &data_name_ptrs[0],
                     spec->getDataOffsets(),
                     &hdf_data_types[0],
                     40,
                     nullptr,
                     0,
                     nullptr);
    }
    else if ( in_target && !exists_okay )
    {
      GEOSX_ERROR( "HDFTableIO: A table with the same hdf_id already exists in the write target!");
    }
    else
    {
      GEOSX_ERROR_IF( ! CheckCompatible( target, spec ), "HDFTableIO: A table with the same hdf_id already exists in the write target, but is not compatible with the specification.");
    }
  }

  virtual void Write( string const &  m_write_target_name, DataSpec const * spec ) override
  {
    HDFFile target( m_write_target_name );
    // MPI::Reduce(m_need_file_realloc)
    // if (m_need_file_realloc)
    // m_target_row_limit *= 2;
    // H5TBreserve(m_target_row_limit)
    // if ( do_verify ) hdf_tbl->Verify( target );
    H5TBappend_records(target,spec->getID().c_str(),m_buffered_count,spec->getTotalDataSize(),spec->getDataOffsets(),spec->getDataSizes(),&m_data_buffer[0]);
    EmptyBuffer();
  }

  virtual void ClearAfter( string const &  m_write_target_name, DataSpec const * spec, localIndex last_good ) override
  {
    HDFFile target( m_write_target_name );
    hsize_t num_cols = 0;
    hsize_t num_rows = 0;
    char const * hdf_id = spec->getID().c_str();
    H5TBget_table_info(target,hdf_id,&num_cols,&num_rows);
    H5TBdelete_record(target,hdf_id,last_good, num_rows - last_good);
  }

  inline void Verify( HDFTarget & target, DataSpec const * spec) const
  {
    GEOSX_ERROR_IF( ! (CheckInTarget(target, spec) && CheckCompatible(target, spec)), "HDFTable: Compatible table not found in the write target. Make sure to CreateInTarget().");
  }

  inline bool CheckInTarget( HDFTarget & target, DataSpec const * spec ) const
  {
    htri_t exists = 0;
    H5E_BEGIN_TRY {
      exists = H5Gget_objinfo(target, spec->getID().c_str(), 0, NULL);
    } H5E_END_TRY
    return ( exists  == 0 );
  }

  inline bool CheckCompatible( HDFTarget & target, DataSpec const * spec ) const
  {
    hsize_t o_col_count = 0;
    char const * hdf_id = spec->getID().c_str();
    H5TBget_table_info(target,hdf_id,&o_col_count,NULL);
    if ( integer_conversion<localIndex>(o_col_count) != spec->getDiscreteDataCount() ) return false;
    std::vector<size_t> o_col_sizes(o_col_count);
    H5TBget_field_info(target,hdf_id,NULL,&o_col_sizes[0],NULL,NULL);
    size_t const * col_sizes = spec->getDataSizes();
    for( size_t col = 0; col < o_col_count; ++col)
    {
      if( o_col_sizes[col] != col_sizes[col] ) return false;
    }
    return true;
  }
};



}

#endif