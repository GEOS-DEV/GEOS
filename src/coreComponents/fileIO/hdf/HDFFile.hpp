#ifndef GEOSX_HDFFILE_HPP_
#define GEOSX_HDFFILE_HPP_

#include "cxx-utilities/src/Array.hpp"
#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"

#include <hdf5.h>
#include <hdf5_hl.h>

#include <string>

namespace geosx
{

using namespace traits;

template < typename T >
constexpr bool can_hdf_io = std::is_same<std::remove_const_t<T>, char>::value ||
                            std::is_same<std::remove_const_t<T>, signed char>::value ||
                            std::is_same<std::remove_const_t<T>, real32>::value ||
                            std::is_same<std::remove_const_t<T>, real64>::value ||
                            std::is_same<std::remove_const_t<T>, integer>::value ||
                            std::is_same<std::remove_const_t<T>, localIndex>::value ||
                            std::is_same<std::remove_const_t<T>, globalIndex>::value;

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

inline hid_t GetHDFDataType(std::type_info const & type)
{
  if ( type == typeid(char) )
  {
    return GetHDFDataType<char>();
  }
  else if ( type == typeid(signed char) )
  {
    return GetHDFDataType<signed char>();
  }
  else if ( type == typeid(real32) )
  {
    return GetHDFDataType<real32>();
  }
  else if ( type == typeid(real64) )
  {
    return GetHDFDataType<real64>();
  }
  else if ( type == typeid(integer) )
  {
    return GetHDFDataType<integer>();
  }
  else if ( type == typeid(localIndex) )
  {
    return GetHDFDataType<localIndex>();
  }
  else if ( type == typeid(localIndex) )
  {
    return GetHDFDataType<globalIndex>();
  }
  else
  {
    return GetHDFDataType<char>();
  }
}

inline hid_t GetHDFArrayDataType(std::type_info const & type, hsize_t const rank, hsize_t const * dims)
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

class HDFTable
{
public:

  HDFTable( const string & title,
            const string & hdf_id )
  : m_title(title)
  , m_hdf_id(hdf_id)
  , m_is_final(false)
  , m_row_size(0)
  , m_col_count(0)
  , m_col_names()
  , m_col_name_ptrs()
  , m_col_offsets()
  , m_col_sizes()
  , m_col_types()
  {}

  void AddCol(localIndex const num_items,
              localIndex const units_per_item,
              size_t const unit_byte_size,
              std::type_info const & type,
              string const & name )
  {
    prepColAdd(1);
    localIndex col_size = num_items * units_per_item * unit_byte_size;
    m_col_names.push_back( name );
    m_col_name_ptrs.push_back( m_col_names.back().c_str() );
    m_col_offsets.push_back( m_row_size );
    m_col_sizes.push_back( col_size );
    hsize_t dims[2] = {integer_conversion<hsize_t>(num_items), integer_conversion<hsize_t>(units_per_item) };
    m_col_types.push_back( col_size == 1 ? GetHDFDataType(type) : GetHDFArrayDataType(type,2,&dims[0]) );
    m_row_size += col_size;
  }

  inline void Finalize()
  {
    m_is_final = true;
  }

  inline bool isFinal() { return m_is_final; }

  inline void CreateInTarget( HDFTarget & target, bool exists_okay = true )
  {
    warn_not_final();
    bool in_target = CheckInTarget( target );
    if ( !in_target )
    {
      H5TBmake_table(m_title.c_str(),
                     target,
                     m_hdf_id.c_str(),
                     m_col_count,
                     0,
                     m_row_size,
                     &m_col_name_ptrs[0],
                     &m_col_offsets[0],
                     &m_col_types[0],
                     40,
                     nullptr,
                     0,
                     nullptr);
    }
    else if ( in_target && !exists_okay )
    {
      GEOSX_ERROR( "HDFTable: A table with the same hdf_id already exists in the write target!");
    }
    else
    {
      GEOSX_ERROR_IF( ! CheckCompatible( target ), "HDFTable: A table with the same hdf_id already exists in the write target, but is not compatible with the specification.");
    }
  }

  inline void Verify( HDFTarget & target )
  {
    GEOSX_ERROR_IF( ! (CheckInTarget(target) && CheckCompatible(target)), "HDFTable: Compatible table not found in the write target. Make sure to CreateInTarget().");
  }

  inline bool CheckInTarget( HDFTarget & target )
  {
    htri_t exists = 0;
    H5E_BEGIN_TRY {
      exists = H5Gget_objinfo(target, m_hdf_id.c_str(), 0, NULL);
    } H5E_END_TRY
    return ( exists  == 0 );
  }

  inline bool CheckCompatible( HDFTarget & target )
  {
    warn_not_final();
    hsize_t o_col_count = 0;
    //hsize_t o_num_rows = 0;
    H5TBget_table_info(target,m_hdf_id.c_str(),&o_col_count,NULL);
    if ( o_col_count != m_col_count ) return false;
    std::vector<size_t> o_col_sizes(m_col_count);
    H5TBget_field_info(target,m_hdf_id.c_str(),NULL,&o_col_sizes[0],NULL,NULL);
    for( size_t col = 0; col < m_col_count; ++col)
    {
      if( o_col_sizes[col] != m_col_sizes[col] ) return false;
    }
    return true;
  }

  inline size_t getRowSize() const { warn_not_final(); return m_row_size;}
  inline hsize_t getColCount() const { warn_not_final(); return m_col_count; }
  inline size_t const * getColSizes() const { warn_not_final(); return &m_col_sizes[0]; }
  inline size_t const * getColOffsets() const { warn_not_final(); return &m_col_offsets[0]; }
  inline const string & getTitle() const { return m_title; }
  inline const string & getHDFID() const { return m_hdf_id; }

  void warn_not_final() const
  {
    if(!m_is_final)
    {
      GEOSX_WARNING("HDFTable: The table being operated on has not been finalized.");
    }
  }


  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim == 1) && can_hdf_io<typename ARRAY_T::value_type>, void >::type
  AddArrayCol( ARRAY_T const & arr, string const & name, localIndex vals_per_col = 1 )
  {
    AddCol(arr.size( ), vals_per_col, sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type), name);
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim > 1) && can_hdf_io<typename ARRAY_T::value_type>, void >::type
  AddArrayCol( ARRAY_T const & arr, string const & name)
  {
    AddCol(arr.size( ) / arr.size( 0 ), arr.size( 0 ), sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type), name);
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim == 1) && can_hdf_io<typename ARRAY_T::value_type>, void >::type
  AddArrayIndicesCol( ARRAY_T const & GEOSX_UNUSED_PARAM( arr ), string const & name, localIndex const num_indices, localIndex vals_per_col = 1 )
  {
    AddCol(num_indices, vals_per_col, sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type), name);
  }

  template < typename ARRAY_T >
  inline
  typename std::enable_if < is_array<ARRAY_T> && (ARRAY_T::ndim > 1) && can_hdf_io<typename ARRAY_T::value_type>, void >::type
  AddArrayIndicesCol( ARRAY_T const & arr, string const & name, localIndex const num_indices )
  {
    AddCol(num_indices, arr.size( 0 ), sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type), name);
  }

protected:

  void prepColAdd(localIndex num_cols)
  {
    m_col_count += num_cols;
    m_col_names.reserve(m_col_count);
    m_col_name_ptrs.reserve(m_col_count);
    m_col_offsets.reserve(m_col_count);
    m_col_sizes.reserve(m_col_count);
    m_col_types.reserve(m_col_count);
  }

  string m_title;
  string m_hdf_id;
  size_t m_prealloc_rows;

  bool m_is_final;
  size_t m_row_size;
  hsize_t m_col_count;
  std::vector<string> m_col_names;
  std::vector<const char*> m_col_name_ptrs;
  std::vector<size_t> m_col_offsets;
  std::vector<size_t> m_col_sizes;
  std::vector<hid_t> m_col_types;
};

class HDFTableIO
{
public:
  HDFTableIO( HDFTable const & table_spec, localIndex buffer_default = 4 )
  : m_spec(table_spec)
  , m_target_row_head(0)
  // , m_need_file_realloc(false)
  , m_buffered_count(0)
  , m_row_buffer( buffer_default * m_spec.getRowSize() )
  {
    if (!m_spec.isFinal()) m_spec.warn_not_final();
  }

  void BufferRow( buffer_unit_type const * row )
  {
    size_t row_size = m_spec.getRowSize(); //bytes
    m_row_buffer.resize(m_row_buffer.size() + row_size);
    memcpy(&m_row_buffer[m_buffered_count*row_size],row,row_size);
    m_buffered_count++;
    // if ( m_target_row_head + m_buffered_count > m_target_row_limit )
    // {
    //   m_need_file_realloc = true;
    // }
  }

  virtual void CreateInTarget( HDFTarget & target, bool exists_okay = true )
  {
    m_spec.CreateInTarget( target, exists_okay );
  }

  virtual void WriteBuffered( HDFTarget & target, bool do_verify = false )
  {
    // MPI::Reduce(m_need_file_realloc)
    // if (m_need_file_realloc)
    // m_target_row_limit *= 2;
    // H5TBreserve(m_target_row_limit)
    if ( do_verify ) m_spec.Verify( target );
    H5TBappend_records(target,m_spec.getHDFID().c_str(),m_buffered_count,m_spec.getRowSize(),m_spec.getColOffsets(),m_spec.getColSizes(),&m_row_buffer[0]);
    m_target_row_head += m_buffered_count;
    EmptyBuffer();
  }

  virtual void ClearAfterWriteHead( HDFTarget & target )
  {
    hsize_t num_cols = 0;
    hsize_t num_rows = 0;
    const string hdf_id = m_spec.getHDFID();
    H5TBget_table_info(target,hdf_id.c_str(),&num_cols,&num_rows);
    H5TBdelete_record(target,hdf_id.c_str(),m_target_row_head, num_rows - m_target_row_head);
  }

private:

  void EmptyBuffer( bool dealloc = false )
  {
    m_buffered_count = 0;
    if( dealloc )
    {
      m_row_buffer.resize(0);
    }
  }

  HDFTable m_spec;
  // in the data store
  size_t m_target_row_head;
  // size_t m_target_row_limit;
  // bool m_need_file_realloc;

  // not in the data store
  localIndex m_buffered_count;
  array1d<buffer_unit_type> m_row_buffer;
};

// inline HDFTable TimeSeries( HDFTable const & base_table )
// {
//   string time_title = base_table.getTitle() + string(" Time");
//   string time_id = base_table.getHDFID() + string("_t");
//   HDFTable time_table(time_title,time_id);
//   time_table.AddCol(1,1,sizeof(real64),typeid(real64),"Time");
//   time_table.Finalize();
//   return time_table;
// }

inline HDFTable InitHistoryTable( string const & title, string const & hdf_id )
{
  HDFTable time_history( title, hdf_id );
  time_history.AddCol( 1, 1, sizeof(real64), typeid(real64), "time" );
  return time_history;
}


}

#endif