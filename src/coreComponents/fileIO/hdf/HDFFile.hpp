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
constexpr bool can_hdf_io = std::is_same<std::remove_const_t<T>, real32>::value ||
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
    if( H5Fis_hdf5(filename.c_str() ) < 0 )
    {
      file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else
    {
      // this will fail if the file exists already
      file_id = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    }
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

namespace impl
{
  template < typename OITER, typename LAMBDA >
  void GenColNames(OITER out, localIndex const num_cols, LAMBDA && title_gen, localIndex const * const col_idxs = NULL )
  {
    if ( col_idxs != NULL )
    {
      for( localIndex lidx = 0; lidx < num_cols; ++lidx )
      {
        *out++ = title_gen(col_idxs[lidx]);
      }
    }
    else
    {
      for( localIndex lidx = 0; lidx < num_cols; ++lidx )
      {
        *out++ = title_gen(lidx);
      }
    }
  }

  template < typename OITER >
  // sfinae oiter can act as an output_iterator
  void GenColOffsets(OITER out, localIndex const num_cols, localIndex const units_per_col, size_t const unit_size, localIndex const offset_base = 0 )
  {
    for( localIndex lidx = 0; lidx < num_cols; ++lidx )
    {
      *out++ = offset_base + lidx * (unit_size * units_per_col);
    }
  }

  template < typename OITER >
  // sfinae oiter can act as an output_iterator
  void GenColTypes(OITER out, localIndex const num_cols, localIndex const units_per_col, std::type_info const & tid )
  {
    hsize_t upc = integer_conversion<hsize_t>(units_per_col);
    hid_t tp = units_per_col == 1 ? GetHDFDataType(tid) : GetHDFArrayDataType(tid,1,&upc);
    for( localIndex lidx = 0; lidx < num_cols; ++lidx )
    {
      *out++ = tp;
    }
  }

}// impl namespace

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

  inline void AddCols( localIndex const num_cols,
                localIndex const units_per_col,
                size_t const unit_byte_size,
                std::type_info const & type,
                localIndex const * idxs = NULL)
  {
    AddCols(num_cols,units_per_col,unit_byte_size,type,[](localIndex idx) { return std::to_string(idx); },idxs);
  }

  template < typename LAMBDA >
  void AddCols( localIndex const num_cols,
                localIndex const units_per_col,
                size_t const unit_byte_size,
                std::type_info const & type,
                LAMBDA && title_gen,
                localIndex const * idxs = NULL)
  {
    prepColAdd(num_cols);
    impl::GenColNames(std::back_inserter(m_col_names),num_cols,title_gen,idxs);
    for( localIndex idx = 0; idx < num_cols; ++idx)
    {
      m_col_name_ptrs[m_col_count + idx] = m_col_names[m_col_count + idx].c_str();
    }
    impl::GenColOffsets(std::back_inserter(m_col_offsets),num_cols,units_per_col,unit_byte_size,m_row_size);
    impl::GenColTypes(std::back_inserter(m_col_types),num_cols,units_per_col,type);
    m_row_size += num_cols * units_per_col * unit_byte_size;
  }

  inline void Finalize()
  {
    m_is_final = true;
  }

  inline bool isFinal() { return m_is_final; }

  inline void CreateInTarget( HDFTarget & target, localIndex prealloc_rows = 10)
  {
    warn_not_final();
    H5TBmake_table(m_title.c_str(),
                   target,
                   m_hdf_id.c_str(),
                   m_col_count,
                   prealloc_rows,
                   m_row_size,
                   const_cast<const char**>(&m_col_name_ptrs[0]),
                   &m_col_offsets[0],
                   &m_col_types[0],
                   0,
                   nullptr,
                   0,
                   nullptr);
  }

  inline void Verify( HDFTarget & target )
  {
    GEOSX_ERROR_IF( ! (VerifyInTarget(target) && VerifyCompatible(target)), "HDFTable: Compatible table not found in target. Make sure to CreateInTarget().");
  }

  inline bool VerifyInTarget( HDFTarget & target )
  {
    return ( H5Gget_objinfo(target, m_hdf_id.c_str(), 0, NULL) == 0 );
  }

  inline bool VerifyCompatible( HDFTarget & target )
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

  inline size_t getRowSize() { warn_not_final(); return m_row_size;}
  inline hsize_t getColCount() { warn_not_final(); return m_col_count; }
  inline const string & getTitle() { return m_title; }
  inline const string & getHDFID() { return m_hdf_id; }

  void warn_not_final()
  {
    if(!m_is_final)
    {
      GEOSX_WARNING("HDFTable: The table being operated on has not been finalized.");
    }
  }

protected:

  void prepColAdd(localIndex num_cols)
  {
    m_col_count += num_cols;
    m_col_names.reserve(m_col_count);
    m_col_offsets.reserve(m_col_count);
    m_col_sizes.reserve(m_col_count);
    m_col_types.reserve(m_col_count);
  }

  const string m_title;
  const string m_hdf_id;

  bool m_is_final;
  size_t m_row_size;
  hsize_t m_col_count;
  std::vector<string> m_col_names;
  std::vector<const char*> m_col_name_ptrs;
  std::vector<size_t> m_col_offsets;
  std::vector<hsize_t> m_col_sizes;
  std::vector<hid_t> m_col_types;
};

class HDFTableIO
{
public:
  HDFTableIO( HDFTable const & table_spec, localIndex buffer_default = 4 )
  : m_spec(table_spec)
  , m_file_row_head(0)
  , m_active_target()
  , m_is_open(false)
  , m_buffered_count(0)
  , m_row_buffer( buffer_default * m_spec.getRowSize() )
  , m_buffered_offsets( buffer_default * m_spec.getRowSize() )
  , m_buffered_sizes( buffer_default * m_spec.getRowSize() )
  {
    if (!m_spec.isFinal()) m_spec.warn_not_final();
    // add items to the data store
  }

  void BufferRow( buffer_unit_type const * row )
  {
    /// only supporting flat arrays at the moment
    size_t row_size = m_spec.getRowSize(); //bytes
    hsize_t col_count = m_spec.getColCount(); //# cells in a row
    m_row_buffer.resize(m_row_buffer.size() + row_size);
    memcpy(&m_row_buffer[m_buffered_count*row_size],row,row_size);


    localIndex prev_row_head = col_count * (m_buffered_count-1);
    localIndex new_row_head = col_count * m_buffered_count;

    size_t osize = m_buffered_sizes.size();
    m_buffered_sizes.resize(osize + col_count);
    memcpy(&m_buffered_sizes[new_row_head],&m_buffered_sizes[prev_row_head],col_count * sizeof(decltype(m_buffered_sizes)::value_type));

    osize = m_buffered_offsets.size();
    m_buffered_offsets.resize(osize + col_count);
    memcpy(&m_buffered_offsets[new_row_head],&m_buffered_offsets[prev_row_head], col_count * sizeof(decltype(m_buffered_offsets)::value_type));

    for( hsize_t idx = 0 ; idx < col_count; ++idx )
    {
      m_buffered_offsets[new_row_head + idx] += row_size;
    }

    m_buffered_count++;
  }

  void Open( HDFTarget & target, bool do_verify = false )
  {
    if ( do_verify ) m_spec.Verify( target );
    m_is_open = true;
    m_active_target = target;
  }

  virtual void WriteBuffered( )
  {
    assert(m_is_open);
    H5TBappend_records(m_active_target,m_spec.getHDFID().c_str(),m_buffered_count,m_spec.getRowSize(),&m_buffered_offsets[0],&m_buffered_sizes[0],&m_row_buffer);
    m_file_row_head += m_buffered_count;
    EmptyBuffer();
  }

  virtual void ClearAfterWriteHead( )
  {
    //assert(m_is_open);
    hsize_t num_cols = 0;
    hsize_t num_rows = 0;
    const string hdf_id = m_spec.getHDFID();
    H5TBget_table_info(m_active_target,hdf_id.c_str(),&num_cols,&num_rows);
    H5TBdelete_record(m_active_target,hdf_id.c_str(),m_file_row_head, num_rows - m_file_row_head);
  }

  virtual void Close( )
  {
    m_is_open = false;
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

  // in the data store
  HDFTable m_spec;
  hsize_t m_file_row_head;

  // not in the data store
  HDFTarget m_active_target;
  bool m_is_open;
  localIndex m_buffered_count;
  array1d<buffer_unit_type> m_row_buffer;
  std::vector<size_t> m_buffered_offsets;
  std::vector<size_t> m_buffered_sizes;
};

template < typename ARRAY_T >
typename std::enable_if < is_array<ARRAY_T>, void >::type
SpecFromArray( HDFTable & tbl, ARRAY_T const & arr )
{
  tbl.AddCols(arr.size( ) / arr.size( 0 ),arr.size( 0 ), sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type));
}

template < typename ARRAY_T, typename LAMBDA >
typename std::enable_if < is_array<ARRAY_T>, void >::type
SpecFromArray( HDFTable & tbl, ARRAY_T const & arr, LAMBDA && title_gen )
{
  tbl.AddCols(arr.size( ) / arr.size( 0 ),arr.size( 0 ), sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type),title_gen);
}

template < typename ARRAY_T >
typename std::enable_if < is_array<ARRAY_T>, void >::type
SpecFromArrayIndices( HDFTable & tbl, ARRAY_T const & arr, localIndex const num_indices, localIndex const * indices )
{
  tbl.AddCols(num_indices, arr.size( 0 ), sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type), indices);
}

template < typename ARRAY_T, typename LAMBDA >
typename std::enable_if < is_array<ARRAY_T>, void >::type
SpecFromArrayIndices( HDFTable & tbl, ARRAY_T const & arr, LAMBDA && title_gen, localIndex const num_indices, localIndex const * indices )
{
  tbl.AddCols(num_indices, arr.size( 0 ), sizeof(typename ARRAY_T::value_type), typeid(typename ARRAY_T::value_type),title_gen, indices);
}

}

#endif