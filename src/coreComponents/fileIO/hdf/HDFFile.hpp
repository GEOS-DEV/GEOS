#include "cxx-utilities/src/Array.hpp"
#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"

#include <hdf5.h>
#include <hdf5_hl.h>

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
hid_t GetHDFDataType();
template <>
hid_t GetHDFDataType<real32>() { return H5T_NATIVE_FLOAT; }
template <>
hid_t GetHDFDataType<real64>() { return H5T_NATIVE_DOUBLE; }
template <>
hid_t GetHDFDataType<integer>() { return H5T_NATIVE_INT; }
template <>
hid_t GetHDFDataType<localIndex>() { return H5T_NATIVE_LONG; }
template <>
hid_t GetHDFDataType<globalIndex>() { return H5T_NATIVE_LLONG; }

template < typename T >
hid_t GetHDFArrayDataType(localIndex rnk, hsize_t * dims)
{
  return H5Tarray_create(GetHDFDataType<T>(), rnk, dims);
}

template < typename T >
hid_t GetHDFArrayDataType(hsize_t dim)
{
  return GetHDFArrayDataType<T>(1,&dim);
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
    file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }
  ~HDFFile()
  {
    H5Fclose(file_id);
  }
  // implicit conversion operator so we can use this in place of the usual hid_t
  virtual operator hid_t() final { return file_id; }
private:
  string filename;
  hid_t file_id;
};

namespace impl
{
  template < typename OITER, typename IDX_TYPE >
  // sfinae: oiter must be output iter that can accept strings, idx must be integral
  typename std::enable_if< std::is_integral<IDX_TYPE>::value, void>::type
  GenColNames(OITER out, string const & prefix, localIndex num_cols, IDX_TYPE const * col_idxs )
  {
    for( localIndex lidx = 0; lidx < num_cols; ++lidx )
    {
      *out++ = (prefix + std::to_string(col_idxs[lidx])).c_str();
    }
  }

  template < typename OITER >
  void GenColNames(OITER out, string const & prefix, localIndex num_cols )
  {
    for( localIndex lidx = 0; lidx < num_cols; ++lidx )
    {
      *out++ = (prefix + std::to_string(lidx)).c_str();
    }
  }

  template < typename T, typename OITER, typename IDX_TYPE >
  // sfinae: oiter must be output iter that can accept Ts (trivially T*), idx must be integral
  typename std::enable_if< std::is_integral<IDX_TYPE>::value, void>::type
  GenColOffsets(OITER out, localIndex type_per_cell, localIndex num_cols, IDX_TYPE const * col_idxs)
  {
    for( localIndex lidx = 0; lidx < num_cols; ++lidx )
    {
      *out++ = (col_idxs[num_cols]) * (sizeof(T) * type_per_cell);
    }
  }

  template < typename T, typename OITER >
  // sfinae: oiter must be output iter that can accept Ts (trivially T*), idx must be integral
  void GenColOffsets(OITER out, localIndex type_per_cell, localIndex num_cols )
  {
    for( localIndex lidx = 0; lidx < num_cols; ++lidx )
    {
      *out++ = lidx * (sizeof(T) * type_per_cell);
    }
  }

  template <typename T>
  size_t inline CalcRowSize( size_t unit_per_cell, size_t cells_per_row )
  {
    return sizeof(T) * unit_per_cell * cells_per_row;
  }

  bool TryOpen( hid_t target, string const & id )
  {
    return ( H5Gget_objinfo(target, id.c_str(), 0, NULL) == 0 );
  }

  bool VerifyTable( hid_t target, string const & id, hsize_t num_cols, hsize_t cell_size )
  {
    hsize_t o_num_cols = 0;
    hsize_t o_num_rows = 0;
    H5TBget_table_info(target,id.c_str(),&o_num_cols,&o_num_rows);
    if( o_num_cols != num_cols ) return false;
    std::vector<size_t> col_sizes(num_cols);
    H5TBget_field_info(target,id.c_str(),NULL,&col_sizes[0],NULL,NULL);
    for( auto cs : col_sizes )
    {
      if( cs != cell_size ) return false;
    }
    return true;
  }

  template < typename DATA_TYPE >
  typename std::enable_if< can_hdf_io<DATA_TYPE>, void>::type
  CreateTable(hid_t target,
              string const & title,
              string const & id,
              localIndex const unit_per_cell,
              localIndex const num_cols,
              string const & record_prefix = "")
  {
    size_t row_size = CalcRowSize<DATA_TYPE>(unit_per_cell,num_cols);
    std::vector<const char*> col_names(num_cols);
    GenColNames(std::back_insert_iterator<decltype(col_names)>(col_names), record_prefix, num_cols);
    std::vector<size_t> col_offsets(num_cols);
    GenColOffsets<DATA_TYPE>(std::back_insert_iterator<decltype(col_offsets)>(col_offsets), unit_per_cell, num_cols);
    hid_t hdf_data_type = unit_per_cell == 1 ? GetHDFDataType<DATA_TYPE>() : GetHDFArrayDataType<DATA_TYPE>(unit_per_cell);
    std::vector<hid_t> col_types(num_cols, hdf_data_type);
    H5TBmake_table(title.c_str(), target, id.c_str(), num_cols, 0, row_size, &col_names[0], &col_offsets[0], &col_types[0], 0, nullptr, 0, nullptr);
  }

  template < typename DATA_TYPE, typename IDX_TYPE >
  typename std::enable_if< can_hdf_io<DATA_TYPE> && std::is_integral<IDX_TYPE>::value, void>::type
  CreateTable(hid_t target,
              string const & title,
              string const & id,
              localIndex const unit_per_cell,
              localIndex const num_cols,
              IDX_TYPE const * idxs,
              string const & record_prefix = "")
  {
    size_t row_size = CalcRowSize<DATA_TYPE>(unit_per_cell,num_cols);
    std::vector<const char*> col_names(num_cols);
    GenColNames(std::back_insert_iterator<decltype(col_names)>(col_names), record_prefix, num_cols, idxs);
    std::vector<size_t> col_offsets(num_cols);
    GenColOffsets<DATA_TYPE>(std::back_insert_iterator<decltype(col_offsets)>(col_offsets), unit_per_cell, num_cols, idxs);
    hid_t hdf_data_type = unit_per_cell == 1 ? GetHDFDataType<DATA_TYPE>() : GetHDFArrayDataType<DATA_TYPE>(unit_per_cell);
    std::vector<hid_t> col_types(num_cols, hdf_data_type);
    H5TBmake_table(title.c_str(), target, id.c_str(), num_cols, 0, row_size, &col_names[0], &col_offsets[0], &col_types[0], 0, nullptr, 0, nullptr);
  }

  template < typename DATA_TYPE >
  void AppendRow( hid_t target,
                  string const & id,
                  localIndex const unit_per_cell,
                  localIndex const num_cols,
                  DATA_TYPE const * data )
  {
    size_t row_size = CalcRowSize<DATA_TYPE>(unit_per_cell,num_cols);
    std::vector<size_t> col_offsets(num_cols);
    GenColOffsets<DATA_TYPE>(std::back_insert_iterator<decltype(col_offsets)>(col_offsets), unit_per_cell, num_cols);
    std::vector<size_t> col_sizes(num_cols,sizeof(DATA_TYPE));
    H5TBappend_records(target,id.c_str(),1,row_size,&col_offsets[0],&col_sizes[0],data);
  }

  template < typename DATA_TYPE, typename IDX_TYPE >
  void AppendRow( hid_t target,
                  string const & id,
                  localIndex const unit_per_cell,
                  localIndex const num_cols,
                  DATA_TYPE const * data,
                  IDX_TYPE const * idxs = nullptr)
  {
    size_t row_size = CalcRowSize<DATA_TYPE>(unit_per_cell,num_cols);
    std::vector<size_t> col_offsets(num_cols);
    GenColOffsets<DATA_TYPE>(std::back_insert_iterator<decltype(col_offsets)>(col_offsets), unit_per_cell, num_cols, idxs);
    std::vector<size_t> col_sizes(num_cols,sizeof(DATA_TYPE));
    H5TBappend_records(target,id.c_str(),1,row_size,&col_offsets[0],&col_sizes[0],data);
  }

  template < typename DATA_TYPE, typename COND > //, typename IDX_TYPE = localIndex >
  // sfinae lamda callable, two input args, bool return
  localIndex ColSearch( hid_t target,
                        string const & id,
                        integer const col_idx,
                        localIndex const unit_per_cell,
                        COND && cond)
                        //IDX_TYPE const * idxs = nullptr )
  {
    hsize_t num_cols = 0;
    hsize_t num_rows = 0;
    H5TBget_table_info(target,id.c_str(),&num_cols,&num_rows);
    size_t row_size = CalcRowSize<DATA_TYPE>(unit_per_cell,num_cols);
    std::vector<size_t> col_offsets(num_cols);
    GenColOffsets<DATA_TYPE>(std::back_insert_iterator<decltype(col_offsets)>(col_offsets), unit_per_cell, num_cols );
    size_t cell_size = sizeof(DATA_TYPE) * unit_per_cell;
    Array<DATA_TYPE,1> col_data(num_rows);
    // the offsets are contiguous in this case regardless of whether the table was created from an indexed array or whole array,
    //  the hope is that since this is a read, we can simply ignore the typical "correct" specification of the offsets and read in
    //  the column contigously regardless.. might not work in practice, needs testing (and the hdf spec doesn't make it clear why the
    //  argument is even required since it is supposed to specify the offsets into the passed in data afaik, not the offests into table)
    H5TBread_fields_index(target,id.c_str(),1,&col_idx,0,num_rows,row_size,&col_offsets[col_idx],&cell_size,col_data.data());
    hsize_t row = 0;
    bool found = false;
    for(; row < num_rows; ++row)
    {
      if ( cond(&col_data[row]) )
      {
        found = true;
        break;
      }
    }
    return found ? static_cast<localIndex>(row) : -1;
  }

  void ClearAfter( hid_t target,
                   string const & id,
                   localIndex const row_idx )
  {
    hsize_t num_cols = 0;
    hsize_t num_rows = 0;
    H5TBget_table_info(target,id.c_str(),&num_cols,&num_rows);
    H5TBdelete_record(target,id.c_str(),row_idx,num_rows-row_idx);
  }
}

class HDFIO
{
public:
  HDFIO(HDFTarget & target) : io_target(target) {}
  virtual ~HDFIO() {}
protected:
  HDFTarget io_target;
};

template < class DATA_TYPE, class IDX_TYPE = std::nullptr_t, typename ENABLE = void >
class HDFTabularIO;

// writing out scalar types
template < class DATA_TYPE >
class HDFTabularIO<DATA_TYPE, std::nullptr_t, typename std::enable_if< can_hdf_io< DATA_TYPE > >::type> : public HDFIO
{
public:
  HDFTabularIO(HDFTarget & target) : HDFIO(target), m_cell_size_map() {}
  virtual ~HDFTabularIO() {}

  virtual void CreateTable( string const & title,
                            string const & id,
                            localIndex const unit_per_cell,
                            localIndex const num_cols,
                            string const & record_prefix = "" )
  {
    impl::CreateTable<DATA_TYPE>(this->io_target,title,id,unit_per_cell,num_cols,record_prefix);
    m_cell_size_map[id] = unit_per_cell;
  }
  virtual void AppendRow( string const & id,
                          localIndex const num_cols,
                          DATA_TYPE const * data)
  {
    impl::AppendRow(this->io_target,id,m_cell_size_map[id],num_cols,data);
  }

  template <typename SEARCH_LAMBDA>
  // sfinae lamda callable, one input arg, bool return
  localIndex ColSearch( string const & id,
                        integer const col_idx,
                        SEARCH_LAMBDA && cond )
  {
    return impl::ColSearch<DATA_TYPE>(this->io_target,id,col_idx,m_cell_size_map[id],cond);
  }

  virtual void ClearAfter( string const & id,
                           localIndex const row_idx )
  {
    impl::ClearAfter(this->io_target,id,row_idx);
  }
private:
  map<string,localIndex> m_cell_size_map;
};

template < class DATA_TYPE, class IDX_TYPE >
class HDFTabularIO<DATA_TYPE, IDX_TYPE, typename std::enable_if< can_hdf_io< DATA_TYPE > && std::is_integral< IDX_TYPE >::value >::type> : public HDFIO
{
public:
  HDFTabularIO(HDFTarget & target) : HDFIO(target), m_cell_size_map() {}
  virtual ~HDFTabularIO() {}

  virtual void CreateTable( string const & title,
                            string const & id,
                            localIndex const unit_per_cell,
                            localIndex const num_cols,
                            IDX_TYPE const * idxs,
                            string const & record_prefix = "" )
  {
    impl::CreateTable<DATA_TYPE,IDX_TYPE>(this->io_target,title,id,unit_per_cell,num_cols,idxs,record_prefix);
    m_cell_size_map[id] = unit_per_cell;
  }

  virtual void AppendRow( string const & id,
                          DATA_TYPE const * data,
                          localIndex const num_cols,
                          IDX_TYPE const * idxs = nullptr )
  {
    impl::AppendRow(this->io_target,id,m_cell_size_map[id],num_cols,data,idxs);
  }

  template <typename SEARCH_LAMBDA>
  // sfinae lamda callable, one input arg, bool return
  localIndex ColSearch( string const & id,
                        integer const col_idx,
                        SEARCH_LAMBDA && cond )
  {
    return impl::ColSearch<DATA_TYPE>(this->io_target,id,col_idx,m_cell_size_map[id],cond);
  }

  virtual void ClearAfter( string const & id,
                           localIndex const row_idx )
  {
    impl::ClearAfter(this->io_target,id,row_idx);
  }

protected:
  map<string,localIndex> m_cell_size_map;
};

// if we want to take whole chunks of data, we must want all indices < the cell stride index

#define HDF_SLICE -2

template < typename ARRAY_TYPE, typename ENABLE = void >
class HDFTableIO;


template < typename T, int NDIM, typename PERMUTATION, typename INDEX_TYPE, template< typename > class DATA_VECTOR_TYPE >
class HDFTableIO< LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE >, typename std::enable_if< can_hdf_io< T > >::type > : public HDFIO
{
public:
  using array_type = LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE >;
  using index_type = array1d< INDEX_TYPE >;
  using component_type = array2d< INDEX_TYPE >;
  using value_type = typename array_type::value_type;

  HDFTableIO(HDFTarget & target,
              string const & title,
              string const & hdf_id,
              localIndex const csd,
              index_type const & idxs,
              component_type const & cmps,
              string const & record_prefix = "") :
    HDFIO(target),
    m_is_open(false),
    m_title(title),
    m_hdf_id(hdf_id),
    m_csd(csd),
    m_idxs(idxs),
    m_cmps(cmps),
    m_num_cells(0),
    m_cell_size(1),
    m_record_prefix(record_prefix),
    m_internal_copy()
  {
    // how many cells do we want to dump?
    m_num_cells = m_idxs.size( );

    // num indices speced for each cell
    m_cell_size = m_cmps.size( 0 );

    m_internal_copy.reserve(m_num_cells * m_cell_size);
  }

  virtual ~HDFTableIO() {}

  virtual void OpenTable( )
  {
    if ( !impl::TryOpen(this->io_target,m_hdf_id) )
    {
      impl::CreateTable<value_type>(this->io_target,m_title,m_hdf_id,m_cell_size,m_num_cells,m_record_prefix);
    }
    else
    {
      impl::VerifyTable(this->io_target,m_hdf_id,m_num_cells,m_cell_size);
    }
    m_is_open = true;
  }
  virtual void AppendRow( array_type const & row )
  {
    // assert(m_is_open);
    _update_internal_copy(row);
    impl::AppendRow<value_type>(this->io_target,m_hdf_id,m_cell_size,m_num_cells,m_internal_copy.data());
  }
  virtual void ClearAfter( localIndex first_to_delete )
  {
    //assert(m_is_open);
    impl::ClearAfter(this->io_target,m_hdf_id,first_to_delete);
  }
  template < typename LAMBDA >
  // sfinae lambda returns bool, accepts value_type*
  void ColSearch( localIndex col_idx , LAMBDA && cond)
  {
    return impl::ColSearch<value_type>(this->io_target,m_hdf_id,col_idx,m_cell_size,cond);
  }
  virtual void CloseTable( )
  {
    m_is_open = false;
  }

private:

  template < int U, class SLICE, class LAMBDA >
  inline
  typename std::enable_if< (U < NDIM-1), T const & >::type
  _index(SLICE & slice, LAMBDA indexer )
  {
    return _index<U+1>(slice[indexer(U)],indexer);
  }

  template < int U, class SLICE, class LAMBDA >
  inline
  typename std::enable_if< (U >= NDIM-1), T const & >::type
  _index(SLICE & slice, LAMBDA indexer )
  {
    return slice[indexer(U)];
  }

  void _update_internal_copy( array_type const & row )
  {
    localIndex offset = 0;
    for( hsize_t cell = 0; cell < m_num_cells; ++cell )
    {
      for ( localIndex cmp = 0; cmp < m_cmps.size( 0 ); ++cmp)
      {
        m_internal_copy[offset++] = _index<0>(row,[&](localIndex dim) { return (dim == m_csd) ? m_idxs[cell] : m_cmps[cmp][dim]; });
      }
    }
  }

  bool m_is_open;

  string const m_title;
  string const m_hdf_id;
  localIndex const m_csd;
  index_type m_idxs;
  component_type m_cmps;
  hsize_t m_num_cells;
  hsize_t m_cell_size;
  string const m_record_prefix;
  array1d< T > m_internal_copy;
};

template < class DATA_ARR_T >
class HDFTabularIO<DATA_ARR_T, std::nullptr_t, typename std::enable_if< is_array<DATA_ARR_T> && can_hdf_io<typename DATA_ARR_T::value_type> >::type > : public HDFIO
{
public:
  using value_type = typename DATA_ARR_T::value_type;

  HDFTabularIO(HDFTarget & target) : HDFIO(target), m_cell_size_map() {}
  virtual ~HDFTabularIO() {}

  virtual void CreateTable(string const & title,
                           string const & id,
                           localIndex const num_cols,
                           localIndex const unit_per_cell = 1,
                           string const & record_prefix = "" )
  {
    impl::CreateTable<value_type>(this->io_target,title,id,unit_per_cell,num_cols,record_prefix);
    m_cell_size_map[id] = unit_per_cell;
  }

  void AppendRow(string const & id,
                 DATA_ARR_T const & array)
  {
    localIndex num_cols = array.size();
    impl::AppendRow<value_type>(this->io_target,id,m_cell_size_map[id],num_cols,array.data());
  }

  void ClearAfter(string const & id,
                  localIndex const row_idx)
  {
    impl::ClearAfter(this->io_target,id,row_idx);
  }

private:
  map<string,localIndex> m_cell_size_map;
};

// currently assumes standard c array ordering
template < class DATA_ARR_T, class IDX_ARR_T >
//sfinae todo: idx_arr_t dim = 1 for single-level indexing, dim = 2 for sub-indexing
class HDFTabularIO<DATA_ARR_T, IDX_ARR_T, typename std::enable_if< is_array<DATA_ARR_T> &&
                                                                   can_hdf_io<typename DATA_ARR_T::value_type> &&
                                                                   is_array<IDX_ARR_T> &&
                                                                   std::is_integral< typename IDX_ARR_T::value_type >::value >::type > : public HDFIO
{
public:
  using value_type = typename DATA_ARR_T::value_type;
  using index_type = typename IDX_ARR_T::value_type;

  HDFTabularIO(HDFTarget & target) : HDFIO(target), m_cell_size_map() {}
  virtual ~HDFTabularIO() {}

  virtual void CreateTable( string const & title,
                            string const & id,
                            localIndex const unit_per_cell,
                            localIndex const num_cols,
                            IDX_ARR_T const & idxs,
                            string const & record_prefix = "" )
  {
    impl::CreateTable<value_type,index_type>(this->io_target,title,id,unit_per_cell,num_cols,idxs.data(),record_prefix);
    m_cell_size_map[id] = unit_per_cell;
  }

  void AppendRow(string const & id,
                 DATA_ARR_T const & array,
                 IDX_ARR_T const & idx)
  {
    localIndex num_cols = array.size();
    impl::AppendRow<value_type,index_type>(this->io_target,id,m_cell_size_map[id],num_cols,array.data(),idx.data());
  }

  void ClearAfter(string const & id,
                  localIndex const row_idx)
  {
    impl::ClearAfter(this->io_target,id,row_idx);
  }

private:
  map<string,localIndex> m_cell_size_map;
};


template < class DATA_ARR_T, class IDX_ARR_T = std::nullptr_t >
// the sfinae from the superclass should handle the sfinae here
// execpt this could potentially be either the pointer or the array version
// above, and the functions in those versions have different parameters, so we might need two
// specializations of this class toooo
class HDFTimeHistoryTabular : public HDFTabularIO<DATA_ARR_T,IDX_ARR_T>
{
public:
  using Super = HDFTabularIO<DATA_ARR_T,IDX_ARR_T>;
  using typename Super::value_type;

  HDFTimeHistoryTabular( HDFTarget & target ) :
    Super(target),
    time_table(target)
  { }

  virtual void CreateTable( string const & title,
                            string const & id,
                            localIndex const num_cols,
                            localIndex const type_per_cell = 1,
                            string const & record_prefix = "" )
  {
    Super::CreateTable(title,id,num_cols,type_per_cell,record_prefix);
    string ttitle = title + string(" Time");
    string tid = GenTimeID(id);
    time_table.CreateTable(ttitle,tid,1,1);
    m_data_2_time_map[id] = tid;
  }
  // virtual void LoadTableMeta();

  virtual void AppendRow(string const & id,// localIndex const num_cols,
                         real64 const time,
                         DATA_ARR_T const & array)//,
                         //localIndex const type_per_cell = 1)
  {
    Super::AppendRow(id,array);
    time_table.AppendRow(m_data_2_time_map[id],1,&time);
  }

  bool ClearAfterTime(string const & id, real64 const time)
  {
    localIndex first_to_clear = time_table.ColSearch(m_data_2_time_map[id],0,[&](real64 * col_time) { return time >= *col_time; });
    bool do_clear = first_to_clear >= 0;
    if ( do_clear )
    {
      Super::ClearAfter(id,first_to_clear);
      time_table.ClearAfter(m_data_2_time_map[id],first_to_clear);
    }
    return do_clear;
  }

private:
  string GenTimeID(string const & data_id) { return data_id + string("_t"); }
  map<string,string> m_data_2_time_map;
  HDFTabularIO<real64> time_table;
};

// (X) table that contains a whole array in a row
// ( ) table that contains specific indices of an array in a row
// ( ) table that contains specific subindices (e.g. the nth component of every entry) of an array in a row
// ( ) table that contains specific subindices of specific indices of an array in a row



}