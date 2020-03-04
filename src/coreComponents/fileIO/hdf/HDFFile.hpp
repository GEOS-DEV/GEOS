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

template <>
hid_t GetHDFDataType<R1Tensor>()
{
  return GetHDFArrayDataType<real64>(3);
}

template <>
hid_t GetHDFDataType<R2Tensor>()
{
  hsize_t dims[2] = {3,3};
  return GetHDFArrayDataType<real64>(2,&dims[0]);
}

template <>
hid_t GetHDFDataType<R2SymTensor>()
{
  return GetHDFArrayDataType<real64>(6);
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

  void ClearFrom( hid_t target,
                   string const & id,
                   localIndex const row_idx )
  {
    hsize_t num_cols = 0;
    hsize_t num_rows = 0;
    H5TBget_table_info(target,id.c_str(),&num_cols,&num_rows);
    H5TBdelete_record(target,id.c_str(),row_idx,num_rows-row_idx);
  }

/// herelper classes for iteratively indexing into lvarrays at runtime
  // template < class INDEX_TYPE, int NDIM >
  // class Indexer;

  template <class INDEX_TYPE, int NDIM >
  class DenseIndexer
  {
  public:
    using index_array_type = array2d<INDEX_TYPE>;

    template < typename ARRAY_TYPE, typename = typename std::enable_if< is_array<ARRAY_TYPE> >::type >
    DenseIndexer( ARRAY_TYPE const & arr ) :
      m_cmp(0),
      m_bounds(0),
      m_num_indices(0),
      m_dim_index(0)
    {
      for( localIndex dd = 0; dd < NDIM; ++dd)
      {
        m_bounds[dd] = arr.size( dd );
      }
      m_num_indices = std::accumulate(m_bounds,m_bounds+NDIM,1,std::multiplies<INDEX_TYPE>());
    }

    template < typename ... BOUNDS >
    DenseIndexer( BOUNDS const ... bounds ) :
      m_cmp(0),
      m_bounds({0}),
      m_num_indices(0),
      m_dim_index({0})
    {
      static_assert( sizeof ... (BOUNDS) == NDIM, "Error: calling DenseIndexer::DenseIndexer with incorrect number of arguments." );
      LvArray::dimUnpack(&m_bounds[0],bounds...);
      m_num_indices = std::accumulate(m_bounds.begin(),m_bounds.end(),1,std::multiplies<INDEX_TYPE>());
    }
    void next()
    {
      this->operator++();
    }
    bool done()
    {
      return m_cmp != getNumIndices();
    }
    void reset()
    {
      m_cmp = 0;
      memset(&m_dim_index[0],0,sizeof(localIndex)*NDIM);
    }
    inline INDEX_TYPE operator()(localIndex dim)
    {
      return m_dim_index[dim];
    }
    inline localIndex operator++()
    {
      localIndex dim = NDIM - 1;
      m_dim_index[dim]++;
      while(m_dim_index[dim] >= m_bounds[0])
      {
        m_dim_index[dim] = 0;
        dim--;
        m_dim_index[dim]++;
      }
      return dim;
    }
    INDEX_TYPE getNumIndices()
    {
      return m_num_indices;
    }
  private:
    localIndex m_cmp;
    std::array<INDEX_TYPE, NDIM> m_bounds;
    localIndex m_num_indices;
    std::array<localIndex, NDIM> m_dim_index;
  };

  template <class INDEX_TYPE >
  class DenseIndexer<INDEX_TYPE,0>
  {
  public:
    template < typename ARRAY_TYPE, typename = typename std::enable_if< is_array<ARRAY_TYPE> >::type >
    DenseIndexer( ARRAY_TYPE const & ) : m_done(false) {}
    template < typename ... BOUNDS >
    DenseIndexer( BOUNDS const ... ) : m_done(false) {}
    void next() { this->operator++(); };
    bool done() { return m_done; }
    void reset() { m_done = false; }
    // should never be called
    inline INDEX_TYPE operator()(localIndex) { return -1; }
    inline localIndex operator++() { m_done = true; return 0; }
    INDEX_TYPE getNumIndices() { return 0; }
  private:
    bool m_done;
  };

  template < class INDEX_TYPE, int NDIM >
  class EnumeratedIndexer
  {
  public:
    using index_array_type = array2d<INDEX_TYPE>;
    EnumeratedIndexer(index_array_type const & idxs) :
      m_idxs(idxs),
      m_cmp(0)
    { }
    void next()
    {
      this->operator++();
    }
    bool done()
    {
      return m_cmp != getNumIndices();
    }
    inline INDEX_TYPE operator()(localIndex dim)
    {
      return m_idxs[m_cmp][dim];
    }
    inline localIndex operator++()
    {
      return m_cmp++;
    }
    void reset()
    {
      m_cmp = 0;
    }
    INDEX_TYPE getNumIndices() { return m_idxs.size( 0 ); }
  private:
    array2d<INDEX_TYPE> const & m_idxs;
    localIndex m_cmp;
  };

  template < class INDEX_TYPE >
  class EnumeratedIndexer<INDEX_TYPE,1>
  {
  public:
    using index_array_type = array1d<INDEX_TYPE>;
    EnumeratedIndexer(index_array_type const & idxs) :
      m_idxs(idxs),
      m_cmp(0)
    {}
    void next()
    {
      this->operator++();
    }
    bool done()
    {
      return m_cmp != getNumIndices();
    }
    inline INDEX_TYPE operator()(localIndex)
    {
      return m_idxs[m_cmp];
    }
    inline localIndex operator++()
    {
      return m_cmp++;
    }
    void reset()
    {
      m_cmp = 0;
    }
    INDEX_TYPE getNumIndices() { return m_idxs.size( ); }
  private:
    array1d<INDEX_TYPE> const & m_idxs;
    localIndex m_cmp;
  };

template< int U, typename T, int NDIM, int USD, typename INDEX_TYPE, typename LAMBDA >
  inline
  typename std::enable_if< ( U >= 0 ) && (U < NDIM-1), T const & >::type
  _index(LvArray::ArraySlice<T,NDIM,USD,INDEX_TYPE> const & slice, LAMBDA indexer )
  {
    return _index<U+1>(slice[indexer(U)],indexer);
  }

template< int U, typename T, int NDIM, int USD, typename INDEX_TYPE, typename LAMBDA >
  inline
  typename std::enable_if< ( U >= 0 ) && (U >= NDIM-1), T const & >::type
  _index(LvArray::ArraySlice<T,NDIM,USD,INDEX_TYPE> const & slice, LAMBDA indexer )
  {
    return slice[indexer(U)];
  }

template< int U, typename T, int NDIM, typename PERMUTATION, typename INDEX_TYPE, template< typename > class DATA_VECTOR_TYPE, typename LAMBDA >
  inline
  typename std::enable_if< ( U >= 0 ) && (U < NDIM-1), T const & >::type
  _index(LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE > const & slice, LAMBDA indexer )
  {
    return _index<U+1>(slice[indexer(U)],indexer);
  }

template< int U, typename T, int NDIM, typename PERMUTATION, typename INDEX_TYPE, template< typename > class DATA_VECTOR_TYPE, typename LAMBDA >
  inline
  typename std::enable_if< ( U >= 0 ) && (U >= NDIM-1), T const & >::type
  _index(LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE > const & slice, LAMBDA indexer )
  {
    return slice[indexer(U)];
  }

template < typename T, int NDIM, typename PERMUTATION, typename INDEX_TYPE, template< typename > class DATA_VECTOR_TYPE , class FLAT_TYPE, class CELL_INDEXER, class COMP_INDEXER >
  inline void _process_flatten( LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE > const & arr,
                                FLAT_TYPE & flat_arr,
                                CELL_INDEXER & cell_indexer,
                                COMP_INDEXER & comp_indexer,
                                localIndex csd = 0)
  {
    localIndex offset = 0;
    for( cell_indexer.reset(); !cell_indexer.done(); ++cell_indexer )
    {
      for(comp_indexer.reset(); !comp_indexer.done(); ++comp_indexer )
      {
        flat_arr[offset++] = _index<0>(arr,[&](localIndex dim) { return (dim == csd) ? cell_indexer(0) : comp_indexer(dim > csd ? dim-1 : dim);});
      }
    }
  }

  template < typename T, int NDIM, typename PERMUTATION, typename INDEX_TYPE, template< typename > class DATA_VECTOR_TYPE >
  void FlattenArray( LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE > const & arr,
                     array1d< T > & flat_arr,
                     array1d<INDEX_TYPE> & idxs,
                     array2d<INDEX_TYPE> & cmps,
                     localIndex const csd = 0)
  {
    if ( idxs.size( ) > 0 && cmps.size( ) == 0 )
    {
      EnumeratedIndexer<INDEX_TYPE,1> cell_indexer(idxs);
      DenseIndexer<INDEX_TYPE,NDIM-1> cmp_indexer(arr);
      _process_flatten(arr,flat_arr,cell_indexer,cmp_indexer,csd);
    }
    else if( idxs.size( ) > 0 && cmps.size( ) > 0 )
    {
      EnumeratedIndexer<INDEX_TYPE,1> cell_indexer(idxs);
      EnumeratedIndexer<INDEX_TYPE,NDIM-1> cmp_indexer(cmps);
      _process_flatten(arr,flat_arr,cell_indexer,cmp_indexer,csd);
    }
    else if( idxs.size( ) == 0 && cmps.size( ) > 0 )
    {
      DenseIndexer<INDEX_TYPE,1> cell_indexer(arr.size(csd));
      EnumeratedIndexer<INDEX_TYPE,NDIM-1> cmp_indexer(cmps);
      _process_flatten(arr,flat_arr,cell_indexer,cmp_indexer,csd);
    }
    else if( idxs.size( ) == 0 && cmps.size( ) == 0)
    {
      DenseIndexer<INDEX_TYPE,1> cell_indexer(arr.size(csd));
      DenseIndexer<INDEX_TYPE,NDIM-1> cmp_indexer(cmps);
      _process_flatten(arr,flat_arr,cell_indexer,cmp_indexer,csd);
    }
  }

}// impl namespace

template < typename DATA_TYPE, typename ENABLE = void >
class HDFTableIO;


// this can only handle dense contiguous output at the moment, which is fine for time series
template < typename DATA_TYPE >
class HDFTableIO< DATA_TYPE, typename std::enable_if< can_hdf_io< DATA_TYPE > >::type >
{
public:
  using value_type = DATA_TYPE;

  HDFTableIO( string const & title,
              string const & hdf_id,
              localIndex const & num_cells,
              localIndex const & cell_size,
              string const & record_prefix  = "" ) :
    m_is_open(false),
    m_title(title),
    m_hdf_id(hdf_id),
    m_num_cells(num_cells),
    m_cell_size(cell_size),
    m_record_prefix(record_prefix)
  { }

  virtual void OpenTable( HDFTarget & target )
  {
    if ( !impl::TryOpen(target,m_hdf_id) )
    {
      impl::CreateTable<value_type>(target,m_title,m_hdf_id,m_cell_size,m_num_cells,m_record_prefix);
    }
    else
    {
      impl::VerifyTable(target,m_hdf_id,m_num_cells,m_cell_size);
    }
    m_is_open = true;
    m_active_target = target;
  }
  virtual void AppendRow( value_type const * row )
  {
    // assert(m_is_open);
    /// only supporting flat arrays at the moment
    impl::AppendRow<value_type>(m_active_target,m_hdf_id,m_cell_size,m_num_cells,row);
  }
  virtual void ClearFrom( localIndex first_to_delete )
  {
    //assert(m_is_open);
    impl::ClearFrom(m_active_target,m_hdf_id,first_to_delete);
  }
  template < typename LAMBDA >
  // sfinae lambda returns bool, accepts value_type*
  localIndex ColSearch( localIndex col_idx , LAMBDA && cond)
  {
    //assert(m_is_open);
    return impl::ColSearch<value_type>(m_active_target,m_hdf_id,col_idx,m_cell_size,cond);
  }
  virtual void CloseTable( )
  {
    m_is_open = false;
  }

private:

  bool m_is_open;
  string const m_title;
  string const m_hdf_id;
  hsize_t m_num_cells;
  hsize_t m_cell_size;
  string const m_record_prefix;
  HDFTarget m_active_target;
};

// specialization for LvArrays containing types we're able to output with hdf
template < typename T, int NDIM, typename PERMUTATION, typename INDEX_TYPE, template< typename > class DATA_VECTOR_TYPE >
class HDFTableIO< LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE >, typename std::enable_if< can_hdf_io< T > >::type >
{
public:
  using array_type = LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE >;
  using index_array_type = array1d< INDEX_TYPE >;
  using comp_array_type = array2d< INDEX_TYPE >;
  using flat_array_type = array1d< T >;
  using value_type = typename array_type::value_type;

  HDFTableIO( string const & title,
              string const & hdf_id,
              INDEX_TYPE const & num_cells,
              localIndex const & cell_size,
              localIndex const csd = 0,
              string const & record_prefix  = "" ) :
  m_is_open(false),
  m_title(title),
  m_hdf_id(hdf_id),
  m_csd(csd),
  m_idxs(),
  m_cmps(),
  m_num_cells(num_cells),
  m_cell_size(cell_size),
  m_record_prefix(record_prefix),
  m_internal_copy()
  {
    _reserve();
  }

  HDFTableIO( string const & title,
              string const & hdf_id,
              index_array_type const & idxs,
              localIndex const & cell_size,
              localIndex const csd = 0 ,
              string const & record_prefix = "" ) :
  m_is_open(false),
  m_title(title),
  m_hdf_id(hdf_id),
  m_csd(csd),
  m_idxs(idxs),
  m_cmps(),
  m_num_cells(0),
  m_cell_size(cell_size),
  m_record_prefix(record_prefix),
  m_internal_copy()
  {
    m_num_cells = m_idxs.size( );

    _reserve();
  }

  HDFTableIO( string const & title,
              string const & hdf_id,
              INDEX_TYPE const & num_cells,
              comp_array_type const & cmps,
              localIndex const csd = 0,
              string const & record_prefix = "" ) :
    m_is_open(false),
    m_title(title),
    m_hdf_id(hdf_id),
    m_csd(csd),
    m_idxs(),
    m_cmps(cmps),
    m_num_cells(num_cells),
    m_cell_size(1),
    m_record_prefix(record_prefix),
    m_internal_copy()
  {
    // num indices speced for each cell
    m_cell_size = m_cmps.size( 0 );

    _reserve();
  }

  HDFTableIO( string const & title,
              string const & hdf_id,
              localIndex const csd,
              index_array_type const & idxs,
              comp_array_type const & cmps,
              string const & record_prefix = "" ) :
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

    _reserve();
  }

  virtual ~HDFTableIO() {}

  virtual void OpenTable( HDFTarget & target )
  {
    if ( !impl::TryOpen(target,m_hdf_id) )
    {
      impl::CreateTable<value_type>(target,m_title,m_hdf_id,m_cell_size,m_num_cells,m_record_prefix);
    }
    else
    {
      impl::VerifyTable(target,m_hdf_id,m_num_cells,m_cell_size);
    }
    m_is_open = true;
    m_active_target = target;
  }
  virtual void AppendRow( array_type const & row )
  {
    // assert(m_is_open);
    impl::FlattenArray(row,m_internal_copy,m_idxs,m_cmps,m_csd);
    impl::AppendRow<value_type>(m_active_target,m_hdf_id,m_cell_size,m_num_cells,m_internal_copy.data());
  }
  virtual void ClearFrom( localIndex first_to_delete )
  {
    //assert(m_is_open);
    impl::ClearFrom(m_active_target,m_hdf_id,first_to_delete);
  }
  template < typename LAMBDA >
  // sfinae lambda returns bool, accepts value_type*
  localIndex ColSearch( localIndex col_idx , LAMBDA && cond)
  {
    //assert(m_is_open);
    return impl::ColSearch<value_type>(m_active_target,m_hdf_id,col_idx,m_cell_size,cond);
  }
  virtual void CloseTable( )
  {
    m_is_open = false;
  }

protected:

  void _reserve()
  {
    m_internal_copy.reserve(m_num_cells * m_cell_size);
  }

  void _release()
  {
    m_internal_copy.swap(flat_array_type());
  }

  bool m_is_open;
  string const m_title;
  string const m_hdf_id;
  localIndex const m_csd;
  index_array_type m_idxs;
  comp_array_type m_cmps;
  hsize_t m_num_cells;
  hsize_t m_cell_size;
  string const m_record_prefix;
  flat_array_type m_internal_copy;
  HDFTarget m_active_target;
};

template < typename DATA_TYPE >
class HDFTimeHistoryTable;

template < typename T, int NDIM, typename PERMUTATION, typename INDEX_TYPE, template< typename > class DATA_VECTOR_TYPE >
class HDFTimeHistoryTable< LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE > > : public HDFTableIO< LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE >, typename std::enable_if< can_hdf_io< T > >::type >
{
public:
  using Super = HDFTableIO< LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, DATA_VECTOR_TYPE >, typename std::enable_if< can_hdf_io< T > >::type >;
  using typename Super::array_type;
  using typename Super::index_array_type;
  using typename Super::comp_array_type;
  using typename Super::flat_array_type;
  using typename Super::value_type;

  HDFTimeHistoryTable( string const & title,
              string const & hdf_id,
              INDEX_TYPE const & num_cells,
              localIndex const & cell_size,
              localIndex const csd = 0,
              string const & record_prefix  = "" ) :
    Super(title,hdf_id,num_cells,cell_size,csd,record_prefix),
    m_time_table(GenTimeTitle(title),GenTimeID(hdf_id),1,1)
  {}

  HDFTimeHistoryTable( string const & title,
              string const & hdf_id,
              index_array_type const & idxs,
              localIndex const & cell_size,
              localIndex const csd = 0 ,
              string const & record_prefix = "" ) :
    Super(title,hdf_id,idxs,cell_size,csd,record_prefix),
    m_time_table(GenTimeTitle(title),GenTimeID(hdf_id),1,1)
  {}

  HDFTimeHistoryTable( string const & title,
              string const & hdf_id,
              INDEX_TYPE const & num_cells,
              comp_array_type const & cmps,
              localIndex const csd = 0,
              string const & record_prefix = "" ) :
    Super(title,hdf_id,num_cells,cmps,csd,record_prefix),
    m_time_table(GenTimeTitle(title),GenTimeID(hdf_id),1,1)
  {}

  HDFTimeHistoryTable( string const & title,
              string const & hdf_id,
              localIndex const csd,
              index_array_type const & idxs,
              comp_array_type const & cmps,
              string const & record_prefix = "" ) :
    Super(title,hdf_id,idxs,cmps,csd,record_prefix),
    m_time_table(GenTimeTitle(title),GenTimeID(hdf_id),1,1)
  {}

  virtual void OpenTable( HDFTarget & target )
  {
    Super::OpenTable(target);
    m_time_table.OpenTable(target);
  }
  virtual void AppendRow( real64 time, array_type const & row )
  {
    Super::AppendRow(row);
    m_time_table.AppendRow(&time);
  }

  virtual void ClearFromRow( localIndex first_to_clear )
  {
    //assert(m_is_open);
    Super::ClearFrom(first_to_clear);
    m_time_table.ClearFrom(first_to_clear);
  }

  localIndex FirstRowTimeGE( real64 const time )
  {
    return m_time_table.ColSearch(0,[&](real64 * col_time) { return time >= *col_time; });
  }

  bool ClearFromTime( real64 const time )
  {
    localIndex first = FirstRowTimeGE(time);
    bool did_clear = (first >= 0);
    if ( did_clear ) ClearFromRow(first);
    return did_clear;
  }

  virtual void CloseTable( )
  {
    Super::CloseTable( );
    m_time_table.CloseTable( );
  }

private:
  //disallow appending to the data table without also appending to the time table
  using Super::AppendRow;
  static string GenTimeID(string const & id) { return id + string("_t"); }
  static string GenTimeTitle(string const & title) { return title + string(" Time"); }
  HDFTableIO<real64> m_time_table;
};




// template < class DATA_ARR_T, class IDX_ARR_T = std::nullptr_t >
// // the sfinae from the superclass should handle the sfinae here
// // execpt this could potentially be either the pointer or the array version
// // above, and the functions in those versions have different parameters, so we might need two
// // specializations of this class toooo
// class HDFTimeHistoryTabular : public HDFTabularIO<DATA_ARR_T,IDX_ARR_T>
// {
// public:
//   using Super = HDFTabularIO<DATA_ARR_T,IDX_ARR_T>;
//   using typename Super::value_type;

//   HDFTimeHistoryTabular( HDFTarget & target ) :
//     Super(target),
//     time_table(target)
//   { }

//   virtual void CreateTable( string const & title,
//                             string const & id,
//                             localIndex const num_cols,
//                             localIndex const type_per_cell = 1,
//                             string const & record_prefix = "" )
//   {
//     Super::CreateTable(title,id,num_cols,type_per_cell,record_prefix);
//     string ttitle = title + string(" Time");
//     string tid = GenTimeID(id);
//     time_table.CreateTable(ttitle,tid,1,1);
//     m_data_2_time_map[id] = tid;
//   }
//   // virtual void LoadTableMeta();

//   virtual void AppendRow(string const & id,// localIndex const num_cols,
//                          real64 const time,
//                          DATA_ARR_T const & array)//,
//                          //localIndex const type_per_cell = 1)
//   {
//     Super::AppendRow(id,array);
//     time_table.AppendRow(m_data_2_time_map[id],1,&time);
//   }

//   bool ClearFromTime(string const & id, real64 const time)
//   {
//     localIndex first_to_clear = time_table.ColSearch(m_data_2_time_map[id],0,[&](real64 * col_time) { return time >= *col_time; });
//     bool do_clear = first_to_clear >= 0;
//     if ( do_clear )
//     {
//       Super::ClearFrom(id,first_to_clear);
//       time_table.ClearFrom(m_data_2_time_map[id],first_to_clear);
//     }
//     return do_clear;
//   }

// private:
//   string GenTimeID(string const & data_id) { return data_id + string("_t"); }
//   map<string,string> m_data_2_time_map;
//   HDFTabularIO<real64> time_table;
// };

// (X) table that contains a whole array in a row
// ( ) table that contains specific indices of an array in a row
// ( ) table that contains specific subindices (e.g. the nth component of every entry) of an array in a row
// ( ) table that contains specific subindices of specific indices of an array in a row



}