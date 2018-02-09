/**
 * @file ViewWrapper.hpp
 *
 * @date Created on: Jun 8, 2016
 * @authors settgast
 */

#ifndef GEOSX_DATAREPOSITORY_WRAPPERVIEW_HPP_
#define GEOSX_DATAREPOSITORY_WRAPPERVIEW_HPP_

#include "ViewWrapperBase.hpp"
#include "KeyNames.hpp"
#include "common/DataTypes.hpp"
#include "SFINAE_Macros.hpp"
#include <type_traits>
#include "StringUtilities.hpp"
#include "Macros.hpp"

#include "MPI_Communications/CommBufferOps.hpp"


#ifdef USE_ATK
#include "sidre/sidre.hpp"
#include "sidre/SidreTypes.hpp"
#endif


namespace geosx
{
namespace dataRepository
{

/**
 * Templated class to serve as a wrapper to arbitrary objects.
 * @tparam T is any object that is to be wrapped by ViewWrapper
 */
template< typename T >
class ViewWrapper : public ViewWrapperBase
{

public:
  /**
   * @param name name of the object
   * @param parent parent group which owns the ViewWrapper
   */
  explicit ViewWrapper( std::string const & name,
                        ManagedGroup * const parent ):
    ViewWrapperBase(name,parent),
    m_data( std::make_unique<T>() )
  {}

  /**
   * @param name name of the object
   * @param parent parent group that owns the ViewWrapper
   * @param object object that is being wrapped by the ViewWrapper
   */
  explicit ViewWrapper( std::string const & name,
                        ManagedGroup * const parent,
                        std::unique_ptr<T> object ):
    ViewWrapperBase(name,parent),
    m_data( std::move( object ) )
  {}

  /**
   * @param name name of the object
   * @param parent parent group that owns the ViewWrapper
   * @param object object that is being wrapped by the ViewWrapper
   */
  explicit ViewWrapper( std::string const & name,
                        ManagedGroup * const parent,
                        T * object ):
    ViewWrapperBase(name,parent),
    m_data( std::move( std::unique_ptr<T>(object) ) )
  {}

  /**
   * default destructor
   */
  virtual ~ViewWrapper() noexcept override final {}

  /**
   * Copy Constructor
   * @param source source for the copy
   */
  ViewWrapper( ViewWrapper const & source ):
    ViewWrapperBase("test", nullptr),
    m_data(source.m_data)
  {}

  /**
   * Move Constructor
   * @param source source to be moved
   */
  ViewWrapper( ViewWrapper&& source ):
    ViewWrapperBase("test", nullptr),
    m_data( std::move(source.m_data) )
  {}

  /**
   * Copy Assignment Operator
   * @param source rhs
   * @return *this
   */
  ViewWrapper& operator=( ViewWrapper const & source )
  {
    m_data = source.m_data;
    return *this;
  }

  /**
   * Move Assignment Operator
   * @param source
   * @return *this
   */
  ViewWrapper& operator=( ViewWrapper && source )
  {
    m_data = std::move(source.m_data);
    return *this;
  }


  /**
   * Factory Method to make a new ViewWrapper<T>, allocating a new T. Only is
   * going to work if T has a default constructor.
   * Perhaps this is worthless in the general case.
   * @param name name of the object
   * @param parent group that owns the ViewWrapper
   * @return A std::unique_ptr<ViewWrapperBase> that holds the newly allocated
   * ViewWrapper.
   */
  template<typename TNEW>
  static std::unique_ptr<ViewWrapperBase> Factory( std::string const & name,
                                                   ManagedGroup * const parent )
  {
    std::unique_ptr<TNEW> newObject = std::move( std::make_unique<TNEW>() );
    return std::move(std::make_unique<ViewWrapper<T> >( name, parent, std::move(newObject) ) );
  }


  /**
   * Virtual function to return the typeid of T. Not so sure this does what we
   * want?? TODO
   * @return typeid(T)
   */
  virtual const std::type_info& get_typeid() const noexcept override final
  {
    return typeid(T);
  }

  /**
   * static function to cast a ViewWrapper base to a derived ViewWrapper<T>
   * @param base
   * @return casted ViewWrapper<T>
   */
  static ViewWrapper<T>& cast( ViewWrapperBase& base )
  {
    if( base.get_typeid() != typeid(T) )
    {
#ifdef USE_ATK
      SLIC_ERROR("invalid cast attempt");
#endif
    }
    return static_cast< ViewWrapper<T>& >(base);
  }

  template< bool DO_PACKING >
  struct pack_wrapper
  {
    HAS_ALIAS(pointer)


    template<class U = T>
    static typename std::enable_if< has_alias_pointer<U>::value  ||
                                    std::is_pod<U>::value, int >::type
    PackT( char *& buffer, U const & m_data )
    {
      return CommBufferOps::Pack<DO_PACKING>( buffer, m_data );
    }

    template<class U = T>
    static typename std::enable_if< !has_alias_pointer<U>::value &&
                                    !std::is_pod<U>::value, int >::type
    PackT( char *& buffer, U const & m_data )
    {
      return 0;
    }




    template<class U = T>
    static typename std::enable_if< has_alias_pointer<U>::value  &&
                                    !std::is_pod<U>::value &&
                                    !std::is_same<U,string>::value, int >::type
    PackT( char *& buffer, U const & m_data, localIndex_array const & packList )
    {
      return CommBufferOps::Pack<DO_PACKING>( buffer, m_data, packList );
    }

    template<class U = T>
    static typename std::enable_if< ( std::is_pod<U>::value ||
                                      std::is_same<U,string>::value ),int >::type
    PackT( char *& buffer, U const & m_data, localIndex_array const & packList )
    {
      return CommBufferOps::Pack<DO_PACKING>( buffer, m_data );
    }

    template<class U = T>
    static typename std::enable_if< !has_alias_pointer<U>::value &&
                                    !std::is_pod<U>::value &&
                                    !std::is_same<U,string>::value, int >::type
    PackT( char *& buffer, U const & m_data, localIndex_array const & packList )
    {
      return 0;
    }


  };

  virtual int Pack( char *& buffer ) const override final
  {
    int packedSize = 0;

    packedSize += CommBufferOps::Pack<true>( buffer, this->getName() );
    packedSize += pack_wrapper<true>::PackT( buffer, *m_data);

    return packedSize;
  }

  virtual int Pack( char *& buffer, localIndex_array const & packList ) const override final
  {
    int packedSize = 0;

    packedSize += CommBufferOps::Pack<true>( buffer, this->getName() );
    packedSize += pack_wrapper<true>::PackT( buffer, *m_data, packList);

    return packedSize;
  }

  virtual int PackSize( ) const override final
  {
    char * buffer = nullptr;
    int packedSize = 0;

    packedSize += CommBufferOps::Pack<false>( buffer, this->getName() );
    packedSize += pack_wrapper<false>::PackT( buffer, *m_data);

    return packedSize;
  }

  virtual int PackSize( localIndex_array const & packList ) const override final
  {

    char * buffer = nullptr;
    int packedSize = 0;

    packedSize += CommBufferOps::Pack<false>( buffer, this->getName() );
    packedSize += pack_wrapper<false>::PackT( buffer, *m_data, packList);

    return packedSize;
  }




  struct empty_wrapper
  {
    HAS_MEMBER_FUNCTION(empty,bool,const,,)
    template<class U = T>
    static typename std::enable_if<has_memberfunction_empty<U>::value, bool>::type
    empty(ViewWrapper const * parent)
    {
      return parent->m_data->empty();
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_empty<U>::value, bool>::type
    empty(ViewWrapper const * parent)
    {
      return parent;
    }
  };
  virtual bool empty() const override final
  {
    return empty_wrapper::empty(this);
  }

  struct size_wrapper
  {
    HAS_MEMBER_FUNCTION_VARIANT(size,0,localIndex,const,,)
    template<class U = T>
    static typename std::enable_if< has_memberfunction_v0_size<U>::value,localIndex >::type
    size(ViewWrapper const * parent)
    {
      return static_cast<localIndex>(parent->m_data->size());
    }
    template<class U = T>
    static typename std::enable_if< !(has_memberfunction_v0_size<U>::value), localIndex>::type
    size(ViewWrapper const * )
    {
      return 1;//parent->m_data;
    }
  };
  virtual localIndex size() const override final
  {
    return size_wrapper::size(this);
  }


  struct reserve_wrapper
  {
    HAS_MEMBER_FUNCTION(reserve, void, ,VA_LIST(std::size_t),VA_LIST(std::size_t(1)) )
    template<class U = T>
    static typename std::enable_if<has_memberfunction_reserve<U>::value, void>::type
    reserve(ViewWrapper * const parent, std::size_t new_cap)
    {
      return parent->m_data->reserve(new_cap);
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_reserve<U>::value, void>::type
    reserve(ViewWrapper * const, std::size_t )
    {
      return; //parent->m_data;
    }
  };
  virtual void reserve( std::size_t new_cap ) override final
  {
    reserve_wrapper::reserve(this, new_cap);
  }
//  CONDITIONAL_VIRTUAL_FUNCTION( Wrapper<T>,reserve , void,,
// VA_LIST(std::size_t a), VA_LIST(a) )


  HAS_MEMBER_FUNCTION(capacity,std::size_t,const,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(ViewWrapper<T>,capacity,std::size_t,const)

  HAS_MEMBER_FUNCTION(max_size,std::size_t,const,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(ViewWrapper<T>,max_size,std::size_t,const)

  HAS_MEMBER_FUNCTION(clear,void,,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(ViewWrapper<T>,clear,void,)

  HAS_MEMBER_FUNCTION(insert,void,,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(ViewWrapper<T>,insert,void,)


  struct resize_wrapper
  {
    HAS_MEMBER_FUNCTION_VARIANT(resize,0,void,,VA_LIST(int), VA_LIST(static_cast<int>(1)))
    HAS_MEMBER_FUNCTION_VARIANT(resize,1,void,,VA_LIST(unsigned int), VA_LIST( static_cast<unsigned int>(1)))
    HAS_MEMBER_FUNCTION_VARIANT(resize,2,void,,VA_LIST(long), VA_LIST( static_cast<long int>(1)))
    HAS_MEMBER_FUNCTION_VARIANT(resize,3,void,,VA_LIST(unsigned long), VA_LIST(static_cast<unsigned long int>(1)))
    HAS_MEMBER_FUNCTION_VARIANT(resize,4,void,,VA_LIST(long long int), VA_LIST(static_cast<long long int>(1)))
    HAS_MEMBER_FUNCTION_VARIANT(resize,5,void,,VA_LIST(unsigned long long), VA_LIST(static_cast<unsigned long long>(1)))


    template<class U = T>
    static typename std::enable_if< has_memberfunction_v0_resize<U>::value ||
                                    has_memberfunction_v1_resize<U>::value ||
                                    has_memberfunction_v2_resize<U>::value ||
                                    has_memberfunction_v3_resize<U>::value ||
                                    has_memberfunction_v4_resize<U>::value, void>::type
    resize(ViewWrapper * const parent, localIndex const new_size)
    {
      return parent->m_data->resize(new_size);
    }


    template<class U = T>
    static typename std::enable_if<!(has_memberfunction_v0_resize<U>::value)&&
                                   !(has_memberfunction_v1_resize<U>::value)&&
                                   !(has_memberfunction_v2_resize<U>::value)&&
                                   !(has_memberfunction_v3_resize<U>::value)&&
                                   !(has_memberfunction_v4_resize<U>::value), void>::type
    resize(ViewWrapper * const, localIndex )
    {
      return;
    }
  };
  using ViewWrapperBase::resize;
  void resize( localIndex new_size ) override final
  {
    resize_wrapper::resize(this, new_size);
  }

  /**
   * @name Structure to determine return types for data access functions
   */
  ///@{

  /// Invoke macro to generate test to see if type has an alias named "pointer".
  /// This will be used to determine if the
  /// type is to be treated as an "array" or a single object.
  HAS_ALIAS(pointer)

  /**
   * SFINAE specialized structure to control return type based on properties of
   * T.
   * The default template returns a pointer for all calls to data().
   */
  template< class U=T,
            bool HASPOINTERTYPE = has_alias_pointer<U>::value,
            bool ISSTRING = std::is_same<U,std::string>::value >
  struct Get_Type
  {
    typedef U*       type;
    typedef U const * const_type;

    typedef U *       pointer;
    typedef U const * const_pointer;
  };

  /**
   *  Specialization for case when T has a pointer alias, and it is NOT a
   * string.
   *  In this case, we assume that we are storing an array type. The return type
   * is then a reference, unless the
   *  compilation flag is set such that we require a pointer back (good for
   * speed, but no array class convenience).
   *  The resulting types can both be dereferenced with operator[], so no code
   * changes required
   *  unless array member functions have been called.
   */
  template<class U>
  struct Get_Type<U, true, false>
  {

#if CONTAINERARRAY_RETURN_PTR == 1
    typedef typename U::pointer       type;
    typedef typename U::const_pointer const_type;
#else
    typedef U &       type;
    typedef U const & const_type;
#endif
    typedef typename U::pointer       pointer;
    typedef typename U::const_pointer const_pointer;
  };


  /// Specialization for string. Always return a reference.
  template<class U>
  struct Get_Type<U, true, true>
  {
    typedef U &       type;
    typedef U const & const_type;

    typedef U *       pointer;
    typedef U const * const_pointer;
  };

  using rtype       = typename Get_Type<T>::type;
  using rtype_const = typename Get_Type<T>::const_type;

  using pointer       = typename Get_Type<T>::pointer;
  using const_pointer = typename Get_Type<T>::const_pointer;
  ///@}


  HAS_MEMBER_FUNCTION(data,pointer,,,)
  HAS_MEMBER_FUNCTION_VARIANT(data,_const, pointer,const,,)

  /// Case for if m_data has a member function called "data()", and is not a
  // string
  template<class U = T>
  typename std::enable_if<has_memberfunction_data<U>::value && !std::is_same<U,std::string>::value, rtype>::type
  data()
  {
#if CONTAINERARRAY_RETURN_PTR == 1
    return m_data->data();
#else
    return *m_data;
#endif
  }


  template<class U = T>
  typename std::enable_if<has_memberfunction_data<U>::value && !std::is_same<U,string>::value, rtype_const>::type
  data() const
  {
#if CONTAINERARRAY_RETURN_PTR == 1
    return m_data->data();
#else
    return *m_data;
#endif
  }


  /// Case for if m_data is a string
  template<class U = T>
  typename std::enable_if<std::is_same<U,std::string>::value, rtype>::type
  data()
  {
    /// return the object...or a reference to the object
    return *m_data;
  }


  template<class U = T>
  typename std::enable_if<std::is_same<U,std::string>::value, rtype_const>::type
  data() const
  {
    return *m_data;
  }


  /// case for if m_data does NOT have a member function "data()", and is not a
  // string
  template<class U = T>
  typename std::enable_if<!has_memberfunction_data<U>::value && !std::is_same<U,std::string>::value, rtype>::type
  data()
  {
    /// return a c-pointer to the object
    return m_data.get();
  }


  template<class U = T>
  typename std::enable_if<!has_memberfunction_data<U>::value && !std::is_same<U,std::string>::value, rtype_const>::type
  data() const
  {
    return m_data.get();
  }


  T& reference()
  { return *m_data; }

  T const & reference() const
  { return *m_data; }



  /// Case for if m_data has a member function called "data()"
  template<class U = T>
  typename std::enable_if< ( has_memberfunction_data<U>::value || has_memberfunction_v_const_data<U>::value ) &&
                           has_alias_pointer<U>::value && !std::is_same<U,string>::value,typename U::pointer >::type
  dataPtr()
  {
    return m_data->data();
  }

  template<class U = T>
  typename std::enable_if< ( has_memberfunction_data<U>::value || has_memberfunction_v_const_data<U>::value ) &&
                           has_alias_pointer<U>::value && !std::is_same<U,string>::value,typename U::const_pointer >::type
  dataPtr() const
  {
    return m_data->data();
  }


  /// Case for if m_data is a string"
  template<class U = T>
  typename std::enable_if< std::is_same<U,string>::value, char const * >::type
  dataPtr()
  {
    return m_data->data();
  }

  template<class U = T>
  typename std::enable_if< std::is_same<U,string>::value, char const * >::type
  dataPtr() const
  {
    return m_data->data();
  }


  /// case for if m_data does NOT have a member function "data()"
  template<class U = T>
  typename std::enable_if<!( has_memberfunction_data<U>::value || has_memberfunction_v_const_data<U>::value )&&
                          !std::is_same<U,string>::value, U * >::type
  dataPtr()
  {
    return m_data.get();
  }

  template<class U = T>
  typename std::enable_if<!( has_memberfunction_data<U>::value || has_memberfunction_v_const_data<U>::value )&&
                          !std::is_same<U,string>::value, U const *>::type
  dataPtr() const
  {
    return m_data.get();
  }

  HAS_ALIAS(value_type)


  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, size_t>::type
  sizeOfValueType() const
  {
    return sizeof(typename T::value_type);
  }

  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, size_t>::type
  sizeOfValueType() const
  {
    return sizeof(T);
  }

  virtual size_t sizeOfType() const override final
  {
    return sizeOfValueType();
  }


  /// case for if U::value_type exists. Returns the size of dataPtr
  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, localIndex>::type
  dataSize() const
  {
    return size() * sizeof(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the size of dataPtr
  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, localIndex>::type
  dataSize() const
  {
    return size() * sizeof(T);
  }


  /// case for if U::value_type exists. Returns the number of elements given a
  // byte size
  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, localIndex>::type
  numElementsFromDataSize(localIndex d_size) const
  {
    return d_size / sizeof(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the number of elements
  // given a byte size
  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, localIndex>::type
  numElementsFromDataSize(localIndex d_size) const
  {
    return d_size / sizeof(T);
  }

#ifdef USE_ATK
  /* Register the pointer to data with the associated sidre::View. */
  virtual void registerDataPtr() override final
  {
    localIndex d_size = dataSize();
    if (d_size > 0)
    {
      void * ptr = const_cast<void*>( static_cast<void const *>( dataPtr() ) );
      getSidreView()->setExternalDataPtr(axom::sidre::TypeID::INT8_ID, d_size, ptr);
    }
  }


  virtual void storeSizedFromParent() override final
  {
    getSidreView()->setAttributeScalar("__sizedFromParent__", sizedFromParent());
  }


  virtual void loadSizedFromParent() override final
  {
    setSizedFromParent(getSidreView()->getAttributeScalar("__sizedFromParent__"));
    getSidreView()->setAttributeToDefault("__sizedFromParent__");
  }


  virtual void unregisterDataPtr() override final
  {
    getSidreView()->setExternalDataPtr(AXOM_NULLPTR);
  }


  virtual void resizeFromSidre() override final
  {
    if (getSidreView()->isExternal())
    {
      localIndex d_size = getSidreView()->getTotalBytes();
      localIndex numElements = numElementsFromDataSize(d_size);
      resize(numElements);
    }

  }
  #endif /* ATK_FOUND */

  std::unique_ptr<T> m_data;

  ViewWrapper() = delete;
};

template< typename T >
using view_rtype = typename ViewWrapper<T>::rtype;

template< typename T >
using view_rtype_const = typename ViewWrapper<T>::rtype_const;

}
} /* namespace geosx */

#endif /* CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_ */
