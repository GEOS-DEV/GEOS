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
#include "common/Logger.hpp"
#include "SFINAE_Macros.hpp"
#include <type_traits>
#include "StringUtilities.hpp"
#include "Macros.hpp"
#include "Buffer.hpp"



#ifdef USE_ATK
#include "sidre/sidre.hpp"
#include "sidre/SidreTypes.hpp"
#include "SidreWrapper.hpp"
#endif

#include <cstdlib>


namespace geosx
{

/* Forward declarations */
class BasisBase;
class QuadratureBase;
class SimpleGeometricObjectBase;
class PartitionBase;


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
      GEOS_ERROR("invalid cast attempt");
#endif
    }
    return static_cast< ViewWrapper<T>& >(base);
  }


  struct empty_wrapper
  {
    HAS_MEMBER_FUNCTION(empty,bool,const,,)
    template<class U = T>
    static typename std::enable_if<has_memberfunction_empty<U>::value, bool>::type empty(ViewWrapper<T> const * parent)
    {
      return parent->m_data->empty();
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_empty<U>::value, bool>::type empty(ViewWrapper<T> const * parent)
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
    static typename std::enable_if< has_memberfunction_v0_size<U>::value ||
                                    has_memberfunction_v1_size<U>::value ||
                                    has_memberfunction_v2_size<U>::value ||
                                    has_memberfunction_v3_size<U>::value, localIndex>::type size(ViewWrapper<T> const * parent)
    {
      return static_cast<localIndex>(parent->m_data->size());
    }
    template<class U = T>
    static typename std::enable_if< !(has_memberfunction_v0_size<U>::value ||
                                      has_memberfunction_v1_size<U>::value ||
                                      has_memberfunction_v2_size<U>::value ||
                                      has_memberfunction_v3_size<U>::value ), localIndex>::type size(ViewWrapper<T> const * )
    {
      return 1;//parent->m_data;
    }
  };
  virtual localIndex size() const override final
  {
    return size_wrapper::size(this);
  }


  struct num_dimensions_wrapper
  {
    HAS_MEMBER_FUNCTION(numDimensions,int,const,,)

    template<class U = T>
    static typename std::enable_if<has_memberfunction_numDimensions<U>::value, localIndex>::type
    numDimensions(ViewWrapper<T> const * parent)
    { return static_cast<localIndex>(parent->m_data->numDimensions()); }
    
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_numDimensions<U>::value, localIndex>::type
    numDimensions(ViewWrapper<T> const * parent)
    { return 1; }
  };
  virtual localIndex numDimensions() const override final
  {
    return num_dimensions_wrapper::numDimensions(this);
  }


  struct dimension_wrapper
  {
    HAS_MEMBER_FUNCTION(Dimension,long,const, VA_LIST(long), VA_LIST(long(1)))

    template<class U = T>
    static typename std::enable_if<has_memberfunction_Dimension<U>::value, localIndex>::type
    dimension(ViewWrapper<T> const * parent, localIndex i)
    { return static_cast<localIndex>(parent->m_data->Dimension(i)); }
    
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_Dimension<U>::value, localIndex>::type
    dimension(ViewWrapper<T> const * parent, localIndex i)
    { 
      if (i != 0) 
      {
        GEOS_ERROR("Data is only 1D");
        return 0;
      }
      return parent->size(); 
    }
  };
  virtual localIndex dimension(localIndex i) const override final
  {
    return dimension_wrapper::dimension(this, i);
  }


  struct set_dimension_wrapper
  {
    HAS_MEMBER_FUNCTION(setDimensions, void, , VA_LIST(long, const long *), VA_LIST(long(1), nullptr))

    template<class U=T>
    static typename std::enable_if<has_memberfunction_setDimensions<U>::value, void>::type
    setDimensions(ViewWrapper<T> * parent, const long num_dims, const long * dims)
    { parent->m_data->setDimensions(num_dims, dims); }

    template<class U=T>
    static typename std::enable_if<!has_memberfunction_setDimensions<U>::value, void>::type
    setDimensions(ViewWrapper<T> * parent, const long num_dims, const long * dims)
    {
      if (num_dims != 1)
      {
        GEOS_ERROR("Data is only 1D");
        return;
      }
      parent->resize(dims[0]);
    }
  };
  virtual void setDimensions(long num_dims, const long * dims)
  { set_dimension_wrapper::setDimensions(this, num_dims, dims); }


  struct reserve_wrapper
  {
    HAS_MEMBER_FUNCTION(reserve, void, ,VA_LIST(std::size_t),VA_LIST(std::size_t(1)) )
    template<class U = T>
    static typename std::enable_if<has_memberfunction_reserve<U>::value, void>::type reserve(ViewWrapper<T> * const parent, std::size_t new_cap)
    {
      return parent->m_data->reserve(new_cap);
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_reserve<U>::value, void>::type reserve(ViewWrapper<T> * const, std::size_t )
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
    static typename std::enable_if<has_memberfunction_v0_resize<U>::value, void>::type
    resize(ViewWrapper<T> * const parent, int32 const new_size)
    {
      return parent->m_data->resize(new_size);
    }

    template<class U = T>
    static typename std::enable_if<has_memberfunction_v1_resize<U>::value, void>::type
    resize(ViewWrapper<T> * const parent, uint32 const new_size)
    {
      return parent->m_data->resize(new_size);
    }

    template<class U = T>
    static typename std::enable_if<has_memberfunction_v2_resize<U>::value, void>::type
    resize(ViewWrapper<T> * const parent, int64 const new_size)
    {
      return parent->m_data->resize(new_size);
    }

    template<class U = T>
    static typename std::enable_if<has_memberfunction_v3_resize<U>::value, void>::type
    resize(ViewWrapper<T> * const parent, uint64 const new_size)
    {
      return;
    }
  };
  using ViewWrapperBase::resize;
  void resize( localIndex new_size ) override final
  {
    resize_wrapper::resize(this, new_size);
  }


  struct should_resize_wrapper
  {
    HAS_MEMBER_FUNCTION(isSorted,bool,const,,)
    template<class U = T>
    static typename std::enable_if<has_memberfunction_isSorted<U>::value, bool>::type shouldResize()
    { return false;  }

    template<class U = T>
    static typename std::enable_if<!has_memberfunction_isSorted<U>::value, bool>::type shouldResize()
    { return true; }
  };
  virtual bool shouldResize() const override final
  {
    return should_resize_wrapper::shouldResize();
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

  /// case for if U::value_type exists. Returns the size of dataPtr
  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, localIndex>::type
  byteSize() const
  {
    return size() * sizeof(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the size of dataPtr
  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, localIndex>::type
  byteSize() const
  {
    return size() * sizeof(T);
  }


  /// case for if U::value_type exists. Returns the size of an element of dataPtr
  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, localIndex>::type
  elementSize() const
  {
    return sizeof(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the size of an element of dataPtr
  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, localIndex>::type
  elementSize() const
  {
    return sizeof(T);
  }


  /// case for if U::value_type exists. Returns the typeid of an element of dataPtr
  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, const std::type_info&>::type
  elementTypeID() const
  {
    return typeid(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the typeid of an element of dataPtr
  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, const std::type_info&>::type
  elementTypeID() const
  {
    return typeid(T);
  }


  /// case for if U::value_type exists. Returns the number of elements given a byte size
  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, localIndex>::type
  numElementsFromByteSize(localIndex d_size) const
  {
    return d_size / sizeof(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the number of elements
  // given a byte size
  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, localIndex>::type
  numElementsFromByteSize(localIndex d_size) const
  {
    return d_size / sizeof(T);
  }


  struct register_to_write_wrapper
  {
    template<class U = T>
    static typename std::enable_if<(std::is_base_of<SimpleGeometricObjectBase, U>::value ||
                             std::is_base_of<QuadratureBase, U>::value ||
                             std::is_base_of<PartitionBase, U>::value ||
                             std::is_base_of<BasisBase, U>::value), void>::type
    registerToWrite(const ViewWrapper<T> * parent, axom::sidre::View * view)
    { 
#ifdef USE_ATK
      view = (view != nullptr) ? view : parent->getSidreView();
      parent->storeSizedFromParent(view);
      parent->unregisterDataPtr(view); 
#endif
    }


    /* Register the pointer to data with the associated sidre::View. */
    template<class U = T>
    static typename std::enable_if<!(std::is_base_of<SimpleGeometricObjectBase, U>::value ||
                                     std::is_base_of<QuadratureBase, U>::value ||
                                     std::is_base_of<PartitionBase, U>::value ||
                                     std::is_base_of<BasisBase, U>::value), void>::type
    registerToWrite(const ViewWrapper<T> * parent, axom::sidre::View * view)
    {
#ifdef USE_ATK
      view = (view != nullptr) ? view : parent->getSidreView();
      parent->storeSizedFromParent(view);

      localIndex num_elements = parent->size();
      if (num_elements > 0) 
      {
        std::type_index type_index = std::type_index(parent->elementTypeID());
        axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType(type_index);
        if (sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID)
        {
          localIndex byte_size;
          void * ptr = Buffer::pack(parent->reference(), byte_size);
          view->setExternalDataPtr(axom::sidre::TypeID::INT8_ID, byte_size, ptr);
          return;
        }

        localIndex sidre_size = rtTypes::getSidreSize(type_index);
        localIndex byte_size = parent->byteSize();
        localIndex element_size = parent->elementSize();

        int ndims = parent->numDimensions();
        axom::sidre::SidreLength dims[ndims + 1];
        for (localIndex dim = 0; dim < ndims; ++dim)
        {
          dims[dim] = parent->dimension(dim);
        }

        if ( byte_size > num_elements * sidre_size )
        {
          dims[ndims++] = element_size / sidre_size;
        }
        
        void * ptr = const_cast<void*>((void const *) parent->dataPtr());
        view->setExternalDataPtr(sidre_type_id, ndims, dims, ptr);
      }
      else
      {
        parent->unregisterDataPtr(view);
      }
#endif
    }
  };
  void registerToWrite(axom::sidre::View * view=nullptr) const override
  { register_to_write_wrapper::registerToWrite(this, view); }


  struct finish_writing_wrapper
  {
    template<class U = T>
    static typename std::enable_if<(std::is_base_of<SimpleGeometricObjectBase, U>::value ||
                             std::is_base_of<QuadratureBase, U>::value ||
                             std::is_base_of<PartitionBase, U>::value ||
                             std::is_base_of<BasisBase, U>::value), void>::type
    finishWriting(const ViewWrapper<T> * parent, axom::sidre::View * view)
    { 
#ifdef USE_ATK
      view = (view != nullptr) ? view : parent->getSidreView();
      view->setAttributeToDefault("__sizedFromParent__");
      parent->unregisterDataPtr(view);
#endif
    }


    /* Register the pointer to data with the associated sidre::View. */
    template<class U = T>
    static typename std::enable_if<!(std::is_base_of<SimpleGeometricObjectBase, U>::value ||
                              std::is_base_of<QuadratureBase, U>::value ||
                              std::is_base_of<PartitionBase, U>::value ||
                              std::is_base_of<BasisBase, U>::value), void>::type
    finishWriting(const ViewWrapper<T> * parent, axom::sidre::View * view)
    {
#ifdef USE_ATK
      view = (view != nullptr) ? view :parent->getSidreView();
      view->setAttributeToDefault("__sizedFromParent__");

      if (!view->isExternal() || view->getTotalBytes() == 0)
      {
        return;
      }
      
      std::type_index type_index = std::type_index(parent->elementTypeID());
      axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType(type_index);
      if (sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID)
      {
        std::free(view->getVoidPtr());
      }

      parent->unregisterDataPtr(view);
#endif
    }
  };
  void finishWriting(axom::sidre::View * view=nullptr) const override
  { finish_writing_wrapper::finishWriting(this, view); }



  struct register_to_read_wrapper
  {
    template<class U = T>
    static typename std::enable_if<(std::is_base_of<SimpleGeometricObjectBase, U>::value ||
                             std::is_base_of<QuadratureBase, U>::value ||
                             std::is_base_of<PartitionBase, U>::value ||
                             std::is_base_of<BasisBase, U>::value), void>::type
    registerToRead(ViewWrapper<T> * parent, axom::sidre::View * view=nullptr)
    { 
#ifdef USE_ATK
      parent->loadSizedFromParent(view);
      parent->unregisterDataPtr(view);
#endif
    }


    /* Register the pointer to data with the associated sidre::View. */
    template<class U = T>
    static typename std::enable_if<!(std::is_base_of<SimpleGeometricObjectBase, U>::value ||
                              std::is_base_of<QuadratureBase, U>::value ||
                              std::is_base_of<PartitionBase, U>::value ||
                              std::is_base_of<BasisBase, U>::value), void>::type
    registerToRead(ViewWrapper<T> * parent, axom::sidre::View * view = nullptr)
    {
#ifdef USE_ATK
      view = (view != nullptr) ? view : parent->getSidreView();
      parent->loadSizedFromParent(view);
      if (!view->isExternal() || view->getTotalBytes() == 0)
      {
        return;
      }
      
      std::type_index type_index = std::type_index(parent->elementTypeID());
      axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType(type_index);
      if (sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID)
      {
        localIndex byte_size = view->getTotalBytes();
        void * ptr = std::malloc(byte_size);
        view->setExternalDataPtr(axom::sidre::TypeID::INT8_ID, byte_size, ptr);
        return;
      }

      parent->resizeFromSidre(view);
      void * ptr = const_cast<void*>((void const *) parent->dataPtr());
      localIndex sidre_size = rtTypes::getSidreSize(type_index);
      view->setExternalDataPtr(sidre_type_id, parent->byteSize() / sidre_size, ptr);
#endif
    }
  };
  void registerToRead(axom::sidre::View * view=nullptr) override
  { register_to_read_wrapper::registerToRead(this, view); }



  struct finish_reading_wrapper
  {
    template<class U = T>
    static typename std::enable_if<(std::is_base_of<SimpleGeometricObjectBase, U>::value ||
                                    std::is_base_of<QuadratureBase, U>::value ||
                                    std::is_base_of<PartitionBase, U>::value ||
                                    std::is_base_of<BasisBase, U>::value), void>::type
    finishReading(ViewWrapper<T> * parent, axom::sidre::View * view)
    { 
#ifdef USE_ATK
      view = (view != nullptr) ? view : parent->getSidreView();
      parent->unregisterDataPtr(view);
#endif
    }


    /* Register the pointer to data with the associated sidre::View. */
    template<class U = T>
    static typename std::enable_if<!(std::is_base_of<SimpleGeometricObjectBase, U>::value ||
                                     std::is_base_of<QuadratureBase, U>::value ||
                                     std::is_base_of<PartitionBase, U>::value ||
                                     std::is_base_of<BasisBase, U>::value), void>::type
    finishReading(ViewWrapper<T> * parent, axom::sidre::View * view)
    {
#ifdef USE_ATK
      view = (view != nullptr) ? view : parent->getSidreView();
      if (!view->isExternal() || view->getTotalBytes() == 0)
      {
        return;
      }
      
      std::type_index type_index = std::type_index(parent->elementTypeID());
      axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType(type_index);
      if (sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID)
      {
        localIndex byte_size = view->getTotalBytes();
        void * ptr = view->getVoidPtr();
        Buffer::unpack(parent->reference(), ptr, byte_size);
        std::free(ptr);
      }

      parent->unregisterDataPtr(view);
#endif
    }
  };
  void finishReading(axom::sidre::View * view = nullptr) override
  { finish_reading_wrapper::finishReading(this, view); }


  void unregisterDataPtr(axom::sidre::View* view = nullptr) const
  {
#if ATK_FOUND
    view = (view != nullptr) ? view : getSidreView();
    view->setExternalDataPtr(AXOM_NULLPTR);
#endif
  }

  void storeSizedFromParent(axom::sidre::View* view = nullptr) const
  {
#ifdef USE_ATK
    if (SidreWrapper::dataStore().hasAttribute("__sizedFromParent__"))
    {
      view = (view != nullptr) ? view : getSidreView();
      view->setAttributeScalar("__sizedFromParent__", sizedFromParent());
    }
#endif
  }


  void loadSizedFromParent(axom::sidre::View* view = nullptr)
  {
#ifdef USE_ATK
    if (SidreWrapper::dataStore().hasAttribute("__sizedFromParent__"))
    {
      view = (view != nullptr) ? view : getSidreView();
      setSizedFromParent(view->getAttributeScalar("__sizedFromParent__"));
      view->setAttributeToDefault("__sizedFromParent__");
    }
#endif
  }


  void resizeFromSidre(axom::sidre::View* view = nullptr)
  {
#ifdef USE_ATK
    view = (view != nullptr) ? view : getSidreView();
    if (view->isExternal()) 
    { 
      std::type_index type_index = std::type_index(elementTypeID());
      localIndex sidre_size = rtTypes::getSidreSize(type_index);

      localIndex byte_size = view->getTotalBytes();
      localIndex num_elements = numElementsFromByteSize(byte_size);

      int ndims = view->getNumDimensions();
      axom::sidre::SidreLength dims[ndims];
      view->getShape(ndims, dims);

      if ( byte_size > num_elements * sidre_size )
      {
        ndims--;
      }

      localIndex num_elems_recorded = 1;
      for (localIndex i = 0; i < ndims; ++i)
      {
        num_elems_recorded *= dims[i];
      }

      if (num_elems_recorded != num_elements)
      {
        GEOS_ERROR("Number of elements recorded not equal to the calculated number: " << 
                   num_elems_recorded << " " << num_elements);
      }

      setDimensions(ndims, dims);
    }
#endif
  }


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
