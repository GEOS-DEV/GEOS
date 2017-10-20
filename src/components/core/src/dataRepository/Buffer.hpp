#ifndef DATAREPOSITORY_BUFFERHPP
#define DATAREPOSITORY_BUFFERHPP


#include "common/DataTypes.hpp"
#include "common/Logger.hpp"
#include <vector>
#include <cstdlib>
#include <string>
#include <utility>
#include <type_traits>



namespace geosx
{
namespace dataRepository
{

class Buffer
{


  template <typename T>
  static T* allocBuffer(localIndex byte_size)
  {
    void* buff = std::malloc(byte_size);
    if (buff == nullptr)
    {
      GEOS_ERROR("Allocation error");
    }

    m_buffers.push_back(buff);
    return static_cast<T*>(buff);
  }


  static void clear()
  {
    for (void * buff : m_buffers)
    {
      std::free(buff);
    }
    m_buffers.clear();
  }

  HAS_ALIAS(value_type)



  /* Arithmetic type
   * Format:
   *   T data
   */
  template <typename T>
  static typename std::enable_if<std::is_arithmetic<T>::value, localIndex>::type
  packed_size(const T & value)
  {
    return sizeof(T);
  }
  

  template <typename T>
  static typename std::enable_if<std::is_arithmetic<T>::value, void *>::type
  pack(const T & value, localIndex & byte_size, void * buffer=nullptr)
  {
    byte_size = packed_size(arr);
    T * buff = static_cast<T *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<T>(byte_size);
    }

    buff[0] = value;
    return buff;
  }
  

  template <typename T>
  static typename std::enable_if<std::is_arithmetic<T>::value, localIndex>::type
  unpack(T & value, const void * buffer, localIndex byte_size=-1)
  {
    T * buff = static_cast<T *>(buffer);
    value = buff[0];
    return sizeof(T);
  }




  /* TensorBaseT<T_DIM> type
   * Format:
   *   realT[T_DIM] data
   */
  template <int T_DIM>
  static localIndex packed_size(const TensorBaseT<T_DIM> & t)
  {
    return T_DIM * sizeof(realT);
  }
  

  template <int T_DIM>
  static void * pack(const TensorBaseT<T_DIM> & t, localIndex & byte_size, void * buffer=nullptr)
  {
    byte_size = packed_size(arr);
    realT * buff = static_cast<realT *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<realT>(byte_size);
    }

    std::memcpy(buffer, t.Data(), byte_size);
    return buff;
  }
  

  template <int T_DIM>
  static localIndex unpack(TensorBaseT<T_DIM> & t, const void * buffer, localIndex byte_size=-1)
  {
    const realT * buff = static_cast<const realT *>(buffer);
    byte_size = packed_size(t)
    std::memcpy(t.Data(), buff, byte_size);
    return byte_size;
  }




  /* String
   * Format:
   *   char[] null terminated string
   */
  static localIndex packed_size(const string & s)
  {
    return (s.length() + 1) * sizeof(char);
  }
  

  static void * pack(const string & s, localIndex & byte_size, void * buffer=nullptr)
  {
    byte_size = packed_size(s);
    T * buff = static_cast<char *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<char>(byte_size);
    }

    std::memcpy(buff, s.data(), s.length() * sizeof(char));
    buff[s.length()] = '\0';
    return buff;
  }
  

  template <typename T>
  static localIndex unpack(string & s, const void * buffer, localIndex byte_size=-1)
  {
    s = string(static_cast<char *>(buffer));
    return s.length() + 1;
  }




  /* 1D array of objects.
   * Format:
   *   localIndex the total number of bytes
   *   localIndex the length of the array
   *   pack(a[0])
   *   pack(a[1])
   *      .
   *      .
   *      .
   *   pack(a[N-1])
   */
  template <typename T>
  static localIndex packed_size(const array<T> & arr)
  {
    localIndex byte_size = 2 * sizeof(localIndex);
    for (const T & elem : arr)
    {
      byte_size += packed_size(a);
    }

    return byte_size;
  }


  template <typename T>
  static void * pack(const array<T> & arr, localIndex & byte_size, void * buffer=nullptr)
  {
    byte_size = packed_size(arr);
    localIndex * buff = static_cast<localIndex *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<localIndex>(byte_size);
    }

    buff[0] = byte_size;
    buff[1] = arr.size();
    void * v_buff = static_cast<void *>(buff + 2);

    localIndex offset = 0;
    localIndex bytes_written;
    for (const T & elem : arr)
    {
      pack(elem, bytes_written, v_buff + offset);
      offset += bytes_written;
    }

    return buff;
  }


  template <typename T>
  static localIndex unpack(array<T> & arr, const void * buffer, localIndex byte_size=-1)
  {
    const localIndex * buff = static_cast<localIndex *> buffer;
    localIndex bytes_recorded = buff[0];
    localIndex num_arrays = buff[1];
    arr.resize(num_arrays);

    if (bytes_recorded != byte_size && byte_size >= 0)
    {
      GEOS_ERROR("Number of bytes recorded not equal to number of bytes in buffer: " <<
                 bytes_recorded << " " << byte_size);
    }

    const void* v_buff = static_cast<void *>(buff + 2);
    localIndex offset = 0;
    for (localIndex i = 0 ; i < num_arrays; ++i)
    {
      offset += unpack(arr[i], v_buff + offset);      
    }

    localIndex bytes_read = offset + 2 * sizeof(localIndex);
    if (bytes_read != bytes_recorded)
    {
      GEOS_ERROR("Number of bytes read not equal to number of bytes in recorded: " <<
                 bytes_read << " " << bytes_recorded);
    }

    return bytes_read;
  }




  /* 2D array of plain old data.
   * Format:
   *   localIndex size of the first dimension
   *   localIndex size of the second dimension
   *   T[] data
   */
  template <typename T>
  static typename std::enable_if<!has_alias_value_type<T>::value, localIndex>::type
  packed_size(const Array2dT<T> & arr)
  {
    localIndex byte_size = 2 * sizeof(localIndex);
    byte_size += arr.size() * sizeof(T);
    return byte_size;
  }


  template <typename T>
  static typename std::enable_if<!has_alias_value_type<T>::value, void *>::type
  pack(const Array2dT<T> & arr, localIndex & byte_size, void * buffer=nullptr)
  {
    byte_size = packed_size(arr);
    localIndex * buff = static_cast<localIndex *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<localIndex>(byte_size);
    }

    buff[0] = arr.Dimension(0);
    buff[1] = arr.Dimension(1);

    std::memcpy(buff + 2, arr.data(), arr.size() * sizeof(T));
    return buff;
  }


  template <typename T>
  static typename std::enable_if<!has_alias_value_type<T>::value, localIndex>::type
  unpack(Array2dT<T> & arr, const void * buffer, localIndex byte_size=-1)
  {
    const localIndex * buff = static_cast<localIndex *> buffer;
    const localIndex dim0 = buff[0];
    const localIndex dim1 = buff[0];
    arr.resize2(dim0, dim1);

    localIndex bytes_recorded = 2 * sizeof(localIndex);
    bytes_recorded += dim0 * dim1 * sizeof(T);
    if (bytes_recorded != byte_size && byte_size >= 0)
    {
      GEOS_ERROR("Number of bytes recorded not equal to number of bytes in buffer: " <<
                 bytes_recorded << " " << byte_size);
    }

    std::memcpy(arr.data(), buff + 2, dim0 * dim1 * sizeof(T));
    return bytes_read;
  }




  /* Pair
   * Format:
   *   localIndex the number of bytes in pack(T)
   *   localIndex the number of bytes in pack(V)
   *   pack(T)
   *   pack(V)
   */
  template <typename T, V>
  static localIndex packed_size(const std::pair<T, V> & p)
  {
    localIndex byte_size = 2 * sizeof(localIndex);
    byte_size += packed_size(p.first) + packed_size(p.second);
    return byte_size;
  }


  template <typename T, V>
  static void * pack(const std::pair<T, V> & p, localIndex & byte_size, void * buffer=nullptr)
  {
    localIndex byte_size_first = packed_size(p.first);
    localIndex byte_size_second = packed_size(p.second);
    byte_size = 2 * sizeof(localIndex) + byte_size_first + byte_size_second;

    localIndex * buff = static_cast<localIndex *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<localIndex>(byte_size);
    }

    buff[0] = byte_size_first;
    buff[1] = byte_size_second;
    void * v_buff = static_cast<void *>(buff + 2);

    pack(p.first, byte_size_first, v_buff);
    pack(p.second, byte_size_second, v_buff + byte_size_first);
    return buff;
  }


  template <typename T, V>
  static localIndex unpack(std::pair<T, V> & p, const void * buffer, localIndex byte_size=-1)
  {
    const localIndex * buff = static_cast<localIndex *> buffer;
    localIndex byte_size_first = buff[0];
    localIndex byte_size_second = buff[1];
    localIndex bytes_recorded = 2 * sizeof(localIndex) + byte_size_first + byte_size_second;

    if (bytes_recorded != byte_size && byte_size >= 0)
    {
      GEOS_ERROR("Number of bytes read not equal to number of bytes in buffer: " <<
                 bytes_recorded << " " << byte_size);
    }

    const void * v_buff = static_cast<const void *>(buff + 2);
    unpack(p.first, v_buff, byte_size_first);
    unpack(p.second, v_buff + byte_size_first, byte_size_second);

    return bytes_recorded;
  }



  /* Map
   * Format:
   *   localIndex the total number of bytes
   *   localIndex the number of key-value pairs
   *   pack(K0)
   *   pack(V0)
   *   pack(K1)
   *   pack(V1)
   *      .
   *      .
   *      .
   *   pack(KN)
   *   pack(VN)
   */
  template <typename K, V>
  static localIndex packed_size(const map<K, V> & m)
  {
    localIndex byte_size = 2 * sizeof(localIndex);
    for (map<K, V>::value_type & p : m)
    {
      byte_size += packed_size(p.first);
      byte_size += packed_size(p.second);
    }
    return byte_size;
  }


  template <typename K, V>
  static void * pack(const map<K, V> & m, localIndex & byte_size, void * buffer=nullptr)
  {
    byte_size = packed_size(m);
    localIndex * buff = static_cast<localIndex *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<localIndex>(byte_size);
    }

    buff[0] = byte_size;
    buff[1] = m.size();
    void * v_buff = static_cast<void *>(buff + 2);
    
    localIndex offset = 0;
    localIndex prev_size;
    for (map<K, V>::value_type & p : m)
    {
      pack(p.first, prev_size, v_buff + offset);
      offset += prev_size;
      pack(p.second, prev_size, v_buff + offset);
      offset += prev_size;
    }

    return buff;
  }


  template <typename K, V>
  static localIndex unpack(map<K, V> & m, const void * buffer, localIndex byte_size=-1)
  {
    const localIndex * buff = static_cast<localIndex *> buffer;
    localIndex bytes_recorded = buff[0];
    localIndex num_pairs = buff[1];

    if (bytes_recorded != byte_size && byte_size >= 0)
    {
      GEOS_ERROR("Number of bytes recorded not equal to number of bytes in buffer: " <<
                  bytes_recorded << " " << byte_size);
    }

    localIndex offset = 0;
    const void * v_buff = static_cast<const void *>(buff + 2);
    for (localIndex i = 0; i < num_pairs; ++i)
    {
      K key;
      offset += unpack(key, v_buff + offset);
      V value;
      offset += unpack(value, v_buff + offset);
      m[key] = std::move(value);
    }

    localIndex bytes_read = 2 * sizeof(localIndex) + offset;
    if (bytes_read != bytes_recorded && byte_size >= 0)
    {
      GEOS_ERROR("Number of bytes read not equal to number of bytes recorded: " <<
                  bytes_read << " " << bytes_recorded);
    }

    return bytes_read;
  }

  /* Arithmetic type
   * Format:
   *   T data
   */
  static localIndex packed_size(const SimpleGeometricObjectBase & value)
  { return 0; }
  

  pack(const T & value, localIndex & byte_size, void * buffer=nullptr)
  { 
    byte_size = packed_size;
    return buffer; 
  }
  

  template <typename T>
  static typename std::enable_if<std::is_arithmetic<T>::value, localIndex>::type
  unpack(T & value, const void * buffer, localIndex byte_size=-1)
  {
    T * buff = static_cast<T *>(buffer);
    value = buff[0];
    return sizeof(T);
  }


private:
  Buffer();
  ~Buffer();

  static std::vector<void *> m_buffers;  
}




}     /* end namespace dataRepository */
}     /* end namespace geosx */


#endif /* DATAREPOSITORY_BUFFERHPP */