// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
#ifndef DATAREPOSITORY_BUFFERHPP
#define DATAREPOSITORY_BUFFERHPP


#include "common/DataTypes.hpp"
#include "common/Logger.hpp"
#include "SFINAE_Macros.hpp"
#include <vector>
#include <cstdlib>
#include <string>
#include <utility>
#include <type_traits>
#include <memory>



namespace geosx
{

/* Forward declarations */
class BasisBase;
class QuadratureBase;
class SimpleGeometricObjectBase;
class PartitionBase;
class NeighborCommunicator;

namespace dataRepository
{

class Buffer
{
public:   

  template <typename T>
  static T* allocBuffer(localIndex byte_size)
  {
    void* buff = std::malloc(byte_size);
    if (buff == nullptr)
    {
      GEOS_ERROR("Allocation error");
    }

    return reinterpret_cast<T*>(buff);
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
    byte_size = packed_size(value);
    T * buff = reinterpret_cast<T *>(buffer);
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
    const T * buff = reinterpret_cast<const T *>(buffer);
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
    byte_size = packed_size(t);
    realT * buff = reinterpret_cast<realT *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<realT>(byte_size);
    }

    std::memcpy(buff, t.Data(), byte_size);
    return buff;
  }
  

  template <int T_DIM>
  static localIndex unpack(TensorBaseT<T_DIM> & t, const void * buffer, localIndex byte_size=-1)
  {
    const realT * buff = reinterpret_cast<const realT *>(buffer);
    byte_size = packed_size(t);
    std::memcpy(t.Data(), buff, byte_size);
    return byte_size;
  }




  /* String
   * Format:
   *   char[] null terminated string
   */
  static localIndex packed_size(const string & s)
  {
    return ( static_cast<localIndex>(s.length()) + 1) * sizeof(char);
  }
  

  static void * pack(const string & s, localIndex & byte_size, void * buffer=nullptr)
  {
    byte_size = packed_size(s);
    char * buff = reinterpret_cast<char *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<char>(byte_size);
    }

    std::memcpy(buff, s.data(), s.length() * sizeof(char));
    buff[s.length()] = '\0';
    return buff;
  }
  

  static localIndex unpack(string & s, const void * buffer, localIndex byte_size=-1)
  {
    s = string(reinterpret_cast<const char *>(buffer));
    return static_cast<localIndex>(s.length()) + 1;
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
      byte_size += packed_size(elem);
    }

    return byte_size;
  }


  template <typename T>
  static void * pack(const array<T> & arr, localIndex & byte_size, void * buffer=nullptr)
  {
    byte_size = packed_size(arr);
    localIndex * buff = reinterpret_cast<localIndex *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<localIndex>(byte_size);
    }

    buff[0] = byte_size;
    buff[1] = arr.size();
    char * c_buff = reinterpret_cast<char *>(buff + 2);

    localIndex offset = 0;
    localIndex bytes_written;
    for (const T & elem : arr)
    {
      pack(elem, bytes_written, c_buff + offset);
      offset += bytes_written;
    }

    return buff;
  }


  template <typename T>
  static localIndex unpack(array<T> & arr, const void * buffer, localIndex byte_size=-1)
  {
    const localIndex * buff = reinterpret_cast<const localIndex *>(buffer);
    localIndex bytes_recorded = buff[0];
    localIndex num_arrays = buff[1];
    arr.resize(num_arrays);

    if (bytes_recorded != byte_size && byte_size >= 0)
    {
      GEOS_ERROR("Number of bytes recorded not equal to number of bytes in buffer: " <<
                 bytes_recorded << " " << byte_size);
    }

    const char * c_buff = reinterpret_cast<const char *>(buff + 2);
    localIndex offset = 0;
    for (localIndex i = 0 ; i < num_arrays; ++i)
    {
      offset += unpack(arr[i], c_buff + offset);      
    }

    localIndex bytes_read = offset + 2 * sizeof(localIndex);
    if (bytes_read != bytes_recorded)
    {
      GEOS_ERROR("Number of bytes read not equal to number of bytes in recorded: " <<
                 bytes_read << " " << bytes_recorded);
    }

    return bytes_read;
  }




  /* Set of objects.
   * Format:
   *   localIndex the total number of bytes
   *   localIndex the number of items in the set
   *   pack(a[0])
   *   pack(a[1])
   *      .
   *      .
   *      .
   *   pack(a[N-1])
   */
  template <typename T>
  static localIndex packed_size(const set<T> & s)
  {
    localIndex byte_size = 2 * sizeof(localIndex);
    for (const T & elem : s)
    {
      byte_size += packed_size(elem);
    }

    return byte_size;
  }


  template <typename T>
  static void * pack(const set<T> & s, localIndex & byte_size, void * buffer=nullptr)
  {
    byte_size = packed_size(s);
    localIndex * buff = reinterpret_cast<localIndex *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<localIndex>(byte_size);
    }

    buff[0] = byte_size;
    buff[1] = s.size();
    char * c_buff = reinterpret_cast<char *>(buff + 2);

    localIndex offset = 0;
    localIndex bytes_written;
    for (const T & elem : s)
    {
      pack(elem, bytes_written, c_buff + offset);
      offset += bytes_written;
    }

    return buff;
  }


  template <typename T>
  static localIndex unpack(set<T> & s, const void * buffer, localIndex byte_size=-1)
  {
    const localIndex * buff = reinterpret_cast<const localIndex *>(buffer);
    localIndex bytes_recorded = buff[0];
    localIndex num_arrays = buff[1];
    s.resize(num_arrays);

    if (bytes_recorded != byte_size && byte_size >= 0)
    {
      GEOS_ERROR("Number of bytes recorded not equal to number of bytes in buffer: " <<
                 bytes_recorded << " " << byte_size);
    }

    const char * c_buff = reinterpret_cast<const char *>(buff + 2);
    localIndex offset = 0;
    for (localIndex i = 0 ; i < num_arrays; ++i)
    {
      offset += unpack(s[i], c_buff + offset);      
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
    localIndex * buff = reinterpret_cast<localIndex *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<localIndex>(byte_size);
    }

    buff[0] = arr.size(0);
    buff[1] = arr.size(1);

    std::memcpy(buff + 2, arr.data(), arr.size() * sizeof(T));
    return buff;
  }


  template <typename T>
  static typename std::enable_if<!has_alias_value_type<T>::value, localIndex>::type
  unpack(Array2dT<T> & arr, const void * buffer, localIndex byte_size=-1)
  {
    const localIndex * buff = reinterpret_cast<const localIndex *>(buffer);
    const localIndex dim0 = buff[0];
    const localIndex dim1 = buff[1];
    arr.resize(dim0, dim1);

    localIndex bytes_recorded = 2 * sizeof(localIndex);
    bytes_recorded += dim0 * dim1 * sizeof(T);
    if (bytes_recorded != byte_size && byte_size >= 0)
    {
      GEOS_ERROR("Number of bytes recorded not equal to number of bytes in buffer: " <<
                 bytes_recorded << " " << byte_size);
    }

    std::memcpy(arr.data(), buff + 2, dim0 * dim1 * sizeof(T));
    return bytes_recorded;
  }




  /* Pair
   * Format:
   *   localIndex the number of bytes in pack(T)
   *   localIndex the number of bytes in pack(V)
   *   pack(T)
   *   pack(V)
   */
  template <typename T, typename V>
  static localIndex packed_size(const std::pair<T, V> & p)
  {
    localIndex byte_size = 2 * sizeof(localIndex);
    byte_size += packed_size(p.first) + packed_size(p.second);
    return byte_size;
  }


  template <typename T, typename V>
  static void * pack(const std::pair<T, V> & p, localIndex & byte_size, void * buffer=nullptr)
  {
    localIndex byte_size_first = packed_size(p.first);
    localIndex byte_size_second = packed_size(p.second);
    byte_size = 2 * sizeof(localIndex) + byte_size_first + byte_size_second;

    localIndex * buff = reinterpret_cast<localIndex *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<localIndex>(byte_size);
    }

    buff[0] = byte_size_first;
    buff[1] = byte_size_second; 
    char * c_buff = reinterpret_cast<char *>(buff + 2);

    pack(p.first, byte_size_first, c_buff);
    pack(p.second, byte_size_second, c_buff + byte_size_first);
    return buff;
  }


  template <typename T, typename V>
  static localIndex unpack(std::pair<T, V> & p, const void * buffer, localIndex byte_size=-1)
  {
    const localIndex * buff = reinterpret_cast<const localIndex *>(buffer);
    localIndex byte_size_first = buff[0];
    localIndex byte_size_second = buff[1];
    localIndex bytes_recorded = 2 * sizeof(localIndex) + byte_size_first + byte_size_second;

    if (bytes_recorded != byte_size && byte_size >= 0)
    {
      GEOS_ERROR("Number of bytes read not equal to number of bytes in buffer: " <<
                 bytes_recorded << " " << byte_size);
    }

    const char * c_buff = reinterpret_cast<const char *>(buff + 2);
    unpack(p.first, c_buff, byte_size_first);
    unpack(p.second, c_buff + byte_size_first, byte_size_second);

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
  template <typename K, typename V>
  static localIndex packed_size(const map<K, V> & m)
  {
    localIndex byte_size = 2 * sizeof(localIndex);
    for (typename map<K, V>::value_type const & p : m)
    {
      byte_size += packed_size(p.first);
      byte_size += packed_size(p.second);
    }
    return byte_size;
  }


  template <typename K, typename V>
  static void * pack(const map<K, V> & m, localIndex & byte_size, void * buffer=nullptr)
  {
    byte_size = packed_size(m);
    localIndex * buff = reinterpret_cast<localIndex *>(buffer);
    if (buff == nullptr)
    {
      buff = allocBuffer<localIndex>(byte_size);
    }

    buff[0] = byte_size;
    buff[1] = m.size();
    char * c_buff = reinterpret_cast<char *>(buff + 2);
    
    localIndex offset = 0;
    localIndex prev_size;
    for (typename map<K, V>::value_type const & p : m)
    {
      pack(p.first, prev_size, c_buff + offset);
      offset += prev_size;
      pack(p.second, prev_size, c_buff + offset);
      offset += prev_size;
    }

    return buff;
  }


  template <typename K, typename V>
  static localIndex unpack(map<K, V> & m, const void * buffer, localIndex byte_size=-1)
  {
    const localIndex * buff = reinterpret_cast<const localIndex *>(buffer);
    localIndex bytes_recorded = buff[0];
    localIndex num_pairs = buff[1];

    if (bytes_recorded != byte_size && byte_size >= 0)
    {
      GEOS_ERROR("Number of bytes recorded not equal to number of bytes in buffer: " <<
                  bytes_recorded << " " << byte_size);
    }

    localIndex offset = 0;
    const char * c_buff = reinterpret_cast<const char *>(buff + 2);
    for (localIndex i = 0; i < num_pairs; ++i)
    {
      K key;
      offset += unpack(key, c_buff + offset);
      V value;
      offset += unpack(value, c_buff + offset);
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


  template <typename T>
  static localIndex packed_size(const std::unique_ptr<T> & s)
  {
    GEOS_ERROR("You shouldn't be packing a unique pointer!"); 
    return 0;
  }


  template <typename T>
  static void * pack(const std::unique_ptr<T> & ptr, localIndex & byte_size, void * buffer=nullptr)
  { 
    GEOS_ERROR("You shouldn't be packing a unique pointer!"); 
    byte_size = 0;
    return nullptr;
  }


  template <typename T>
  static localIndex unpack(std::unique_ptr<T> & ptr, const void * buffer, localIndex byte_size=-1)
  { 
    GEOS_ERROR("You shouldn't be unpacking a unique pointer!");
    return 0;
  }


  template <typename T>
  static typename std::enable_if<std::is_same<T, BasisBase>::value ||
                                 std::is_same<T, QuadratureBase>::value ||
                                 std::is_same<T, SimpleGeometricObjectBase>::value ||
                                 std::is_same<T, PartitionBase>::value ||
                                 std::is_same<T, NeighborCommunicator>::value , localIndex>::type
  packed_size(const T & data)
  {
    GEOS_ERROR("You shouldn't be packing a BasisBase, QuadratureBase, SimpleGeometricObjectBase, or PartitionBase!"); 
    return 0;
  }


  template <typename T>
  static typename std::enable_if<std::is_same<T, BasisBase>::value ||
                                 std::is_same<T, QuadratureBase>::value ||
                                 std::is_same<T, SimpleGeometricObjectBase>::value ||
                                 std::is_same<T, PartitionBase>::value ||
                                 std::is_same<T, NeighborCommunicator>::value , void *>::type
  pack(const T & data, localIndex & byte_size, void * buffer=nullptr)
  {
    GEOS_ERROR("You shouldn't be packing a BasisBase, QuadratureBase, SimpleGeometricObjectBase, or PartitionBase!"); 
    byte_size = 0;
    return nullptr;
  }


  template <typename T>
  static typename std::enable_if<std::is_same<T, BasisBase>::value ||
                                 std::is_same<T, QuadratureBase>::value ||
                                 std::is_same<T, SimpleGeometricObjectBase>::value ||
                                 std::is_same<T, PartitionBase>::value ||
                                 std::is_same<T, NeighborCommunicator>::value , localIndex>::type
  unpack(T & data, const void * buffer, localIndex byte_size=-1)
  { 
    GEOS_ERROR("You shouldn't be packing a BasisBase, QuadratureBase, SimpleGeometricObjectBase, or PartitionBase!"); 
    return 0;
  }

private:
  Buffer();
  ~Buffer();

};


}     /* end namespace dataRepository */
}     /* end namespace geosx */


#endif /* DATAREPOSITORY_BUFFERHPP */
