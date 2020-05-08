/*
 * MeshFields.hpp
 *
 *  Created on: May 7, 2020
 *      Author: settgast
 */

#ifndef GEOSX_MESH_MESHFIELDS_HPP_
#define GEOSX_MESH_MESHFIELDS_HPP_


#include "managers/ObjectManagerBase.hpp"

namespace geosx
{
/**
 *
 */
namespace extrinsicMeshData
{

/**
 * @struct ExtrinisicMeshDataBase
 * @brief  A base class for use in a CRT pattern that provides the getters
 *         to all the derived types.
 * @tparam LEAF The Derived Type that uses #ExtrinisicMeshDataBase as a base.
 *
 * More about the crt pattern
 */
template< typename LEAF >
struct ExtrinisicMeshDataBase
{


  static void registerField( ObjectManagerBase & om,
                             string const & nameOfRegisteringObject )
  {
    om.registerWrapper< typename LEAF::Type >( LEAF::key )->
        setApplyDefaultValue( LEAF::defaultValue )->
        setPlotLevel( LEAF::plotLevel )->
        setDescription( LEAF::description )->
        setRegisteringObjects( nameOfRegisteringObject );
  }


  /**
   *
   * @param om
   * @return
   */
  static auto const & get( ObjectManagerBase const & om )
  {
    return om.getReference<typename LEAF::Type>( LEAF::key).toViewConst();
  }

  /**
   *
   * @param om
   * @return
   */
  static auto const & get( ObjectManagerBase & om )
  {
    return om.getReference<typename LEAF::Type>(LEAF::key).toView();
  }
};




/**
 * @struct ParentIndex
 * @brief  Holds the interface for registering and getting access to the
 *         data where the parentIndex will be stored.
 */
struct ParentIndex : public ExtrinisicMeshDataBase< ParentIndex >
{
  /// char[] key for the parentIndex.
  static constexpr auto key = "parentIndex";

  /// Type of data to be registered with the repository.
  using Type = array1d<localIndex>;

  /// The default value for the data.
  static constexpr auto defaultValue = -1;

  /// The plot level for the registered data.
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_2;

  /// The description to place in the sphinx documentation for this entry.
  static constexpr auto description = "Index of parent within the mesh object it is registered on.";
};


/**
 * @struct ChildIndex
 * @brief  Holds the interface for registering and getting access to the
 *         data where the childIndex will be stored.
 */
struct ChildIndex : public ExtrinisicMeshDataBase< ChildIndex >
{
  /// char[] key for the childIndex.
  static constexpr auto key = "childIndex";

  /// Type of data to be registered with the repository.
  using Type = array1d<localIndex>;

  /// The default value for the data.
  static constexpr auto defaultValue = -1;

  /// The plot level for the registered data.
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_2;

  /// The description to place in the sphinx documentation for this entry.
  static constexpr auto description = "Index of child within the  mesh object it is registered on.";

};


} // namespace meshfields
} // namespace geosx

#define USE_EXTRINSIC_MESH_DATA

#endif /* GEOSX_MESH_MESHFIELDS_HPP_ */
