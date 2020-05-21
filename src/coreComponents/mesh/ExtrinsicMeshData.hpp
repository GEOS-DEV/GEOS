/*
 * MeshFields.hpp
 *
 *  Created on: May 7, 2020
 *      Author: settgast
 */

#ifndef GEOSX_EXTRINSIC_MESH_DATA_HPP_
#define GEOSX_EXTRINSIC_MESH_DATA_HPP_


#include "managers/ObjectManagerBase.hpp"

namespace geosx
{
/**
 *
 */
namespace extrinsicMeshData
{

/**
 * @struct ParentIndex
 * @brief  Holds the interface for registering and getting access to the
 *         data where the parentIndex will be stored.
 */
struct ParentIndex
{
  static constexpr auto key = "parentIndex";
  using DataType = localIndex;
  using Type = array1d<DataType>;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_2;
  static constexpr auto description = "Index of parent within the mesh object it is registered on.";
};


/**
 * @struct ChildIndex
 * @brief  Holds the interface for registering and getting access to the
 *         data where the childIndex will be stored.
 */
struct ChildIndex
{
  static constexpr auto key = "childIndex";
  using DataType = localIndex;
  using Type = array1d<DataType>;
  static constexpr auto defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_2;
  static constexpr auto description = "Index of child within the  mesh object it is registered on.";
};


struct DegreeFromCrack
{
  static constexpr auto key = "degreeFromCrack";
  using DataType = integer;
  using Type = array1d<DataType>;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1;
  static constexpr auto description = "Connectivity distance from crack.";
};

struct DegreeFromCrackTip
{
  static constexpr auto key = "degreeFromCrackTip";
  using DataType = integer;
  using Type = array1d<integer>;
  static constexpr DataType defaultValue = 100000;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1;
  static constexpr auto description = "Degree of connectivity separation from crack tip.";
};

struct SIFNode
{
  static constexpr auto key = "SIFNode";
  using DataType = real64;
  using Type = array1d<DataType>;
  static constexpr DataType defaultValue = 0;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;
  static constexpr auto description = "SIF on the node.";
};

struct RuptureTime
{
  static constexpr auto key = "ruptureTime";
  using DataType = real64;
  using Type = array1d<real64>;
  static constexpr DataType defaultValue = 1.0e9;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;
  static constexpr auto description = "Time that the node was ruptured.";
};





struct SIF_I
{
  static constexpr auto key = "SIF_I";
  using DataType = real64;
  using Type = array1d<DataType>;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1;
  static constexpr auto description = "SIF_I of the edge.";
};
struct SIF_II
{
  static constexpr auto key = "SIF_II";
  using DataType = real64;
  using Type = array1d<DataType>;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1;
  static constexpr auto description = "SIF_II of the edge.";
};
struct SIF_III
{
  static constexpr auto key = "SIF_III";
  using DataType = real64;
  using Type = array1d<DataType>;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1;
  static constexpr auto description = "SIF_III of the edge.";
};






} // namespace extrinsicMeshData
} // namespace geosx

#define GEOSX_EXTRINSIC_MESH_DATA_HPP_

#endif /* GEOSX_MESH_MESHFIELDS_HPP_ */
