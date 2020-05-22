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
 * @name Trait definitions for the surface generators.
 */
///@{

/**
 * @struct ParentIndex
 * @brief  Holds the interface for registering and getting access to the
 *         data where the parentIndex will be stored.
 */
struct ParentIndex
{
  static constexpr auto key = "parentIndex";
  using DataType = localIndex;
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_2;
  static constexpr auto description = "Index of parent within the mesh object it is registered on.";
  dataRepository::ViewKey viewKey = { key };
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
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_2;
  static constexpr auto description = "Index of child within the  mesh object it is registered on.";
  dataRepository::ViewKey viewKey = { key };
};


struct DegreeFromCrack
{
  static constexpr auto key = "degreeFromCrack";
  using DataType = integer;
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1;
  static constexpr auto description = "Connectivity distance from crack.";
  dataRepository::ViewKey viewKey = { key };
};

struct DegreeFromCrackTip
{
  static constexpr auto key = "degreeFromCrackTip";
  using DataType = integer;
  using Type = array1d< integer >;
  static constexpr DataType defaultValue = 100000;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1;
  static constexpr auto description = "Degree of connectivity separation from crack tip.";
  dataRepository::ViewKey viewKey = { key };
};

struct SIFNode
{
  static constexpr auto key = "SIFNode";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = 0;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;
  static constexpr auto description = "SIF on the node.";
  dataRepository::ViewKey viewKey = { key };
};

struct RuptureTime
{
  static constexpr auto key = "ruptureTime";
  using DataType = real64;
  using Type = array1d< real64 >;
  static constexpr DataType defaultValue = 1.0e9;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;
  static constexpr auto description = "Time that the object was ruptured.";
  dataRepository::ViewKey viewKey = { key };
};


struct RuptureRate
{
  static constexpr auto key = "ruptureRate";
  using DataType = real64;
  using Type = array1d< real64 >;
  static constexpr DataType defaultValue = 1.0e99;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;
  static constexpr auto description = "Rate of rupture.";
  dataRepository::ViewKey viewKey = { key };
};



struct SIF_I
{
  static constexpr auto key = "SIF_I";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1;
  static constexpr auto description = "SIF_I of the edge.";
  dataRepository::ViewKey viewKey = { key };
};
struct SIF_II
{
  static constexpr auto key = "SIF_II";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1;
  static constexpr auto description = "SIF_II of the edge.";
  dataRepository::ViewKey viewKey = { key };
};
struct SIF_III
{
  static constexpr auto key = "SIF_III";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1;
  static constexpr auto description = "SIF_III of the edge.";
  dataRepository::ViewKey viewKey = { key };
};

struct RuptureState
{
  static constexpr auto key = "ruptureState";
  using DataType = integer;
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = 0;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;
  static constexpr auto description = "Rupture state of the face.0=not ready for rupture. 1=ready for rupture. 2=ruptured.";
  dataRepository::ViewKey viewKey = { key };
};


struct SIFonFace
{
  static constexpr auto key = "SIFonFace";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = 1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;
  static constexpr auto description = "SIF on the face.";
  dataRepository::ViewKey viewKey = { key };
};


struct K_IC
{
  static constexpr auto key = "K_IC";
  using DataType = R1Tensor;
  using Type = array1d< DataType >;
  static constexpr real64 defaultValue = 1e99;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;
  static constexpr auto description = "K_IC in each plane.";
  dataRepository::ViewKey viewKey = { key };
};



struct K_IC_00
{
  static constexpr auto key = "K_IC_00";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr real64 defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT;
  static constexpr auto description = "K_IC on 0-plane, in 0-direction.";
  dataRepository::ViewKey viewKey = { key };
};
struct K_IC_01
{
  static constexpr auto key = "K_IC_01";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr real64 defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT;
  static constexpr auto description = "K_IC on 0-plane, in 1-direction.";
  dataRepository::ViewKey viewKey = { key };
};
struct K_IC_02
{
  static constexpr auto key = "K_IC_02";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr real64 defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT;
  static constexpr auto description = "K_IC on 0-plane, in 2-direction.";
  dataRepository::ViewKey viewKey = { key };
};
struct K_IC_10
{
  static constexpr auto key = "K_IC_10";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr real64 defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT;
  static constexpr auto description = "K_IC on 1-plane, in 0-direction.";
  dataRepository::ViewKey viewKey = { key };
};
struct K_IC_11
{
  static constexpr auto key = "K_IC_11";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr real64 defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT;
  static constexpr auto description = "K_IC on 1-plane, in 1-direction.";
  dataRepository::ViewKey viewKey = { key };
};
struct K_IC_12
{
  static constexpr auto key = "K_IC_12";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr real64 defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT;
  static constexpr auto description = "K_IC on 1-plane, in 2-direction.";
  dataRepository::ViewKey viewKey = { key };
};
struct K_IC_20
{
  static constexpr auto key = "K_IC_20";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr real64 defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT;
  static constexpr auto description = "K_IC on 2-plane, in 0-direction.";
  dataRepository::ViewKey viewKey = { key };
};
struct K_IC_21
{
  static constexpr auto key = "K_IC_21";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr real64 defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT;
  static constexpr auto description = "K_IC on 2-plane, in 1-direction.";
  dataRepository::ViewKey viewKey = { key };
};
struct K_IC_22
{
  static constexpr auto key = "K_IC_22";
  using DataType = real64;
  using Type = array1d< DataType >;
  static constexpr real64 defaultValue = -1;
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT;
  static constexpr auto description = "K_IC on 2-plane, in 2-direction.";
  dataRepository::ViewKey viewKey = { key };
};



struct PrimaryCandidateFace
{
  static constexpr auto key = "primaryCandidateFace";
  using DataType = localIndex;
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = 0;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;
  static constexpr auto description = "SIF_III of the edge.";
  dataRepository::ViewKey viewKey = { key };
};

struct IsFaceSeparable
{
  static constexpr auto key = "isFaceSeparable";
  using DataType = integer;
  using Type = array1d< DataType >;
  static constexpr DataType defaultValue = 0;
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;
  static constexpr auto description = "A flag to mark if the face is separable.";
  dataRepository::ViewKey viewKey = { key };
};


///@}


/**
 * @name Trait definitions for the Solid Mechanics Solvers.
 */
///@{

///@}

/**
 * @name Trait definitions for the Flow Solvers.
 */
///@{

///@}



} // namespace extrinsicMeshData
} // namespace geosx

#define GEOSX_EXTRINSIC_MESH_DATA_HPP_

#endif /* GEOSX_MESH_MESHFIELDS_HPP_ */
