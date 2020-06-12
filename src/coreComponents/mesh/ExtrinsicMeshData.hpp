/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#ifndef GEOSX_EXTRINSIC_MESH_DATA_HPP_
#define GEOSX_EXTRINSIC_MESH_DATA_HPP_

/**
 * @file ExtrinsicMeshData.hpp
 */

#include "managers/ObjectManagerBase.hpp"

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{


/**
 * Trait struct for registering/accessing the ParentIndex.
 */
struct ParentIndex
{
  static constexpr auto key = "parentIndex";                            ///< The key for registration.
  using DataType = localIndex;                                          ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = -1;                          ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_2; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = "Index of parent within the mesh object it is registered on.";
};


/**
 * Trait struct for registering/accessing the  ChildIndex
 */
struct ChildIndex
{
  static constexpr auto key = "childIndex";                             ///< The key for registration.
  using DataType = localIndex;                                          ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = -1;                          ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_2; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = "Index of child within the mesh object it is registered on.";
};


/**
 * Trait struct for registering/accessing the  DegreeFromCrack
 */
struct DegreeFromCrack
{
  static constexpr auto key = "degreeFromCrack";                        ///< The key for registration.
  using DataType = integer;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = -1;                          ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto
    description = "Distance to the crack in terms of topological distance. "
                  "(i.e. how many nodes are along the path to the closest "
                  "node that is on the crack surface.";
};

/// Trait struct for registering/accessing the  DegreeFromCrackTip
struct DegreeFromCrackTip
{
  static constexpr auto key = "degreeFromCrackTip";                     ///< The key for registration.
  using DataType = integer;                                             ///< The base type for registration of a templated container.
  using Type = array1d< integer >;                                      ///< The type to be registered
  static constexpr DataType defaultValue = 100000;                      ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto
    description = "Distance to the crack tip in terms of topological distance. "
                  "(i.e. how many nodes are along the path to the closest "
                  "node that is on the crack surface.";
};

/// Trait struct for registering/accessing the  SIFNode
struct SIFNode
{
  static constexpr auto key = "SIFNode";                                ///< The key for registration.
  using DataType = real64;                                              ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = 0;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto
    description = "Calculated Stress Intensity Factor on the node.";
};

/// Trait struct for registering/accessing the  RuptureTime
struct RuptureTime
{
  static constexpr auto key = "ruptureTime";                            ///< The key for registration.
  using DataType = real64;                                              ///< The base type for registration of a templated container.
  using Type = array1d< real64 >;                                       ///< The type to be registered
  static constexpr DataType defaultValue = 1.0e9;                       ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = "Time that the object was ruptured/split.";
};

/// Trait struct for registering/accessing the  RuptureRate
struct RuptureRate
{
  static constexpr auto key = "ruptureRate";                            ///< The key for registration.
  using DataType = real64;                                              ///< The base type for registration of a templated container.
  using Type = array1d< real64 >;                                       ///< The type to be registered
  static constexpr DataType defaultValue = 1.0e99;                      ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto
    description = "Rate of rupture in terms of number of objects split per time.";
};


/// Trait struct for registering/accessing the  SIF_I
struct SIF_I
{
  static constexpr auto key = "SIF_I";                                  ///< The key for registration.
  using DataType = real64;                                              ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = -1;                          ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto
    description = "Calculated mode 1 Stress Intensity Factor on the node.";
};

/// Trait struct for registering/accessing the  SIF_II
struct SIF_II
{
  static constexpr auto key = "SIF_II";                                 ///< The key for registration.
  using DataType = real64;                                              ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = -1;                          ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto
    description = "Calculated mode 2 Stress Intensity Factor on the node.";

};

/// Trait struct for registering/accessing the  SIF_III
struct SIF_III
{
  static constexpr auto key = "SIF_III";                                ///< The key for registration.
  using DataType = real64;                                              ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = -1;                          ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_1; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto
    description = "Calculated mode 3 Stress Intensity Factor on the node.";
};

/// Trait struct for registering/accessing the  RuptureState
struct RuptureState
{
  static constexpr auto key = "ruptureState";                           ///< The key for registration.
  using DataType = integer;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = 0;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto
    description = "Rupture state of the face: \n 0=not ready for rupture \n 1=ready for rupture \n 2=ruptured.";
};


/// Trait struct for registering/accessing the  SIFonFace
struct SIFonFace
{
  static constexpr auto key = "SIFonFace";                              ///< The key for registration.
  using DataType = real64;                                              ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = 1;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = "Calculated Stress Intensity Factor on the face.";
};


/// Trait struct for registering/accessing the  K_IC
struct K_IC
{
  static constexpr auto key = "K_IC";                                   ///< The key for registration.
  using DataType = R1Tensor;                                            ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr real64 defaultValue = 1e99;                          ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto
    description = "Critical Stress Intensity Factor :math:`K_{IC}` in the plane of the face.";
};

/// Trait struct for registering/accessing the  K_IC_00
struct K_IC_00
{
  static constexpr auto key = "K_IC_00";                               ///< The key for registration.
  using DataType = real64;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                    ///< The type to be registered
  static constexpr real64 defaultValue = -1;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = ":math:`K_{IC}` on 0-plane, in 0-direction.";
};

/// Trait struct for registering/accessing the  K_IC_01
struct K_IC_01
{
  static constexpr auto key = "K_IC_01";                               ///< The key for registration.
  using DataType = real64;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                    ///< The type to be registered
  static constexpr real64 defaultValue = -1;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = ":math:`K_{IC}` on 0-plane, in 1-direction.";
};

/// Trait struct for registering/accessing the  SIF_II
struct K_IC_02
{
  static constexpr auto key = "K_IC_02";                               ///< The key for registration.
  using DataType = real64;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                    ///< The type to be registered
  static constexpr real64 defaultValue = -1;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = ":math:`K_{IC}` on 0-plane, in 2-direction.";
};

/// Trait struct for registering/accessing the  K_IC_10
struct K_IC_10
{
  static constexpr auto key = "K_IC_10";                               ///< The key for registration.
  using DataType = real64;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                    ///< The type to be registered
  static constexpr real64 defaultValue = -1;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = ":math:`K_{IC}` on 1-plane, in 0-direction.";
};

/// Trait struct for registering/accessing the  K_IC_11
struct K_IC_11
{
  static constexpr auto key = "K_IC_11";                               ///< The key for registration.
  using DataType = real64;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                    ///< The type to be registered
  static constexpr real64 defaultValue = -1;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = ":math:`K_{IC}` on 1-plane, in 1-direction.";
};

/// Trait struct for registering/accessing the  K_IC_12
struct K_IC_12
{
  static constexpr auto key = "K_IC_12";                               ///< The key for registration.
  using DataType = real64;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                    ///< The type to be registered
  static constexpr real64 defaultValue = -1;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = ":math:`K_{IC}` on 1-plane, in 2-direction.";
};

/// Trait struct for registering/accessing the  K_IC_20
struct K_IC_20
{
  static constexpr auto key = "K_IC_20";                               ///< The key for registration.
  using DataType = real64;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                    ///< The type to be registered
  static constexpr real64 defaultValue = -1;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = ":math:`K_{IC}` on 2-plane, in 0-direction.";
};

/// Trait struct for registering/accessing the  K_IC_21
struct K_IC_21
{
  static constexpr auto key = "K_IC_21";                               ///< The key for registration.
  using DataType = real64;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                    ///< The type to be registered
  static constexpr real64 defaultValue = -1;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = ":math:`K_{IC}` on 2-plane, in 1-direction.";
};

/// Trait struct for registering/accessing the  K_IC_22
struct K_IC_22
{
  static constexpr auto key = "K_IC_22";                               ///< The key for registration.
  using DataType = real64;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                    ///< The type to be registered
  static constexpr real64 defaultValue = -1;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::NOPLOT; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = ":math:`K_{IC}` on 2-plane, in 2-direction.";
};

/// Trait struct for registering/accessing the  PrimaryCandidateFace
struct PrimaryCandidateFace
{
  static constexpr auto key = "primaryCandidateFace";                   ///< The key for registration.
  using DataType = localIndex;                                          ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = 0;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = "??.";
};

/// Trait struct for registering/accessing the  IsFaceSeparable
struct IsFaceSeparable
{
  static constexpr auto key = "isFaceSeparable";                        ///< The key for registration.
  using DataType = integer;                                             ///< The base type for registration of a templated container.
  using Type = array1d< DataType >;                                     ///< The type to be registered
  static constexpr DataType defaultValue = 0;                           ///< The dataRepository::DefaultValue
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0; ///< The default dataRepository::PlotLevel

  /// Description of the data associated with this trait.
  static constexpr auto description = "A flag to mark if the face is separable.";
};


} // namespace extrinsicMeshData
} // namespace geosx

#define GEOSX_EXTRINSIC_MESH_DATA_HPP_

#endif /* GEOSX_MESH_MESHFIELDS_HPP_ */
