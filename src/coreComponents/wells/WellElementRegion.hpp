/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file WellElementRegion.hpp
 *
 */

#ifndef CORECOMPONENTS_WELLS_WELLELEMENTREGION_HPP_
#define CORECOMPONENTS_WELLS_WELLELEMENTREGION_HPP_

#include "mesh/ElementRegionBase.hpp"
#include "InternalWellGenerator.hpp"

namespace geosx
{

class MeshLevel;

/**
 * @class WellElementRegion
 *
 */
class WellElementRegion : public ElementRegionBase
{
public:

  /**
   * @brief constructor
   * @param name The name of the object in the data hierarchy.
   * @param parent Pointer to the parent group in the data hierarchy.
   */
  WellElementRegion( string const & name, ManagedGroup * const parent );

  WellElementRegion() = delete;
  virtual ~WellElementRegion() override;

  /**
   * @brief The key name for the WellElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string CatalogName()
  { return "WellElementRegion"; }

  /**
   *
   * @return the name of this type in the catalog
   */
  virtual const string getCatalogName() const override final
  { return WellElementRegion::CatalogName(); }

  /**
   * @brief Setter for the name of the InternalWellGenerator object of this well
   * @param name the name of the InternalWellGenerator object
   */
  void SetWellGeneratorName( string const & name ) { m_wellGeneratorName = name; }

  /**
   * @brief Getter for the name of the InternalWellGenerator object of this well
   * @param the name of the InternalWellGenerator object
   */
  string const & GetWellGeneratorName() const { return m_wellGeneratorName; }

  /**
   * @brief Setter for the name of the WellControls object of this well
   * @param name the name of the WellControls object
   */
  void SetWellControlsName( string const & name ) { m_wellControlsName = name; }

  /**
   * @brief Getter for the name of the subRegion
   * @param the name of the subRegion object
   */
  string const & GetSubRegionName() const { return m_subRegionName; }

  
  // not implemented, this task is performed in GenerateWell
  virtual void GenerateMesh( ManagedGroup const * ) override {}

 /**
  * @brief Build the local well elements and perforations from global well geometry
  * @param[in] meshLevel the mesh object (single level only)
  * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
  * @param[in] nodeOffsetGlobal the offset of the first global well node ( = offset of last global mesh node + 1 )
  * @param[in] elemOffsetGlobal the offset of the first global well element ( = offset of last global mesh elem + 1 )
  */
  void GenerateWell( MeshLevel & mesh, 
                     InternalWellGenerator const & wellGeometry,
                     globalIndex nodeOffsetGlobal,
                     globalIndex elemOffsetGlobal );

  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    static constexpr auto wellControlsString  = "wellControlsName";
    static constexpr auto wellGeneratorString = "wellGeneratorName";

    dataRepository::ViewKey wellControlsName  = { wellControlsString };
    dataRepository::ViewKey wellGeneratorName = { wellGeneratorString };

  } viewKeysWellElementRegion;

  struct groupKeyStruct : public ElementRegionBase::groupKeyStruct
  {
  } groupKeysWellElementRegion;

private:

  /// name of the (unique) subregion
  const string m_subRegionName;

  /// name of the well constraint 
  string m_wellControlsName;
 
  /// name of the well generator
  string m_wellGeneratorName;

};

} /* namespace geosx */

#endif /* CORECOMPONENTS_WELLS_WELLELEMENTREGION_HPP_ */
