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

/*
 * @file InternalWellGenerator.hpp
 *
 */

#ifndef GEOSX_WELLS_INTERNALWELLGENERATOR_HPP_
#define GEOSX_WELLS_INTERNALWELLGENERATOR_HPP_

#include "dataRepository/Group.hpp"
#include "WellGeneratorBase.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const nodeCoords       = "polylineNodeCoords";
string const segmentConn      = "polylineSegmentConn";
}
}


/**
 * @class InternalWellGenerator
 *
 * This class processes the data of a single well from the XML and generates the well geometry
 */
class InternalWellGenerator : public WellGeneratorBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  InternalWellGenerator( const std::string & name,
                         Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~InternalWellGenerator() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this type in the catalog
   */
  static string CatalogName() { return "InternalWell"; }

  ///@}

  /**
   * @name Overriding functions defined in MeshGeneratorBase and above
   */
  ///@{

  /// not implemented
  virtual void GetElemToNodesRelationInBox ( std::string const &,
                                             int const *,
                                             int const &,
                                             int *,
                                             int const ) override {}

  /// not implemented
  virtual void RemapMesh ( dataRepository::Group * const ) override {}

protected:
  void PostProcessInput() override final;

private:
  void GeneratePolyLine() override final;
};
}

#endif /* GEOSX_WELLS_INTERNALWELLGENERATOR_HPP_ */
