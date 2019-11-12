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

  InternalWellGenerator( const std::string& name,
                         Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~InternalWellGenerator() override;

  /**
   * @return the name of this type in the catalog
   */  
  static string CatalogName() { return "InternalWell"; }

  struct viewKeyStruct : WellGeneratorBase::viewKeyStruct
  {
    constexpr static auto nodeCoords       = "polylineNodeCoords";
    constexpr static auto segmentConn      = "polylineSegmentConn";
  };

  /// not implemented
  virtual void GenerateElementRegions( DomainPartition & GEOSX_UNUSED_ARG( domain ) ) override {}

  /// not implemented 
  virtual void GetElemToNodesRelationInBox ( std::string const & GEOSX_UNUSED_ARG( elementType ),
                                             int const * GEOSX_UNUSED_ARG( index ),
                                             int const & GEOSX_UNUSED_ARG( iEle ),
                                             int * GEOSX_UNUSED_ARG( nodeIDInBox ),
                                             int const GEOSX_UNUSED_ARG( size )) override {}

  /// not implemented
  virtual void RemapMesh ( dataRepository::Group * const GEOSX_UNUSED_ARG( domain ) ) override {}

protected:
  void PostProcessInput() override final;


private:
  void GeneratePolyLine() override final;

};
}

#endif /* GEOSX_WELLS_INTERNALWELLGENERATOR_HPP_ */
