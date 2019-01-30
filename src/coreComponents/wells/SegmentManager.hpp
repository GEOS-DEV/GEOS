/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * @file SegmentManager.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SEGMENTMANAGER_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SEGMENTMANAGER_HPP

#include "dataRepository/ManagedGroup.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class MeshLevel;

namespace dataRepository
{
namespace keys
{
static constexpr auto segments = "Segments";
}
}

class DomainPartition;
class Segment;

class SegmentManager : public ObjectManagerBase
{
public:

  explicit SegmentManager( string const & name, dataRepository::ManagedGroup * const parent );
  ~SegmentManager() override;

  SegmentManager() = delete;
  SegmentManager( SegmentManager const &) = delete;
  SegmentManager( SegmentManager && ) = delete;

  dataRepository::ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  virtual const string getCatalogName() const override;

  localIndex numSegmentsGlobal() const { return integer_conversion<localIndex>(m_segmentList.size()); }
  localIndex numSegmentsLocal()  const { return integer_conversion<localIndex>(size());         }

  Segment const * getSegment( localIndex iseg ) const;
  Segment *       getSegment( localIndex iseg );

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto segElementRegionString    = "segElementRegion";
    static constexpr auto segElementSubregionString = "segElementSubregion";
    static constexpr auto segElementIndexString     = "segElementIndex";
    static constexpr auto segIndexString            = "segIndex";

    static constexpr auto gravityDepthString        = "gravityDepth";

    using ViewKey = dataRepository::ViewKey;
    
    ViewKey segElementRegion    = { segElementRegionString };
    ViewKey segElementSubregion = { segElementSubregionString };
    ViewKey segElementIndex     = { segElementIndexString };
    ViewKey segIndex            = { segIndexString };

    ViewKey gravityDepth        = { gravityDepthString };

  } viewKeysPerfManager;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto segmentString = "Segment";

    dataRepository::GroupKey segment = { segmentString };

  } groupKeysSegmentManager;

protected:

  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override;

private:

  void PrecomputeData( MeshLevel const * domain );

  array1d<localIndex> m_segElementRegion;
  array1d<localIndex> m_segElementSubregion;
  array1d<localIndex> m_segElementIndex;
  array1d<localIndex> m_segIndex;

  array1d<real64> m_gravityDepth;

  string_array m_segmentList;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SEGMENTMANAGER_HPP
