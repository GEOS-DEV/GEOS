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
 * @file PerforationData.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATIONDATA_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATIONDATA_HPP

#include "dataRepository/ManagedGroup.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class MeshLevel;

namespace dataRepository
{
namespace keys
{
static constexpr auto perforations = "Perforations";
}
}

class DomainPartition;
class Perforation;

class PerforationData : public ObjectManagerBase
{
public:

  explicit PerforationData( string const & name, dataRepository::ManagedGroup * const parent );
  ~PerforationData() override;

  PerforationData() = delete;
  PerforationData( PerforationData const &) = delete;
  PerforationData( PerforationData && ) = delete;

  dataRepository::ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  virtual const string getCatalogName() const override;

  localIndex numPerforationsGlobal() const { return integer_conversion<localIndex>(m_perforationList.size()); }
  localIndex numPerforationsLocal()  const { return integer_conversion<localIndex>(size());         }

  Perforation const * getPerforation( localIndex iperf ) const;
  Perforation *       getPerforation( localIndex iperf );

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto reservoirElementRegionString    = "reservoirElementRegion";
    static constexpr auto reservoirElementSubregionString = "reservoirElementSubregion";
    static constexpr auto reservoirElementIndexString     = "reservoirElementIndex";
    static constexpr auto perforationIndexString          = "perforationIndex";

    static constexpr auto gravityDepthString               = "gravityDepth";

    dataRepository::ViewKey reservoirElementRegion    = { reservoirElementRegionString    };
    dataRepository::ViewKey reservoirElementSubregion = { reservoirElementSubregionString };
    dataRepository::ViewKey reservoirElementIndex     = { reservoirElementIndexString     };
    dataRepository::ViewKey perforationIndex          = { perforationIndexString };

    dataRepository::ViewKey gravityDepth              = { gravityDepthString };

  } viewKeysPerforationManager;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto perforationString = "Perforation";

    dataRepository::GroupKey perforation = { perforationString };

  } groupKeysPerforationManager;

protected:

  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override;

private:

  void ConnectToCells( MeshLevel const * domain );
  void PrecomputeData( MeshLevel const * domain );

  array1d<localIndex> m_reservoirElementRegion;
  array1d<localIndex> m_reservoirElementSubregion;
  array1d<localIndex> m_reservoirElementIndex;
  array1d<localIndex> m_perforationIndex;

  array1d<real64> m_gravityDepth;

  string_array m_perforationList;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATIONDATA_HPP
