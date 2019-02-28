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

#ifndef GEOSX_CORECOMPONENTS_WELLS_PERFORATIONDATA_HPP
#define GEOSX_CORECOMPONENTS_WELLS_PERFORATIONDATA_HPP

#include "dataRepository/ManagedGroup.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto perforationData = "PerforationData";
}
}

class DomainPartition;
class MeshLevel;
class Perforation;

class PerforationData : public ObjectManagerBase
{
public:

  explicit PerforationData( string const & name, dataRepository::ManagedGroup * const parent );
  ~PerforationData() override;

  PerforationData() = delete;
  PerforationData( PerforationData const &) = delete;
  PerforationData( PerforationData && ) = delete;

  virtual const string getCatalogName() const override;
  
  dataRepository::ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  localIndex numPerforationsLocal()  const
  { return integer_conversion<localIndex>(size()); }
  localIndex numPerforationsGlobal() const;
  
  Perforation const * getPerforation( localIndex iperf ) const;
  Perforation *       getPerforation( localIndex iperf );

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto reservoirElementRegionString    = "reservoirElementRegion";
    static constexpr auto reservoirElementSubregionString = "reservoirElementSubregion";
    static constexpr auto reservoirElementIndexString     = "reservoirElementIndex";
    static constexpr auto wellElementIndexString          = "wellElementIndex";
    static constexpr auto perforationIndexString          = "perforationIndex";
    static constexpr auto transmissibilityString          = "transmissibility";
    static constexpr auto gravityDepthString              = "gravityDepth";

    dataRepository::ViewKey reservoirElementRegion    = { reservoirElementRegionString };
    dataRepository::ViewKey reservoirElementSubregion = { reservoirElementSubregionString };
    dataRepository::ViewKey reservoirElementIndex     = { reservoirElementIndexString };
    dataRepository::ViewKey wellElementIndex          = { wellElementIndexString };
    dataRepository::ViewKey perforationIndex          = { perforationIndexString };
    dataRepository::ViewKey transmissibility          = { transmissibilityString };
    dataRepository::ViewKey gravityDepth              = { gravityDepthString };

  } viewKeysPerforationData;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
  } groupKeysPerforationData;

protected:

  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override;

private:

  void ConnectToCells( MeshLevel const * domain );

  array1d<localIndex> m_reservoirElementRegion;
  array1d<localIndex> m_reservoirElementSubregion;
  array1d<localIndex> m_reservoirElementIndex;
  array1d<localIndex> m_wellElementIndex;
  array1d<localIndex> m_perforationIndex;
  array1d<real64> m_transmissibility;
  array1d<real64> m_gravityDepth;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_WELLS_PERFORATIONDATA_HPP
