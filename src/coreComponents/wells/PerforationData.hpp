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

/**
 * @class PerforationData
 *
 * This class keeps tracks of all the local perforations on this rank
 */  
class PerforationData : public ObjectManagerBase
{
public:

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  explicit PerforationData( string const & name, 
                            dataRepository::ManagedGroup * const parent );

  /**
   * @brief default destructor
   */
  ~PerforationData() override;

  /// deleted default constructor
  PerforationData() = delete;

  /// deleted copy constructor
  PerforationData( PerforationData const &) = delete;

  /// deleted move constructor
  PerforationData( PerforationData && ) = delete;

  /// deleted assignment operator
  PerforationData & operator=( PerforationData const & ) = delete;

  /// deleted move operator
  PerforationData & operator=( PerforationData && ) = delete;

  virtual const string getCatalogName() const override;
  
  dataRepository::ManagedGroup * CreateChild( string const & childKey, 
                                              string const & childName ) override;

  /**
   * @brief Getter for the number of perforations on this rank
   * @return the number of perforations on this rank
   */
  localIndex numPerforationsLocal()  const
  { return integer_conversion<localIndex>(size()); }

  /**
   * @brief Getter for the total number of perforations
   * @return the number of perforations on this rank
   */
  localIndex numPerforationsGlobal() const;

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto reservoirElementRegionString    = "reservoirElementRegion";
    static constexpr auto reservoirElementSubregionString = "reservoirElementSubregion";
    static constexpr auto reservoirElementIndexString     = "reservoirElementIndex";

    static constexpr auto wellElementIndexString          = "wellElementIndex";

    static constexpr auto locationString                  = "location";
    static constexpr auto transmissibilityString          = "transmissibility";
    static constexpr auto gravityDepthString              = "gravityDepth";

    dataRepository::ViewKey reservoirElementRegion    = { reservoirElementRegionString };
    dataRepository::ViewKey reservoirElementSubregion = { reservoirElementSubregionString };
    dataRepository::ViewKey reservoirElementIndex     = { reservoirElementIndexString };

    dataRepository::ViewKey wellElementIndex          = { wellElementIndexString };

    dataRepository::ViewKey location                  = { locationString };
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

  // indices of the reseroir well connected to this perforation
  array1d<localIndex> m_reservoirElementRegion;
  array1d<localIndex> m_reservoirElementSubregion;
  array1d<localIndex> m_reservoirElementIndex;

  // indices of the well element to which perforations are attached
  array1d<localIndex> m_wellElementIndex;

  // geometric info on this perforation
  array1d<R1Tensor> m_location;
  array1d<real64> m_transmissibility;
  array1d<real64> m_gravityDepth;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_WELLS_PERFORATIONDATA_HPP
