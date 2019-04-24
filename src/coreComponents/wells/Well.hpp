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
 * @file Well.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_WELL_HPP_
#define GEOSX_CORECOMPONENTS_WELLS_WELL_HPP_

#include "managers/ObjectManagerBase.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "WellElementManager.hpp"
#include "PerforationData.hpp"
#include "PerforationManager.hpp"

namespace geosx
{

class Well : public ObjectManagerBase
{
public:

  // define the type of well (producer or injector)
  enum class Type {
                    PRODUCER,
                    INJECTOR
                  };

  // define the well control
  enum class Control {
                       BHP,
                       GASRATE,
                       OILRATE,
                       WATERRATE,
                       LIQUIDRATE
                     };
  
  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */  
  explicit Well( string const & name, 
                 dataRepository::ManagedGroup * const parent );

  /**
   * @brief default destructor
   */
  virtual ~Well() override;

  /// deleted default constructor
  Well() = delete;

  /// deleted copy constructor
  Well( Well const & ) = delete;

  /// deleted move constructor
  Well( Well && ) = delete;

  /// deleted assignment operator
  Well & operator=( Well const & ) = delete;

  /// deleted move operator
  Well & operator=( Well && ) = delete;

  /// Catalog name interface
  static string CatalogName() { return "Well"; }

  using CatalogInterface = cxx_utilities::CatalogInterface< Well, std::string const &, ManagedGroup * const >;

  static CatalogInterface::CatalogType & GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual const string getCatalogName() const override { return CatalogName(); }

  virtual dataRepository::ManagedGroup * CreateChild(string const & childKey, string const & childName) override;

  virtual void InitializePostSubGroups( ManagedGroup * const rootGroup ) override;

  virtual void PostProcessInput() override;

  WellElementSubRegion * getWellElements()             { return &m_wellElementSubRegion; }
  WellElementSubRegion const * getWellElements() const { return &m_wellElementSubRegion; }
  
  PerforationData * getPerforations()             { return &m_perforationData; }
  PerforationData const * getPerforations() const { return &m_perforationData; }
  
  R1Tensor const & getGravityVector() const;

  real64 const & getReferenceDepth() const { return m_referenceDepth; }

  Type getType() const { return m_type; }

  void setControl( Control control, real64 const & val );
  Control getControl() const { return m_currentControl; }

  real64 const & getTargetBHP() const { return m_targetBHP; }

  real64 const & getTargetRate() const { return m_targetRate; }

  real64 const & getInjectionStream( localIndex ic ) const
  { return m_injectionStream[ic]; }
  
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto referenceDepthString  = "referenceDepth";
    static constexpr auto typeString            = "type";
    static constexpr auto controlString         = "control";
    static constexpr auto targetBHPString       = "targetBHP";
    static constexpr auto targetRateString      = "targetRate";
    static constexpr auto injectionStreamString = "injectionStream";
    
    dataRepository::ViewKey referenceDepth  = { referenceDepthString };
    dataRepository::ViewKey type            = { typeString };
    dataRepository::ViewKey control         = { controlString };
    dataRepository::ViewKey targetBHP       = { targetBHPString };
    dataRepository::ViewKey targetRate      = { targetRateString };
    dataRepository::ViewKey injectionStrean = { injectionStreamString };

  } viewKeysWell;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto wellElementsString    = dataRepository::keys::wellElements;
    static constexpr auto wellElementDataString = dataRepository::keys::wellElementData;
    static constexpr auto perforationsString    = dataRepository::keys::perforations;
    static constexpr auto perforationDataString = dataRepository::keys::perforationData;

    dataRepository::GroupKey wellElements    = { wellElementsString };
    dataRepository::GroupKey wellElementData = { wellElementDataString };
    dataRepository::GroupKey perforations    = { perforationsString };
    dataRepository::GroupKey perforationData = { perforationDataString };

  } groupKeysWell;

protected:

  // segments
  WellElementSubRegion m_wellElementSubRegion;
  WellElementManager m_wellElementManager;

  // perforation
  PerforationData    m_perforationData;
  PerforationManager m_perforationManager;

  // depth
  real64 m_referenceDepth;
  
  // well type
  string  m_typeString;
  Type    m_type;

  // well controls
  string  m_controlString;
  Control m_currentControl;
  real64  m_targetBHP;
  real64  m_targetRate;

  // global component fraction at the injector
  array1d<real64>  m_injectionStream;
};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_WELLS_WELL_HPP_
