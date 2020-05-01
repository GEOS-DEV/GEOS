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
 * @file WellControls.hpp
 */


#ifndef GEOSX_WELLS_WELLCONTROLS_HPP
#define GEOSX_WELLS_WELLCONTROLS_HPP

#include "dataRepository/Group.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
static constexpr auto wellControls = "WellControls";
}
}


/**
 * @class WellControls
 * @brief This class describes the controls used to operate a well.
 */
class WellControls : public dataRepository::Group
{
public:

  /** Type of wells
   * Either producer or injector.
   */
  enum class Type
  {
    PRODUCER,  /**< A production well */
    INJECTOR   /**< An injection well */
  };

  
  /** Types of well controls
   * Used to specifiy a well's operating conditions
   */
  enum class Control
  {
    BHP,  /**< The well operates at a specified bottom hole pressure (BHP) */
    GASRATE, /**< The well operates at a specified gas flow rate */
    OILRATE, /**< The well operates at a specified oil flow rate */
    WATERRATE, /**< The well operates at a specified water flow rate */
    LIQUIDRATE /**< The well operates at a specified liquid flow rate (oil + water) */
  };

  
  /**
   * @brief Main constructor for WellControls Objects
   * @param [in] name the name of this instantiation of WellControls in the repository
   * @param [in] parent the parent group of this instantiation of WellControls
   */
  explicit WellControls( string const & name, dataRepository::Group * const parent );

  
  /// default destructor
  ~WellControls() override;

  
  /// deleted default constructor
  WellControls() = delete;

  
  /// deleted copy constructor
  WellControls( WellControls const & ) = delete;

  
  /// deleted move constructor
  WellControls( WellControls && ) = delete;

  
  /// deleted assignment operator
  WellControls & operator=( WellControls const & ) = delete;

  
  /// deleted move operator
  WellControls & operator=( WellControls && ) = delete;

  
  /**
   * @brief Set the reference well elem index where the control will be enforced.
   * @param [in] refIndex the reference well elem index where the control will be enforced
   */
  void SetReferenceWellElementIndex( localIndex refIndex )
  {
    m_refWellElemIndex = refIndex;
  }


  /**
   * @brief Get the reference well elem index where the control will be enforced.
   * @return a localIndex value representing the reference well element index where the control will be enforced
   */
  localIndex const & GetReferenceWellElementIndex() const
  {
    return m_refWellElemIndex;
  }

  
  /**
   * @brief Get the well type (injector or producer).
   * @return a well Type enum
   */
  Type GetType() const { return m_type; }

  
  /**
   * @brief Set the control type and numerical value for a well.
   * @param [in] control a Control enum with the type of control that is enforced
   * @param [in] val value for the control (depending on the control type, can be a maximum bottom hole pressure, a minimum water rate...)
   */
  void SetControl( Control control, real64 const & val );

  
  /**
   * @brief Get the control type for the well.
   * @return the Control enum enforced at the well
   */
  Control GetControl() const { return m_currentControl; }

  
  /**
   * @brief Get the target Bottom Hole Pressure value.
   * @return a value for the target Bottom Hole Pressure
   */
  const real64 & GetTargetBHP() const { return m_targetBHP; }

  
  /**
   * @brief Get the target rate
   * @return the target rate
   */
  const real64 & GetTargetRate() const { return m_targetRate; }

  
  /**
   * @brief Get the composition of the injection rate
   * @return a global component fraction vector
   */
  arrayView1d< real64 const > const & GetInjectionStream() const { return m_injectionStream; }

  
  /**
   * @brief Print a console output with some WellControl variables (for console debugging purposes)
   */
  void Debug() const;

  
  struct viewKeyStruct
  {
    static constexpr auto refWellElemIndexString = "referenceWellElementIndex";
    static constexpr auto typeString             = "type";
    static constexpr auto controlString          = "control";
    static constexpr auto targetBHPString        = "targetBHP";
    static constexpr auto targetRateString       = "targetRate";
    static constexpr auto injectionStreamString  = "injectionStream";
    dataRepository::ViewKey referenceIndex  = { refWellElemIndexString };
    dataRepository::ViewKey type            = { typeString };
    dataRepository::ViewKey control         = { controlString };
    dataRepository::ViewKey targetBHP       = { targetBHPString };
    dataRepository::ViewKey targetRate      = { targetRateString };
    dataRepository::ViewKey injectionStream = { injectionStreamString };
  } viewKeysWellControls;

protected:

  virtual void PostProcessInput() override;

  virtual void InitializePostInitialConditions_PreSubGroups( Group * const rootGroup ) override;

private:

  /// well type as string
  string m_typeString;

  /// well type (as Type enum)
  Type m_type;

  /// reference index
  localIndex m_refWellElemIndex;

  /// well controls as string
  string m_inputControlString;

  /// well controls as a Control enum
  Control m_currentControl;
  
  /// Target bottom hole pressure value
  real64 m_targetBHP;

  /// Target rate value
  real64 m_targetRate;

  /// vector with global component fractions at the injector
  array1d< real64 >  m_injectionStream;

};

} //namespace geosx

#endif //GEOSX_MANAGERS_WELLS_WELLCONTROLS_HPP
