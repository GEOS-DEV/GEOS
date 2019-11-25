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
 *
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
 *
 * This class describes the controls used to operate a well
 */  
class WellControls : public dataRepository::Group
{
public:

  // define the type of well (producer or injector)
  enum class Type
  {
    PRODUCER,
    INJECTOR
  };

  // define the well control
  enum class Control
  {
    BHP,
    GASRATE,
    OILRATE,
    WATERRATE,
    LIQUIDRATE
  };

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  explicit WellControls( string const & name, dataRepository::Group * const parent );
  
  /**
   * @brief default destructor
   */
  ~WellControls() override;

  /// deleted default constructor
  WellControls() = delete;

  /// deleted copy constructor
  WellControls( WellControls const &) = delete;

  /// deleted move constructor
  WellControls( WellControls && ) = delete;

  /// deleted assignment operator
  WellControls & operator=( WellControls const & ) = delete;

  /// deleted move operator
  WellControls & operator=( WellControls && ) = delete;

  /**
   * @brief Setter for the reference well elem index where the control will be enforced
   * @param refIndex the reference well elem index where the control will be enforced
   */
  void SetReferenceWellElementIndex( localIndex refIndex )
  {
    m_refWellElemIndex = refIndex;
  }


  /**
   * @brief Getter for the reference well elem index where the control will be enforced
   * @return the reference well element index where the control will be enforced
   */
  localIndex const & GetReferenceWellElementIndex() const
  {
    return m_refWellElemIndex;
  }

  /**
   * @brief Getter for the well type (injector or producer)
   * @return the well type
   */
  Type GetType() const { return m_type; }

  /**
   * @brief Setter for the control equation at this well
   * @param control the type of control that is enforced
   * @param val the value of the control (max pressure, min rate, etc)
   */
  void SetControl( Control control, real64 const & val );

  /**
   * @brief Getter for the control equation at this well
   * @return the type of control that is enforced at this well
   */
  Control GetControl() const { return m_currentControl; }

  /**
   * @brief Getter for the target BH pressure
   * @return the target BH pressure
   */
  const real64 & GetTargetBHP() const { return m_targetBHP; }

  /**
   * @brief Getter for the target rate
   * @return the target rate
   */
  const real64 & GetTargetRate() const { return m_targetRate; }

  /**
   * @brief Getter for the composition of the injection rate
   * @param the index of the component
   * @return the global component fraction for component ic
   */
  real64 GetInjectionStream( localIndex ic ) const;

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

  /// well type
  string m_typeString;
  Type m_type;

  /// reference index
  localIndex m_refWellElemIndex;

  /// well controls
  string  m_inputControlString;
  Control m_currentControl;
  real64  m_targetBHP;
  real64  m_targetRate;
 
  /// global component fraction at the injector
  array1d<real64>  m_injectionStream;

};

} //namespace geosx

#endif //GEOSX_MANAGERS_WELLS_WELLCONTROLS_HPP
