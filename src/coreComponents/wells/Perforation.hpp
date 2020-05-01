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
 * @file Perforation.hpp
 */

#ifndef GEOSX_WELLS_PERFORATION_HPP
#define GEOSX_WELLS_PERFORATION_HPP

#include "dataRepository/Group.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto perforation = "Perforation";
}
}


/**
 * @class Perforation
 *
 * This class describes a perforation with its location, well Peaceman index  and corresponding well element
 */
class Perforation : public dataRepository::Group
{
public:

  /**
   * @brief main constructor for Perforation Objects
   * @param [in] name the name of this instantiation of Perforation in the repository
   * @param [in] parent the parent group of this instantiation of Perforation
   */
  explicit Perforation( string const & name, dataRepository::Group * const parent );

  /// default destructor
  ~Perforation() override;

  
  /// deleted default constructor
  Perforation() = delete;

  
  /// deleted copy constructor
  Perforation( Perforation const & ) = delete;

  
  /// deleted move constructor
  Perforation( Perforation && ) = delete;

  
  /// deleted assignment operator
  Perforation & operator=( Perforation const & ) = delete;

  
  /// deleted move operator
  Perforation & operator=( Perforation && ) = delete;

  
  /**
   * @brief Get the linear distance between the well head and the perforation
   * @return distance between the well head and the perforation
   */
  real64 const & GetDistanceFromWellHead() const { return m_distanceFromHead; }

  
  /**
   * @brief Get the well Peaceman index at the perforation
   * @return the well Peaceman index
   */
  real64 GetWellTransmissibility() const { return m_wellTransmissibility; }


  struct viewKeyStruct
  {
    static constexpr auto distanceFromHeadString  = "distanceFromHead";
    static constexpr auto wellTransmissibilityString = "transmissibility";
    dataRepository::ViewKey distanceFromHead  = { distanceFromHeadString };
    dataRepository::ViewKey wellTransmissibility = { wellTransmissibilityString };
  } viewKeysPerforation;

protected:

  void PostProcessInput() override;

private:

  /// linear distance from well head
  real64 m_distanceFromHead;

  /// well index at this perforation
  real64 m_wellTransmissibility;

};

} //namespace geosx

#endif //GEOSX_WELLS_PERFORATION_HPP
