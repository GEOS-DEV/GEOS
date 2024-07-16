/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StubGroup.hpp
 */

#ifndef GEOS_DATAREPOSITORY_STUBGROUP_HPP_
#define GEOS_DATAREPOSITORY_STUBGROUP_HPP_

#include "Group.hpp"
#include "Wrapper.hpp"

namespace geos
{

namespace dataRepository
{

template<typename T>
struct StubWrapperInfo
{
  using type = T;
  string name;
  T defaultValue;
  string description;
  InputFlags inputFlag;
  RestartFlags restartFlag;
  PlotLevel plotLevel;
};

// used to allow terse constructor usage in C++17 due to lack of CTAD for aggregates (available in C++20)
template <typename T>
GEOS_FORCE_INLINE
auto stubWrapper( string const & name,
                  T const & defaultValue,
                  string const & description = "",
                  InputFlags inputFlag = InputFlags::OPTIONAL ,
                  RestartFlags restartFlag = RestartFlags::WRITE,
                  PlotLevel plotLevel = PlotLevel::NOPLOT )
{
  return StubWrapperInfo< T >{ name, description, defaultValue, inputFlag, restartFlag, plotLevel };
}

template < typename Base, typename... StubWrapperInfos >
class StubGroup : public Base
{
public:
  using super = Base;

  StubGroup(std::string const & name, Group * const parent, StubWrapperInfos... ws)
    : Base(name, parent)
  {
    // Register each wrapper using the data in StubWrapperInfos...
    registerWrappersFromInfo(ws...);
  }
private:
  // Helper to register a single wrapper
  template <typename WrapperInfo>
  void registerWrapper(WrapperInfo const & info)
  {
    // we default-initialize a unique_ptr so the wrapper takes ownership of the data
    auto & wrapper = this->template registerWrapper<typename WrapperInfo::type>(info.name, std::unique_ptr< typename WrapperInfo::type >{});
    wrapper.setInputFlag(info.inputFlags)
           .setApplyDefault(info.defaultValue)
           .setDescription(info.description)
           .setRestartFlags(info.restartFlags)
           .setPlotLevel(info.plotLevel);
  }

  // Expand the pack of WrapperInfo and register each one
  template <typename... Ws>
  void registerWrappersFromInfo(Ws const &... ws)
  {
    (registerWrapper(ws), ...);
  }
};

} // namesapce dataRepository

} // namespace geos

#endif