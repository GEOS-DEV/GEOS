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
 * @file MaterialPropertiesUpdate.hpp
 */

#ifndef SRC_CORECOMPONENTS_FIELDSPECIFICATION_MATERIALPROPERTIESUPDATE_HPP_
#define SRC_CORECOMPONENTS_FIELDSPECIFICATION_MATERIALPROPERTIESUPDATE_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geos
{

/**
 * @class MaterialPropertiesUpdate
 *
 *
 */
class MaterialPropertiesUpdate : public TaskBase
{
public:

  /**
   * @brief Constructor for the state reset class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  MaterialPropertiesUpdate( const string & name,
                            Group * const parent );

  /// Destructor for the class
  ~MaterialPropertiesUpdate() override;

  /// Accessor for the catalog name
  static string catalogName() { return "MaterialPropertiesUpdate"; }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**@}*/

private:

  void postProcessInput() override;
};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_FIELDSPECIFICATION_MATERIALPROPERTIESUPDATE_HPP_ */
