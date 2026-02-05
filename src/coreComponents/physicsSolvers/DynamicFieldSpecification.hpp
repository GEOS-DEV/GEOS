/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DynamicFieldSpecification.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_DYNAMIC_FIELDSPECIFCATION_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_DYNAMIC_FIELDSPECIFCATION_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geos
{


class DynamicFieldSpecification : public TaskBase
{
public:

  /**
   * @brief Constructor for the initialization class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  DynamicFieldSpecification( const string & name,
                             Group * const parent );

  /// Destructor for the class
  ~DynamicFieldSpecification() override;

  /// Accessor for the catalog name
  static string catalogName()
  {
    return "DynamicFieldSpecification";
  }

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

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    constexpr static char const * fieldSpecificationNamesString() { return "fieldSpecificationNames"; }
  };

  void postInputInitialization() override;

  /**
   * @brief Return the target field name.
   * @param[in] fieldName Name of the field specification.
   * @return Name of the target field being updated.
   */
  string getTargetFieldName( string const & fieldName ) const;

  stdVector< string > m_fieldSpecificationNames;
};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_DYNAMIC_FIELDSPECIFCATION_HPP_ */
