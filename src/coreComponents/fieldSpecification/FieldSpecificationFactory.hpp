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
 * @file FieldSpecificationFactory.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONFACTORY_HPP
#define GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONFACTORY_HPP


#include "common/TypeDispatch.hpp"
#include "dataRepository/Group.hpp"
#include "FieldSpecificationABC.hpp"

namespace geos
{

/**
 * @brief Generate FieldSpecifications based on the given "higher-level" specification
 * @tparam SPEC_TYPE The type of the high-level specification
 * @param fs The high-level specification used as a blueprint to create FieldSpecification
 * @param manager The parent to store the created FieldSpecifications
 */
template< typename SPEC_TYPE >
void expandFieldSpecification( SPEC_TYPE const & fs,
                               dataRepository::Group & manager ) = delete;

class FieldSpecificationProcessorRegistry
{
public:

  /**
   * @brief Base Processor class for transforming "high-level" specification
   *        into FieldSpecification(s)
   */
  class ProcessorBase
  {
public:
    /**
     * @brief Generate FieldSpecifications based on the given "higher-level" specification
     * @tparam SPEC_TYPE The type of the high-level specification
     * @param fs The high-level specification used as a blueprint to create FieldSpecification
     * @param manager The parent to store the created FieldSpecifications
     */
    virtual void expandFieldSpecification( FieldSpecificationABC const & fs,
                                           dataRepository::Group & GEOS_UNUSED_PARAM( manager ) ) const
    { GEOS_ERROR( GEOS_FMT( "Processor not implemented for field specification of type '{}'.", fs.getCatalogName() ), fs.getDataContext() ); }
protected:
    ProcessorBase() {}
  };

  /**
   * @brief Class for a specific "high-level" specification type to process its objects
   *        into one or multiple equivalent FieldSpecification
   */
  template< typename SPEC_TYPE >
  class Processor final : public ProcessorBase
  {
public:

    /**
     * @brief Add the processors to the static list. Called before main() when a
     *        REGISTER_FIELD_SPECIFICATION_PROCESSOR( SPEC_TYPE ) is in a cpp.
     */
    Processor(): ProcessorBase()
    { s_processors.emplace( SPEC_TYPE::catalogName(), this ); }

    /**
     * @copydoc ProcessorBase::expandFieldSpecification
     */
    void expandFieldSpecification( FieldSpecificationABC const & fs,
                                   dataRepository::Group & manager ) const override
    { geos::expandFieldSpecification( dynamic_cast< SPEC_TYPE const & >( fs ), manager ); }
  };

  /**
   * @return the list of field specification processors.
   */
  static stdMap< string, ProcessorBase const * > const & getProcessors()
  { return s_processors; }


private:

  /**
   * @brief Storage of field specification processors.
   */
  static stdMap< string, ProcessorBase const * > s_processors;

};

#define REGISTER_FIELD_SPECIFICATION_PROCESSOR( SPEC_TYPE ) \
  namespace { \
  GEOS_MAYBE_UNUSED FieldSpecificationProcessorRegistry::Processor< SPEC_TYPE > g_processorOf ## SPEC_TYPE; \
  }

} // namespace geos


#endif //GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONFACTORY_HPP
