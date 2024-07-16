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
 * @file TerseInputRegistry.hpp
 */


#ifndef GEOS_DATAREPOSITORY_TERSEINPUTREGISTRY_HPP_
#define GEOS_DATAREPOSITORY_TERSEINPUTREGISTRY_HPP_

#include "InputExtension.hpp"

namespace geos
{

namespace dataRepository
{

namespace inputExtension
{

template < typename NodeType >
class TerseInputRegistry
{
public:
  using map_type = std::map< string, std::vector< Rule< NodeType > > >;
  // Prevent instantiation
  TerseInputRegistry() = delete;

  // Function to add a new extension rule
  static void addExtensionRule(string const & tag, const Rule< NodeType > & rule)
  {
    getTerseExtensionRules()[tag].push_back( rule );
  }

  // Static function to retrieve the set of currently registered expansion rules
  static const map_type & getInputExtensionRules( )
  {
    return getTerseExtensionRules();
  }
private:
  static map_type & getTerseExtensionRules()
  {
    static map_type extensionRules;
    return extensionRules;
  }
};

} // namespace inputExtension

} // namespace dataRepository

} // namespace geos

#define REGISTER_SUGAR_EXTENSION_RULE( sugar_tag, document_node, rule_name, configure_function ) \
  namespace { \
    using namespace geos::dataRepository::inputExtension; \
    struct AutoRegister_##rule_name { \
      AutoRegister_##rule_name() { \
        Rule< document_node > rule; \
        configure_function( rule ); \
        TerseInputRegistry< document_node >::addExtensionRule( sugar_tag, rule ); \
      } \
    } autoRegister_##rule_name; \
  }

// Macro to define a rule, configure it using a provided function, and add it to the extension rules registry
#define REGISTER_CATALOG_EXTENSION_RULE( group_class, rule_name, configure_function ) \
  REGISTER_SUGAR_EXTENSION_RULE( group_class::CataloName(), default_document_node, rule_name, configure_function )

#endif