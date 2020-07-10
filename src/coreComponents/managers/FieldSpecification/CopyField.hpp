
/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#ifndef COPYFIELD_HPP_
#define COPYFIELD_HPP_



#include <string>
#include <limits>

#include <managers/Tasks/TaskBase.hpp>

namespace geosx
{

/**
 * @class CopyField
 * @brief A Task to copy a field.
 */
class CopyField : public TaskBase
{
public:
  /// @copydoc geosx::dataRepository::Group::Group
  explicit CopyField( std::string const & name,
                      Group * const parent );

  /// Destructor
  virtual ~CopyField() override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string CatalogName() { return "CopyField"; }

  /// @copydoc geosx::TaskBase::Execute
  void Execute( real64 const time_n,
                real64 const dt,
                integer const cycleNumber,
                integer const eventCounter,
                real64 const eventProgress,
                Group * domain ) override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr auto fromString = "from";
    static constexpr auto toString = "to";
    static constexpr auto targetRegionsString = "targetRegions";
  };
  /// @endcond
private:
  /// Name of the field to copy
  string m_from;

  /// Name of the new field
  string m_to;

  /// Regions on which the copy will be done
  string_array m_targetRegions;
};
} /* namespace geosx */


#endif /* COPYFIELD_HPP_ */
