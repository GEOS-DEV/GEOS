/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TableRelativePermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITY_HPP
#define GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITY_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "managers/Functions/TableFunction.hpp"

namespace geosx
{
namespace constitutive
{

class TableRelativePermeabilityUpdate final : public RelativePermeabilityBaseUpdate
{
public:

  TableRelativePermeabilityUpdate( arrayView1d< TableFunctionKernelWrapper const > const & relPermTableKernelWrappers,
                                   arrayView1d< integer const > const & phaseTypes,
                                   arrayView1d< integer const > const & phaseOrder,
                                   arrayView3d< real64 > const & phaseRelPerm,
                                   arrayView4d< real64 > const & dPhaseRelPerm_dPhaseVolFrac )
    : RelativePermeabilityBaseUpdate( phaseTypes,
                                      phaseOrder,
                                      phaseRelPerm,
                                      dPhaseRelPerm_dPhaseVolFrac ),
    m_relPermTableKernelWrappers( relPermTableKernelWrappers )
  {}

  /// Default copy constructor
  TableRelativePermeabilityUpdate( TableRelativePermeabilityUpdate const & ) = default;

  /// Default move constructor
  TableRelativePermeabilityUpdate( TableRelativePermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  TableRelativePermeabilityUpdate & operator=( TableRelativePermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  TableRelativePermeabilityUpdate & operator=( TableRelativePermeabilityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void Compute( arraySlice1d< real64 const > const & phaseVolFraction,
                        arraySlice1d< real64 > const & phaseRelPerm,
                        arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) const override;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void Update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const > const & phaseVolFraction ) const override
  {
    Compute( phaseVolFraction,
             m_phaseRelPerm[k][q],
             m_dPhaseRelPerm_dPhaseVolFrac[k][q] );
  }

private:

  arrayView1d< TableFunctionKernelWrapper const > m_relPermTableKernelWrappers;

};

class TableRelativePermeability : public RelativePermeabilityBase
{
public:

  TableRelativePermeability( std::string const & name, dataRepository::Group * const parent );

  virtual ~TableRelativePermeability() override;

  static std::string CatalogName() { return "TableRelativePermeability"; }

  virtual string getCatalogName() const override { return CatalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = TableRelativePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr auto relPermTableNamesString = "relPermTableNames";

    using ViewKey = dataRepository::ViewKey;
    ViewKey relPermTableNames = { relPermTableNamesString };

  } vieKeysTableRelativePermeability;

protected:

  virtual void PostProcessInput() override;

  virtual void InitializePreSubGroups( Group * const ) override;

private:

  // Validate the relative permeability table provided in input (increasing phase vol frac and rel perm, etc)
  void ValidateRelativePermeabilityTable( TableFunction const & relPermTable ) const;

  /// Relative permeability table names (one for each phase)
  array1d< string > m_relPermTableNames;

  /// Relative permeability kernel wrapper
  array1d< TableFunctionKernelWrapper > m_relPermTableKernelWrappers;

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TableRelativePermeabilityUpdate::
  Compute( arraySlice1d< real64 const > const & phaseVolFraction,
           arraySlice1d< real64 > const & phaseRelPerm,
           arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  real64 volFrac[ 1 ] = { 0.0 };
  real64 relPermDerivative[ 1 ] = { 0.0 };
  for( localIndex ip = 0; ip < phaseVolFraction.size(); ++ip )
  {
    volFrac[0] = phaseVolFraction[ip];
    m_relPermTableKernelWrappers[ip].Compute( volFrac, phaseRelPerm[ip], relPermDerivative );
    dPhaseRelPerm_dPhaseVolFrac[ip][ip] = relPermDerivative[0];
  }
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITY_HPP
