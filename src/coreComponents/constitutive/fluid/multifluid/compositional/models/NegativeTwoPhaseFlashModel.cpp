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
 * @file NegativeTwoPhaseFlashModel.cpp
 */

#include "NegativeTwoPhaseFlashModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

string NegativeTwoPhaseFlashModel::catalogName()
{
  return "TwoPhase";
}

NegativeTwoPhaseFlashModel::NegativeTwoPhaseFlashModel( string const & name,
                                                        ComponentProperties const & componentProperties,
                                                        EquationOfState const & equationOfState ):
  FunctionBase( name, componentProperties, equationOfState )
{}

typename NegativeTwoPhaseFlashModel::KernelWrapper
NegativeTwoPhaseFlashModel::createKernelWrapper() const
{
  return KernelWrapper( m_componentProperties.getNumberOfComponents(), 0, 1 );
}

NegativeTwoPhaseFlashModelUpdate::NegativeTwoPhaseFlashModelUpdate( integer const numComponents,
                                                                    integer const liquidIndex,
                                                                    integer const vapourIndex ):
  m_numComponents( numComponents ),
  m_liquidIndex( liquidIndex ),
  m_vapourIndex( vapourIndex )
{}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
