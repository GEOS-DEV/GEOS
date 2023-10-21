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

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::
NegativeTwoPhaseFlashModel( string const & name,
                            string_array const & componentNames,
                            array1d< real64 > const & componentMolarWeight,
                            ComponentProperties const & componentProperties ):
  FunctionBase( name,
                componentNames,
                componentMolarWeight,
                componentProperties )
{}

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
typename NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::KernelWrapper
NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        m_componentProperties );
}

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
NegativeTwoPhaseFlashModelUpdate< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::
NegativeTwoPhaseFlashModelUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                                  ComponentProperties const & componentProperties ):
  FunctionBaseUpdate( componentMolarWeight,
                      componentProperties )
{}

// Explicit instantiation of the model template.
template class NegativeTwoPhaseFlashModel< PengRobinsonEOS, PengRobinsonEOS >;
template class NegativeTwoPhaseFlashModel< SoaveRedlichKwongEOS, SoaveRedlichKwongEOS >;

REGISTER_CATALOG_ENTRY( FunctionBase, NegativeTwoPhaseFlashPRPR,
                        string const &,
                        string_array const &,
                        array1d< real64 > const &,
                        ComponentProperties const & )

REGISTER_CATALOG_ENTRY( FunctionBase, NegativeTwoPhaseFlashSRKSRK,
                        string const &,
                        string_array const &,
                        array1d< real64 > const &,
                        ComponentProperties const & )

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
