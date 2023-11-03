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

// Naming conventions
template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
string NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::catalogName()
{
  return EOS_TYPE_LIQUID::catalogName();
}

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::
NegativeTwoPhaseFlashModel( string const & name,
                            ComponentProperties const & componentProperties ):
  FunctionBase( name,
                componentProperties )
{}

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
typename NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::KernelWrapper
NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::createKernelWrapper() const
{
  return KernelWrapper( m_componentProperties );
}

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
NegativeTwoPhaseFlashModelUpdate< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::
NegativeTwoPhaseFlashModelUpdate( ComponentProperties const & componentProperties ):
  FunctionBaseUpdate( componentProperties ),
  m_numComponents( componentProperties.getNumberOfComponents())
{}

// Explicit instantiation of the model template.
template class NegativeTwoPhaseFlashModel<
    CubicEOSPhaseModel< PengRobinsonEOS >,
    CubicEOSPhaseModel< PengRobinsonEOS > >;
template class NegativeTwoPhaseFlashModel<
    CubicEOSPhaseModel< SoaveRedlichKwongEOS >,
    CubicEOSPhaseModel< SoaveRedlichKwongEOS > >;

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
