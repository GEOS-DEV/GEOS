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
template< typename EOS_TYPE > struct EOSCatalogName {};
template<> struct EOSCatalogName< PengRobinsonEOS > { static constexpr char const * catalogName() { return "PengRobinson";  } };
template<> struct EOSCatalogName< SoaveRedlichKwongEOS > { static constexpr char const * catalogName() { return "SoaveRedlichKwong";  } };

// Naming conventions
template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
string NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::catalogName()
{
  return EOSCatalogName< EOS_TYPE_LIQUID >::catalogName();
}

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::
NegativeTwoPhaseFlashModel( string const & name,
                            ComponentProperties const & componentProperties ):
  FunctionBase( name, componentProperties )
{}

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
typename NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::KernelWrapper
NegativeTwoPhaseFlashModel< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::createKernelWrapper() const
{
  return KernelWrapper( m_componentProperties.getNumberOfComponents() );
}

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
NegativeTwoPhaseFlashModelUpdate< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::
NegativeTwoPhaseFlashModelUpdate( integer const numComponents ):
  m_numComponents( numComponents )
{}

// Explicit instantiation of the model template.
template class NegativeTwoPhaseFlashModel< PengRobinsonEOS, PengRobinsonEOS >;
template class NegativeTwoPhaseFlashModel< SoaveRedlichKwongEOS, SoaveRedlichKwongEOS >;

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
