/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file SolidBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{

namespace constitutive
{



class SolidBase : public constitutive::ConstitutiveBase
{
public:
  SolidBase( string const & name,
             Group * const parent );

  virtual ~SolidBase() override;

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void StateUpdatePoint( localIndex const k,
                                 localIndex const q,
                                 R2SymTensor const & Dadt,
                                 R2Tensor const & Rot,
                                 real64 const dt,
                                 integer const updateStiffnessFlag ) = 0;

//  virtual void BatchUpdate( arrayView2d<real64 const> const & Dadt,
//                            arrayView2d<real64 const> const & Rot )


  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto defaultDensityString  = "defaultDensity";
    static constexpr auto densityString  = "density";

    static constexpr auto stressString = "stress";

  };

  real64   defaultDensity() const { return m_defaultDensity; }
  real64 & defaultDensity()       { return m_defaultDensity; }

  arrayView2d<real64>       const & density()       { return m_density; }
  arrayView2d<real64 const> const & density() const { return m_density; }

  arrayView2d<R2SymTensor>       const & getStress()       { return m_stress; }
  arrayView2d<R2SymTensor const> const & getStress() const { return m_stress; }

protected:

//  template< typename LEAFCLASS, typename POLICY=materialUpdatePolicy, typename ... ARGS >
//  void BatchUpdateKernel( ARGS && ... args );


  real64 m_defaultDensity;
  array2d<real64> m_density;

  array2d<R2SymTensor> m_stress;

};


//template< typename LEAFCLASS, typename POLICY, typename ... ARGS >
//void SolidBase::BatchUpdateKernel( ARGS && ... args )
//{
//  LaunchKernel<POLICY>( GEOSX_LAMBDA ( localIndex const k, localIndex const q )
//  {
//    LEAFCLASS::Compute( args... );
//  } );
//}

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_ */
