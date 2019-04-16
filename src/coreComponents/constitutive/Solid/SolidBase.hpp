/**
 * @file SolidBase.hpp
 */

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_
#define SRC_CORECOMPONENTS_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{

namespace constitutive
{



class SolidBase : public constitutive::ConstitutiveBase
{
public:
  SolidBase( string const & name,
             ManagedGroup * const parent );

  virtual ~SolidBase() override;

  virtual void DeliverClone( string const & name,
                             ManagedGroup * const parent,
                             std::unique_ptr<ConstitutiveBase> & clone ) const override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void StateUpdatePoint( localIndex const k,
                                 localIndex const q,
                                 R2SymTensor const & Dadt,
                                 R2Tensor const & Rot,
                                 integer const systemAssembleFlag ) = 0;

//  virtual void BatchUpdate( arrayView2d<real64 const> const & Dadt,
//                            arrayView2d<real64 const> const & Rot )


  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto defaultDensityString  = "defaultDensity";
    static constexpr auto densityString  = "density";
    static constexpr auto deviatorStressString = "DeviatorStress";
    static constexpr auto meanStressString = "MeanStress";
  };

  real64   defaultDensity() const { return m_defaultDensity; }
  real64 & defaultDensity()       { return m_defaultDensity; }

  arrayView2d<real64>       const & density()       { return m_density; }
  arrayView2d<real64 const> const & density() const { return m_density; }

  arrayView2d<real64>        const & meanStress()       { return m_meanStress; }
  arrayView2d<real64 const > const & meanStress() const { return m_meanStress; }

  arrayView2d<R2SymTensor>       const & deviatorStress()       { return m_deviatorStress; }
  arrayView2d<R2SymTensor const> const & deviatorStress() const { return m_deviatorStress; }

protected:

//  template< typename LEAFCLASS, typename POLICY=materialUpdatePolicy, typename ... ARGS >
//  void BatchUpdateKernel( ARGS && ... args );


  real64 m_defaultDensity;
  array2d<real64> m_density;

  array2d<real64> m_meanStress;
  array2d<R2SymTensor> m_deviatorStress;

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

#endif /* SRC_CORECOMPONENTS_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_ */
