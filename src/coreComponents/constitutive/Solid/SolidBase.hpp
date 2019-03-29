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
    static constexpr auto density0String  = "density0";
    static constexpr auto densityString  = "density";
    static constexpr auto deviatorStressString = "DeviatorStress";
    static constexpr auto meanStressString = "MeanStress";
  };


  struct StiffnessTensor
  {
    real64 c11;
    real64 c22;
    real64 c33;
    real64 c44;
    real64 c55;
    real64 c66;
    real64 c23;
    real64 c13;
    real64 c12;
    real64 c21;
    real64 c31;
    real64 c32;
  };
protected:

//  template< typename LEAFCLASS, typename POLICY=materialUpdatePolicy, typename ... ARGS >
//  void BatchUpdateKernel( ARGS && ... args );


  real64 m_density0;
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
