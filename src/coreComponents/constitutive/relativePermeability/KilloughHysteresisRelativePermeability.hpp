//
// Created by root on 11/9/22.
//

#ifndef GEOSX_KILLOUGHHYSTERESISRELATIVEPERMEABILITY_HPP
#define GEOSX_KILLOUGHHYSTERESISRELATIVEPERMEABILITY_HPP

#include "constitutive/KilloughHysteresis.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx {
namespace constitutive{

class KilloughHysteresisRelativePermeability : public KilloughHysteresis
{

public:

  KilloughHysteresisRelativePermeability(std::string const& name, Group * const parent);

  class KilloughHysteresisRelativePermeabilityKernel : public KernelKilloughHysteresisBase
  {

  public:

    KilloughHysteresisRelativePermeabilityKernel( real64 const & jerauldParam_a,
                                                  real64 const & jerauldParam_b,
                                                  real64 const & killoughCurvatureParam,
                                                  arrayView1d< real64 const > const & landParam,
                                                  const arrayView1d< const geosx::real64 > & drainageMinPhaseVolFraction,
                                                  const arrayView1d< const geosx::real64 > & imbibitionMinPhaseVolFraction,
                                                  const arrayView1d< const geosx::real64 > & drainageRelPermEndPoint,
                                                  const arrayView1d< const geosx::real64 > & imbibitionRelPermEndPoint,
                                                  const arrayView1d< const geosx::real64 > & drainageMaxPhaseVolFraction,
                                                  const arrayView1d< const geosx::real64 > & imbibitionMaxPhaseVolFraction,
                                                  const arrayView3d< geosx::real64, relperm::USD_RELPERM > & phaseTrapppedVolFrac );

    GEOSX_HOST_DEVICE
    inline void
    computeImbibitionWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                     TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
//                                     real64 const & jerauldParam_a,
//                                     real64 const & jerauldParam_b,
                                     real64 const & landParam,
                                     real64 const & phaseVolFraction,
                                     real64 const & phaseMinHistoricalVolFraction,
                                     real64 const & imbibitionPhaseMinWettingVolFraction,
                                     real64 const & drainagePhaseMaxVolFraction,
                                     real64 const & imbibitionPhaseMaxVolFraction,
                                     real64 const & drainageRelPermEndPoint,
                                     real64 const & imbibitionRelPermEndPoint,
                                     real64 & phaseTrappedVolFrac,
                                     real64 & phaseRelPerm,
                                     real64 & dPhaseRelPerm_dPhaseVolFrac ) const;
    GEOSX_HOST_DEVICE
    inline void
    computeImbibitionNonWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                        TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
//                                        real64 const & jerauldParam_a,
//                                        real64 const & jerauldParam_b,
                                        real64 const & landParam,
                                        real64 const & phaseVolFraction,
                                        real64 const & phaseMaxHistoricalVolFraction,
                                        real64 const & drainagePhaseMinVolFraction,
                                        real64 const & imbibitionPhaseMinVolFraction,
                                        real64 const & drainagePhaseMaxVolFraction,
                                        real64 const & drainageRelPermEndPoint,
                                        real64 & phaseTrappedVolFrac,
                                        real64 & phaseRelPerm,
                                        real64 & dPhaseRelPerm_dPhaseVolFrac ) const;




  private:

    /// Curvature parameter introduced for wetting phase hysteresis in Killough
    real64 const m_killoughCurvatureParam;
  };

struct viewKeyStruct : public KilloughHysteresis::viewKeyStruct
{
  static constexpr char const * killoughCurvatureParameterString() { return "killoughCurvatureParameter"; }
};


  void postProcessInput() override;
//  KernelWrapper createKernelWrapper();
  static std::string catalogName() { return "KilloughHysteresisRelativePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }
private:

  real64 m_killoughCurvatureParam;

};

GEOSX_HOST_DEVICE
inline void
KilloughHysteresisRelativePermeability::KilloughHysteresisRelativePermeabilityKernel::
computeImbibitionWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                 TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
//                                 real64 const & jerauldParam_a,
//                                 real64 const & jerauldParam_b,
                                 real64 const & landParam,
                                 real64 const & phaseVolFraction,
                                 real64 const & phaseMinHistoricalVolFraction,
                                 real64 const & imbibitionPhaseMinWettingVolFraction,
                                 real64 const & drainagePhaseMaxVolFraction,
                                 real64 const & imbibitionPhaseMaxVolFraction,
                                 real64 const & drainageRelPermEndPoint,
                                 real64 const & imbibitionRelPermEndPoint,
                                 real64 & phaseTrappedVolFrac,
                                 real64 & phaseRelPerm,
                                 real64 & dPhaseRelPerm_dPhaseVolFrac ) const
{

  // Step 0: preparing keypoints in the (S,kr) plan
  // if consistent, S should be equal to 1 - imbibitionPhaseMinVolNonWettingFraction for two-phase flow
  // (but wetting and nonwetting phase hysteresis are implemented in a decoupled fashion)
  real64 const S = phaseVolFraction;
  real64 const Smxi = imbibitionPhaseMaxVolFraction;
  real64 const Smxd = drainagePhaseMaxVolFraction;

  // Swc is the common end min endpoint saturation for wetting curves
  real64 const Swc = imbibitionPhaseMinWettingVolFraction;

  if( S <= Swc )
  {
    phaseRelPerm = 0.0;
    dPhaseRelPerm_dPhaseVolFrac = 0.0;
  }
  else if( S >= Smxd )
  {
    phaseRelPerm = drainageRelPermEndPoint;
    dPhaseRelPerm_dPhaseVolFrac = 0.0;
  }
  else
  {
    real64 const krwei = imbibitionRelPermEndPoint;
    real64 const krwedAtSmxi = drainageRelPermKernelWrapper.compute( &Smxi );

    // Step 1: Compute the new end point

    // Step 1.a: get the value at the max non-wetting residual value
    real64 const deltak = krwei - krwedAtSmxi;

    // Step 1.b: get the trapped from wetting data
    real64 const Shy = ( phaseMinHistoricalVolFraction > Swc ) ? phaseMinHistoricalVolFraction : Swc;
    real64 const A = 1 + m_jerauldParam_a * ( Shy - Swc );
    real64 const numerator = Shy - Smxd;
    real64 const denom = A + landParam * pow( ( Smxd - Shy ) / ( Smxd - Swc ), 1 + m_jerauldParam_b/landParam );
    real64 const Scrt = Smxd + numerator / denom;

    // Step 1.c: find the new endpoint
    // this is the saturation for the scanning curve endpoint
    real64 const krwedAtScrt = drainageRelPermKernelWrapper.compute( &Scrt );
    real64 const krwieStar = krwedAtScrt
                             + deltak * pow( ( Smxd - Scrt ) / LvArray::math::max( minScriMinusScrd, ( Smxd - Smxi ) ), m_killoughCurvatureParam );

    // Step 2: get the normalized value of saturation
    real64 const ratio = ( Smxi - Swc ) / ( Scrt - Shy );
    real64 const Snorm = Smxi - ( Scrt - S ) * ratio; // normalized saturation from equation 2.166
    real64 const dSnorm_dS =  ratio;
    real64 dkri_dSnorm = 0.0;
    real64 const krwiAtSnorm = imbibitionRelPermKernelWrapper.compute( &Snorm, &dkri_dSnorm );
    real64 const dkriAtSnorm_dS = dkri_dSnorm * dSnorm_dS;

    // Step 3: Get the final value at evaluated saturation
    real64 const krdAtShy = drainageRelPermKernelWrapper.compute( &Shy );
    real64 const imbibitionRelPermRatio = (krwieStar - krdAtShy) / krwei;

    phaseRelPerm = krdAtShy + krwiAtSnorm * imbibitionRelPermRatio;
    dPhaseRelPerm_dPhaseVolFrac = dkriAtSnorm_dS * imbibitionRelPermRatio;
  }


  // Updating the trapped phase volume fraction
  // Residual water is constant in Killough hysteresis model
  phaseTrappedVolFrac = LvArray::math::min( Swc, phaseVolFraction );
}


GEOSX_HOST_DEVICE
inline void
KilloughHysteresisRelativePermeability::KilloughHysteresisRelativePermeabilityKernel::
computeImbibitionNonWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                    TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
//                                    real64 const & jerauldParam_a,
//                                    real64 const & jerauldParam_b,
                                    real64 const & landParam,
                                    real64 const & phaseVolFraction,
                                    real64 const & phaseMaxHistoricalVolFraction,
                                    real64 const & drainagePhaseMinVolFraction,
                                    real64 const & imbibitionPhaseMinVolFraction,
                                    real64 const & drainagePhaseMaxVolFraction,
                                    real64 const & drainageRelPermEndPoint,
                                    real64 & phaseTrappedVolFrac,
                                    real64 & phaseRelPerm,
                                    real64 & dPhaseRelPerm_dPhaseVolFrac ) const
{
  // note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

  // Step 1: for a given value of the max historical saturation, Shy, compute the trapped critical saturation, Scrt,
  //         using Land's method. The calculation includes the modifications from Jerauld. This is equation 2.162 from
  //         the IX technical description.
  real64 const S = phaseVolFraction;
  real64 const Scri = imbibitionPhaseMinVolFraction;
  real64 const Scrd = drainagePhaseMinVolFraction;
  real64 const Smx = drainagePhaseMaxVolFraction;
  real64 const Shy = phaseMaxHistoricalVolFraction < Smx ? phaseMaxHistoricalVolFraction : Smx; // to make sure that Shy < Smax
  real64 Scrt = 0;
  computeTrappedCriticalPhaseVolFraction( Scrd,
                                          Shy,
                                          Smx,
                                          landParam,
                                          Scrt );

  if( S <= Scrt )  // S is below the trapped critical saturation, so the relperm is zero
  {
    phaseRelPerm = 0.0;
    dPhaseRelPerm_dPhaseVolFrac = 0.0;
  }
  else if( S >= Smx ) // S is above the max saturation, so we just skip the rest and set the relperm to the endpoint
  {
    phaseRelPerm = drainageRelPermEndPoint;
    dPhaseRelPerm_dPhaseVolFrac = 0.0;
  }
  else
  {
    // Step 2: compute the normalized saturation, S_norm, at which the imbibition relperm curve will be evaluated.
    //         This is equation 2.166 from the IX technical description.
    real64 const ratio = ( Smx - Scri ) / ( Shy - Scrt );
    real64 const Snorm = Scri + ( S - Scrt ) * ratio; // normalized saturation from equation 2.166
    real64 const dSnorm_dS = ratio;

    // Step 3: evaluate the imbibition relperm, kri(Snorm), at the normalized saturation, Snorm.
    real64 dkri_dSnorm = 0;
    real64 const kriAtSnorm = imbibitionRelPermKernelWrapper.compute( &Snorm, &dkri_dSnorm );
    real64 const dkriAtSnorm_dS = dkri_dSnorm * dSnorm_dS;

    // Step 4: evaluate the drainage relperm, krd(Shy), at the max hystorical saturation, Shy.
    real64 const krdAtShy = drainageRelPermKernelWrapper.compute( &Shy );

    // Step 5: evaluate the drainage relperm, krd(Smx), at the max drainage saturation, Smx.
    real64 const krdAtSmx = drainageRelPermEndPoint;

    // Step 6: apply the formula blending drainage and imbibition relperms from the Killough model.
    //         This equation 2.165 from the IX technical description.
    real64 const drainageRelPermRatio = krdAtShy / krdAtSmx;
    phaseRelPerm = kriAtSnorm * drainageRelPermRatio;
    dPhaseRelPerm_dPhaseVolFrac = dkriAtSnorm_dS * drainageRelPermRatio;
  }

  //updating trapped vol fraction
  phaseTrappedVolFrac = LvArray::math::min( phaseVolFraction, Scrt );
}

} //end namespace constitutive
} //end namespace geosx

#endif //GEOSX_KILLOUGHHYSTERESISRELATIVEPERMEABILITY_HPP
