#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGE_DECOMPOSITION_HPP_ 
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGE_DECOMPOSITION_HPP_ 

#include <vector>
#include "../DamageSpectralUtilities.hpp"
#include "../PropertyConversions.hpp"
#include "../InvariantDecompositions.hpp"

namespace geosx
{
namespace constitutive
{

//Template class declaration
template< char const *DECOMPOSITION_TYPE >
class Decomposition
{};

//SPECTRAL
template<>
class Decomposition<"Spectral">
{

GEOSX_HOST_DEVICE
static void stressCalculator( real64 const ( &strain )[6],
                               real64 const damageFactor,
                               real64 ( & stress )[6],
                               real64 ( & stiffness )[6][6], 
                               real64 & sed ) const                                
  {
    
    strain[3] = strain[3]/2; // eigen-decomposition below does not use engineering strains
    strain[4] = strain[4]/2;
    strain[5] = strain[5]/2;

    real64 traceOfStrain = strain[0] + strain[1] + strain[2];

    real64 mu = stiffness.m_shearModulus[k];
    real64 lambda = conversions::bulkModAndShearMod::toFirstLame( stiffness.m_bulkModulus[k], mu );

    // get eigenvalues and eigenvectors

    real64 eigenValues[3] = {};
    real64 eigenVectors[3][3] = {};
    LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, strain );

    // tranpose eigenVectors matrix

    real64 temp[3][3] = {};
    LvArray::tensorOps::transpose< 3, 3 >( temp, eigenVectors );
    LvArray::tensorOps::copy< 3, 3 >( eigenVectors, temp );

    // get trace+ and trace-

    real64 tracePlus = fmax( traceOfStrain, 0.0 );
    real64 traceMinus = fmin( traceOfStrain, 0.0 );

    // build symmetric matrices of positive and negative eigenvalues

    real64 eigenPlus[6] = {};
    real64 eigenMinus[6] = {};
    real64 Itensor[6] = {};

    for( int i = 0; i < 3; i++ )
    {
      Itensor[i] = 1;
      eigenPlus[i] = fmax( eigenValues[i], 0.0 );
      eigenMinus[i] = fmin( eigenValues[i], 0.0 );
    }

    real64 positivePartOfStrain[6] = {};
    real64 negativePartOfStrain[6] = {};
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( positivePartOfStrain, eigenVectors, eigenPlus );
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( negativePartOfStrain, eigenVectors, eigenMinus );

    // stress

    real64 positiveStress[6] = {};
    real64 negativeStress[6] = {};
    LvArray::tensorOps::scaledCopy< 6 >( positiveStress, Itensor, lambda*tracePlus );
    LvArray::tensorOps::scaledCopy< 6 >( negativeStress, Itensor, lambda*traceMinus );

    LvArray::tensorOps::scaledAdd< 6 >( positiveStress, positivePartOfStrain, 2*mu );
    LvArray::tensorOps::scaledAdd< 6 >( negativeStress, negativePartOfStrain, 2*mu );

    LvArray::tensorOps::copy< 6 >( stress, negativeStress );
    LvArray::tensorOps::scaledAdd< 6 >( stress, positiveStress, damageFactor );

    // stiffness

    real64 IxITensor[6][6] = {};
    for( int i=0; i < 3; i++ )
    {
      for( int j=0; j < 3; j++ )
      {
        IxITensor[i][j] = 1.0;
      }
    }

    real64 cPositive[6][6] = {};
    real64 positiveProjector[6][6] = {};
    real64 negativeProjector[6][6] = {};

    PositiveProjectorTensor( eigenValues, eigenVectors, positiveProjector );
    NegativeProjectorTensor( eigenValues, eigenVectors, negativeProjector );

    LvArray::tensorOps::scaledCopy< 6, 6 >( cPositive, IxITensor, lambda*heaviside( traceOfStrain ));
    LvArray::tensorOps::scaledCopy< 6, 6 >( stiffness, IxITensor, lambda*heaviside( -traceOfStrain ));

    LvArray::tensorOps::scale< 6, 6 >( positiveProjector, 2*mu );
    LvArray::tensorOps::scale< 6, 6 >( negativeProjector, 2*mu );

    LvArray::tensorOps::add< 6, 6 >( cPositive, positiveProjector );
    LvArray::tensorOps::add< 6, 6 >( stiffness, negativeProjector );

    LvArray::tensorOps::scale< 6, 6 >( cPositive, damageFactor );
    LvArray::tensorOps::add< 6, 6 >( stiffness, cPositive );

    // compute strain energy density

    sed = 0.5 * lambda * tracePlus * tracePlus + mu * doubleContraction( positivePartOfStrain, positivePartOfStrain );

  }
};

template<>
class Decomposition<"VolDev">
{

  //voldev
  GEOSX_HOST_DEVICE
  static void stressCalculator( real64 const ( &strain )[6],
                                real64 const damageFactor,
                                real64 ( & stress )[6],
                                real64 ( & stiffness )[6][6],
                                real64 & sed ) const
  {

    real64 volStrain;
    real64 devStrain;
    real64 deviator[6];

    twoInvariant::strainDecomposition( strain,
                                       volStrain,
                                       devStrain,
                                       deviator );

    // degrade shear stiffness always
    // degrade volumetric stiffness when in tension

    stiffness.m_shearModulus *= damageFactor;

    if( volStrain > 0 )
    {
      stiffness.m_bulkModulus *= damageFactor;
    }

    // compute stress invariants and recompose full stress tensor

    real64 stressP = stiffness.m_bulkModulus * volStrain;
    real64 stressQ = 3 * stiffness.m_shearModulus * devStrain;

    twoInvariant::stressRecomposition( stressP,
                                       stressQ,
                                       deviator,
                                       stress );

    // update strain energy density
    // TODO: refactor as a proper history variable update.  the code below doesn't allow for rewinds.

    real64 sed = 0.5 * (stressQ * devStrain) / damageFactor;

    if( volStrain > 0 )
    {
      sed += 0.5 * (stressP * volStrain) / damageFactor;
    }

  }
};

template<>
class Decomposition<"None">
{
  //nosplit
  GEOSX_HOST_DEVICE
  static void stressCalculator( real64 const ( &strain )[6],
                                real64 const damageFactor,
                                real64 ( & stress )[6],
                                real64 ( & stiffness )[6][6],
                                real64 & sed ) const
  {

    LvArray::tensorOps::scale< 6 >( stress, damageFactor );

    stiffness.scaleParams( damageFactor );

    real64 traceOfStrain = strain[0] + strain[1] + strain[2];

    strain[3] = strain[3]/2; // to compute the energy using double contraction, we dont use the engineering strain
    strain[4] = strain[4]/2;
    strain[5] = strain[5]/2;

    sed = 0.5 * lambda * traceOfStrain * traceOfStrain + mu * doubleContraction( strain, strain );
  }
};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_DECOMPOSITION_HPP_ */