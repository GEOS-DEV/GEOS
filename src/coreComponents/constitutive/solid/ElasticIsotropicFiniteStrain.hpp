/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file ElasticIsotropicFiniteStrain.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_ELASTICISOTROPICFINITESTRAIN_HPP_
#define GEOS_CONSTITUTIVE_SOLID_ELASTICISOTROPICFINITESTRAIN_HPP_

#include "ElasticIsotropic.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullTensor.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos{

namespace constitutive {

/**
 * @class ElasticIsotropicFiniteStrainUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class ElasticIsotropicFiniteStrainUpdates : public ElasticIsotropicUpdates
{
public:
  ElasticIsotropicFiniteStrainUpdates(
    arrayView1d< real64 const > const & bulkModulus,
    arrayView1d< real64 const > const & shearModulus,
    arrayView1d< real64 const > const & thermalExpansionCoefficient,
    arrayView3d< real64, solid::STRESS_USD > const & newStress,
    arrayView3d< real64, solid::STRESS_USD > const & oldStress,
    bool const & disableInelasticity )
  : ElasticIsotropicUpdates(bulkModulus, shearModulus, thermalExpansionCoefficient, newStress, oldStress, disableInelasticity)
  {}

  /// Default copy constructor
  ElasticIsotropicFiniteStrainUpdates( ElasticIsotropicFiniteStrainUpdates const & ) = default;

  /// Default move constructor
  ElasticIsotropicFiniteStrainUpdates( ElasticIsotropicFiniteStrainUpdates && ) = default;

  /// Deleted default constructor
  ElasticIsotropicFiniteStrainUpdates() = delete;

  /// Deleted copy assignment operator
  ElasticIsotropicFiniteStrainUpdates & operator=( ElasticIsotropicFiniteStrainUpdates const & ) = delete;

  /// Deleted move assignment operator
  ElasticIsotropicFiniteStrainUpdates & operator=( ElasticIsotropicFiniteStrainUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullTensor;

  // returns first Piola-Kirchhoff stress which is asymmetric
  GEOS_HOST_DEVICE
  void finiteStrainNoStateUpdate_StressOnly(localIndex const k, localIndex const q,
                                            real64 const (&totalElasticStrain)[6],
                                            real64 const (&fInv)[3][3],
                                            real64 (&firstPiolaStress)[3][3],
                                            real64 (&kirchhoffStress)[6]
                                           ) const;

  GEOS_HOST_DEVICE
  void finiteStrainUpdate_StressOnly(localIndex const k, localIndex const q,
                                     real64 const (&totalElasticStrain)[6],
                                     real64 const (&fInv)[3][3],
                                     real64 (&firstPiolaStress)[3][3]
                                    ) const;

  GEOS_HOST_DEVICE
  void finiteStrainNoStateUpdate(localIndex const k, localIndex const q,
                                 real64 const (&elasticDeformGrad)[3][3],
                                 real64 (&firstPiolaStress)[3][3],
                                 real64 (&kirchhoffStress)[6],
                                 real64 (&stiffness)[9][9]) const;

  GEOS_HOST_DEVICE
  void finiteStrainUpdate(localIndex const k, localIndex const q,
                          real64 const (&elasticDeformGrad)[3][3],
                          real64 (&firstPiolaStress)[3][3],
                          real64 (&stiffness)[9][9]) const;

  GEOS_HOST_DEVICE
  void finiteStrainUpdate(localIndex const k, localIndex const q,
                          real64 const (&elasticDeformGrad)[3][3],
                          real64 (&firstPiolaStress)[3][3],
                          DiscretizationOps & stiffness) const;

  GEOS_HOST_DEVICE
  void computeLogElasticStrain(real64 const (&elasticDeformGrad)[3][3],
                               real64 (&elasticStrain)[6],
                               real64 (&eigenValues)[3],
                               real64 (&eigenVectors)[3][3],
                               real64 (&eigenVectorsT)[3][3]) const;

private:
  GEOS_HOST_DEVICE
  void computeMaterialTangentColumn(real64 const (&deltaCe)[6], real64 const (&M_hat)[6],
                                    real64 const (&Q)[3][3], real64 const (&Q_T)[3][3],
                                    real64 const (&dtau_dEe)[6][6], real64 const (&fInvT)[3][3],
                                    real64 (&deltaP_mat)[3][3]) const;
};

GEOS_HOST_DEVICE
inline
void ElasticIsotropicFiniteStrainUpdates::computeLogElasticStrain(
    real64 const (&elasticDeformGrad)[3][3], real64 (&elasticStrain)[6],
    real64 (&eigenValues)[3], real64 (&eigenVectors)[3][3], real64 (&eigenVectorsT)[3][3]) const
{
  real64 C_e[3][3] = {}; // Right Cauchy tensor
  LvArray::tensorOps::Rij_eq_AkiAkj<3, 3>(C_e, elasticDeformGrad);

  real64 Ce_symmetric[6] = {};
  LvArray::tensorOps::denseToSymmetric<3>(Ce_symmetric, C_e);

  // initial eigen vector are stored in rows
  LvArray::tensorOps::symEigenvectors<3>(eigenValues, eigenVectorsT, Ce_symmetric);
  LvArray::tensorOps::transpose<3, 3>(eigenVectors, eigenVectorsT);

  // matrix logarithm
  real64 logLambda[6] = { {0} };
  for( int i = 0; i < 3; i++ )
  {
    logLambda[i] = 0.5 * log(eigenValues[i]);
  }

  // elastic strain
  // E_e = 0.5 * Q * ln(Lam) * Q^T
  LvArray::tensorOps::Rij_eq_AikSymBklAjl<3>(elasticStrain, eigenVectors, logLambda);

  // scale shear components by 2.0 to use in small strain material model where stiffness is expressed in Voigt notation
  elasticStrain[3] *= 2.0;
  elasticStrain[4] *= 2.0;
  elasticStrain[5] *= 2.0;
}

GEOS_HOST_DEVICE
inline
void ElasticIsotropicFiniteStrainUpdates::computeMaterialTangentColumn(
    real64 const (&deltaCe)[6], real64 const (&M_hat)[6], real64 const (&Q)[3][3],
    real64 const (&Q_T)[3][3], real64 const (&dtau_dEe)[6][6], real64 const (&fInvT)[3][3],
    real64 (&deltaP_mat)[3][3]) const
{
  // Q^T * deltaCe * Q
  real64 deltaCe_hat[6] = {};
  LvArray::tensorOps::Rij_eq_AikSymBklAjl<3>(deltaCe_hat, Q_T, deltaCe);

  // M \circ Q^T * deltaCe * Q
  real64 deltaLnCe_hat[6] = {};
  LvArray::tensorOps::hadamardProduct<6>(deltaLnCe_hat, M_hat, deltaCe_hat);

  // derivative of matrix logarithm deltaEe = dEe / dFe = 1/2 * Q * (M \circ Q^T * deltaCe * Q) * Q^T
  real64 deltaLnCe[6] = {};
  real64 deltaEe[6] = {};
  LvArray::tensorOps::Rij_eq_AikSymBklAjl<3>(deltaLnCe, Q, deltaLnCe_hat);
  LvArray::tensorOps::copy<6>(deltaEe, deltaLnCe);
  LvArray::tensorOps::scale<6>(deltaEe, 0.5);

  // scale shear components for stiffness in voigt notation
  deltaEe[3] *= 2.0;
  deltaEe[4] *= 2.0;
  deltaEe[5] *= 2.0;

  // double contraction dtau / dEe * deltaEe
  real64 dtau_voigt[6] = {};
  real64 deltaTau[3][3] = {};
  LvArray::tensorOps::Ri_eq_AijBj<6, 6>(dtau_voigt, dtau_dEe, deltaEe);
  LvArray::tensorOps::symmetricToDense<3>(deltaTau, dtau_voigt);

  // deltaTau * F^-T
  LvArray::tensorOps::Rij_eq_AikBkj<3, 3, 3>(deltaP_mat, deltaTau, fInvT);
}

GEOS_HOST_DEVICE
inline
void ElasticIsotropicFiniteStrainUpdates::finiteStrainNoStateUpdate_StressOnly(localIndex const k, localIndex const q,
                                                                               real64 const (&totalElasticStrain)[6],
                                                                               real64 const (&fInv)[3][3],
                                                                               real64 (&firstPiolaStress)[3][3],
                                                                               real64 (&kirchhoffStress)[6]
                                                                              ) const
{
  this->smallStrainNoStateUpdate_StressOnly(k, q, totalElasticStrain, kirchhoffStress);

  LvArray::tensorOps::Rij_eq_symAikBjk<3>(firstPiolaStress, kirchhoffStress, fInv);
}

GEOS_HOST_DEVICE
inline
void ElasticIsotropicFiniteStrainUpdates::finiteStrainUpdate_StressOnly(localIndex const k, localIndex const q,
                                                                        real64 const (&totalElasticStrain)[6],
                                                                        real64 const (&fInv)[3][3],
                                                                        real64 (&firstPiolaStress)[3][3]
                                                                       ) const
{
  real64 kirchhoffStress[6] = {};
  finiteStrainNoStateUpdate_StressOnly(k, q, totalElasticStrain, fInv, firstPiolaStress, kirchhoffStress);

  // Can only save the symmetric kirchhoff stress right now
  LvArray::tensorOps::copy<6>(m_oldStress[k][q], m_newStress[k][q]);
  LvArray::tensorOps::copy<6>(m_newStress[k][q], kirchhoffStress);
}

GEOS_HOST_DEVICE
inline
void ElasticIsotropicFiniteStrainUpdates::finiteStrainNoStateUpdate(localIndex const k, localIndex const q,
                                                                   real64 const (&elasticDeformGrad)[3][3],
                                                                   real64 (&firstPiolaStress)[3][3],
                                                                   real64 (&kirchhoffStress)[6],
                                                                   real64 (&stiffness)[9][9]) const
{
  real64 fInv[3][3] = {};
  real64 fInvT[3][3] = {};

  real64 totalElasticStrain[6] = {};

  real64 eigenValues[3] = {};
  real64 eigenVectors[3][3] = {};
  real64 eigenVectorsT[3][3] = {};

  LvArray::tensorOps::invert<3>(fInv, elasticDeformGrad);
  LvArray::tensorOps::transpose<3, 3>(fInvT, fInv);

  computeLogElasticStrain(elasticDeformGrad, totalElasticStrain, eigenValues, eigenVectors, eigenVectorsT);

  finiteStrainNoStateUpdate_StressOnly(k, q, totalElasticStrain, fInv, firstPiolaStress, kirchhoffStress);

  real64 smallStrainStiffness[6][6] = {};
  this->getElasticStiffness(k, q, smallStrainStiffness);

  // derivative of log matrix: \partial ln(C_e) / \partial C_e, C_e = F^T F
  // build symmetric spectral of log matrix derivative
  real64 M_hat[6] = {};
  M_hat[0] = 1.0 / eigenValues[0];
  M_hat[1] = 1.0 / eigenValues[1];
  M_hat[2] = 1.0 / eigenValues[2];

  if (abs(eigenValues[1] - eigenValues[2]) <= 1e-12) {
    M_hat[3] = 1.0 / eigenValues[1];
  } else {
    M_hat[3] = log(eigenValues[1] / eigenValues[2]) / (eigenValues[1] - eigenValues[2]);
  }

  if (abs(eigenValues[0] - eigenValues[2]) <= 1e-12) {
    M_hat[4] = 1.0 / eigenValues[0];
  } else {
    M_hat[4] = log(eigenValues[0] / eigenValues[2]) / (eigenValues[0] - eigenValues[2]);
  }

  if (abs(eigenValues[0] - eigenValues[1]) <= 1e-12) {
    M_hat[5] = 1.0 / eigenValues[0];
  } else {
    M_hat[5] = log(eigenValues[0] / eigenValues[1]) / (eigenValues[0] - eigenValues[1]);
  }

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      real64 deltaFe[3][3] = {0};
      deltaFe[i][j] = 1.0;

      // delCe = delFe^T * Fe + Fe^T * delFe
      real64 dCe[3][3] = {0};
      for (int l = 0; l < 3; ++l) {
        dCe[j][l] += deltaFe[i][j] * elasticDeformGrad[i][l];
        dCe[l][j] += elasticDeformGrad[i][l] * deltaFe[i][j];
      }

      real64 deltaCe[6] = {};
      LvArray::tensorOps::denseToSymmetric<3>(deltaCe, dCe);

      // material tangent
      real64 deltaP_mat[3][3] = {};
      computeMaterialTangentColumn(deltaCe, M_hat, eigenVectors, eigenVectorsT, smallStrainStiffness, fInvT, deltaP_mat);

      // geometric tangent
      real64 deltaP_geo[3][3] = {};
      real64 deltaFeT_FeInvT[3][3] = {0};
      for (int l = 0; l < 3; ++l) {
        deltaFeT_FeInvT[j][l] += deltaFe[i][j] * fInvT[i][l];
      }
      LvArray::tensorOps::Rij_eq_AikBkj<3, 3, 3>(deltaP_geo, firstPiolaStress, deltaFeT_FeInvT);

      real64 deltaP_total[3][3] = {};
      LvArray::tensorOps::copy<3, 3>(deltaP_total, deltaP_mat);
      LvArray::tensorOps::scaledAdd<3, 3>(deltaP_total, deltaP_geo, -1.0);

      // input results into D operator
      int col = 3 * i + j;
      for (int m = 0; m < 3; ++m) {
        for (int n = 0; n < 3; ++n) {
          // row major flattening
          int row = 3 * m + n;
          stiffness[row][col] = deltaP_total[m][n];
        }
      }

    }
  }

}

GEOS_HOST_DEVICE
inline
void ElasticIsotropicFiniteStrainUpdates::finiteStrainUpdate(localIndex const k, localIndex const q,
                                                             real64 const (&elasticDeformGrad)[3][3],
                                                             real64 (&firstPiolaStress)[3][3],
                                                             real64 (&stiffness)[9][9]) const
{
  real64 kirchhoffStress[6] = {};
  finiteStrainNoStateUpdate(k, q, elasticDeformGrad, firstPiolaStress, kirchhoffStress, stiffness);

  // Can only save the symmetric kirchhoff stress right now
  LvArray::tensorOps::copy<6>(m_oldStress[k][q], m_newStress[k][q]);
  LvArray::tensorOps::copy<6>(m_newStress[k][q], kirchhoffStress);
}

GEOS_HOST_DEVICE
inline
void ElasticIsotropicFiniteStrainUpdates::finiteStrainUpdate(localIndex const k, localIndex const q,
                                                             real64 const (&elasticDeformGrad)[3][3],
                                                             real64 (&firstPiolaStress)[3][3],
                                                             DiscretizationOps & stiffness) const
{
  finiteStrainUpdate(k, q, elasticDeformGrad, firstPiolaStress, stiffness.m_c);
}

/**
 * @class ElasticIsotropicFiniteStrain
 *
 * Finite strain isotropic elastic material model.
 */
class ElasticIsotropicFiniteStrain : public ElasticIsotropic
{
public:
  /// @typedef Alias for ElasticIsotropicFiniteStrainUpdates
  using KernelWrapper = ElasticIsotropicFiniteStrainUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  ElasticIsotropicFiniteStrain( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~ElasticIsotropicFiniteStrain() override;

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ElasticIsotropicFiniteStrain";

  /**
   * @brief Static catalog string
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Create a instantiation of the ElasticFiniteStrainIsotropicUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of ElasticFiniteStrainIsotropicUpdate.
   */
  KernelWrapper createKernelUpdates(bool const includeState = true) const
  {
    if (includeState) {
      return KernelWrapper(m_bulkModulus, m_shearModulus, m_thermalExpansionCoefficient,
                           m_newStress, m_oldStress, m_disableInelasticity);
    } else {
      return KernelWrapper(m_bulkModulus, m_shearModulus, m_thermalExpansionCoefficient,
                           arrayView3d< real64, solid::STRESS_USD >(),
                           arrayView3d< real64, solid::STRESS_USD >(),
                           m_disableInelasticity);
    }
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_bulkModulus,
                          m_shearModulus,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_disableInelasticity );
  }

protected:
  virtual void postInputInitialization() override;
};

} // namespace constitutive

} // namespace geos

#endif /* GEOS_CONSTITUTIVE_SOLID_ELASTICISOTROPIC_HPP_ */
