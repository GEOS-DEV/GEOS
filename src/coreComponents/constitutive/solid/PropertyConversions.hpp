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
 * @file PropertyConversions.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_PROPERTYCONVERSIONS_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_PROPERTYCONVERSIONS_HPP_

namespace geosx
{

namespace constitutive
{

/**
 * @class PoissonRatio
 */
struct PoissonRatio
{
public:

  /**
   * @brief Compute Poisson's ratio from other elastic parameters
   * @return Poisson's ratio
   */
  real64 getValue() const
  {
    if( m_K > 0 && m_G > 0 )
    {
      return ( 3.0 * m_K - 2.0 * m_G ) / ( 6.0 * m_K + 2.0 * m_G );
    }
    else if( m_K > 0 && m_E > 0 )
    {
      return ( 3.0 * m_K - m_E ) / ( 6.0 * m_K);
    }
    else if( m_G > 0 && m_E > 0 )
    {
      return 0.5 * m_E / m_G - 1.0;
    }
    else
    {
      return -0.5;
      GEOSX_ERROR( "A specific pair of elastic constants is required: (K,G) or (K,E) or (G,E)" );
    }
  }

  PoissonRatio setBulkModulus( real64 const K )
  {
    m_K = K;
    return *this;
  }

  PoissonRatio setShearModulus( real64 const G )
  {
    m_G = G;
    return *this;
  }

  PoissonRatio setYoungModulus( real64 const E )
  {
    m_E = E;
    return *this;
  }

private:

  /// Bulk modulus
  real64 m_K = 0.0;

  /// Shear modulus
  real64 m_G = 0.0;

  /// Young's modulus
  real64 m_E = 0.0;
};

/**
 * @class YoungModulus
 */
struct YoungModulus
{
public:

  /**
   * @brief Compute Young's modulus from other elastic parameters
   * @return Young's modulus
   */
  real64 getValue() const
  {
    if( m_K > 0 && m_G > 0 )
    {
      return 9.0 * m_K * m_G / ( 3.0 * m_K + m_G );
    }
    else if( m_K > 0 && m_nu > -0.5 && m_nu < 0.5 )
    {
      return m_K * ( 3.0 - 6.0 * m_nu );
    }
    else if( m_G > 0 && m_nu > -0.5 && m_nu < 0.5 )
    {
      return 2.0 * m_G * ( 1.0 + m_nu );
    }
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific pair of elastic constants is required: (K,G) or (K,nu) or (G,nu)" );
    }
  }

  YoungModulus setBulkModulus( real64 const K )
  {
    m_K = K;
    return *this;
  }

  YoungModulus setShearModulus( real64 const G )
  {
    m_G = G;
    return *this;
  }

  YoungModulus setPoissonRatio( real64 const nu )
  {
    m_nu = nu;
    return *this;
  }

private:

  /// Bulk modulus
  real64 m_K = 0.0;

  /// Shear modulus
  real64 m_G = 0.0;

  /// Poisson's ratio
  real64 m_nu = -0.5;
};

/**
 * @class BulkModulus
 */
struct BulkModulus
{
public:

  /**
   * @brief Compute Bulk modulus from other elastic parameters
   * @return Bulk modulus
   */
  real64 getValue() const
  {
    if( m_E > 0 && m_G > 0 )
    {
      return m_E * m_G / ( 9.0 * m_G - 3.0 * m_E );
    }
    else if( m_E > 0 && m_nu > -0.5 && m_nu < 0.5 )
    {
      return m_E / ( 3.0 - 6.0 * m_nu );
    }
    else if( m_G > 0 && m_nu > -0.5 && m_nu < 0.5 )
    {
      return m_G * ( 2.0 + 2.0 * m_nu ) / ( 3.0 - 6.0 * m_nu );
    }
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific pair of elastic constants is required: (E,G) or (E,nu) or (G,nu)" );
    }
  }

  BulkModulus setYoungModulus( real64 const E )
  {
    m_E = E;
    return *this;
  }

  BulkModulus setShearModulus( real64 const G )
  {
    m_G = G;
    return *this;
  }

  BulkModulus setPoissonRatio( real64 const nu )
  {
    m_nu = nu;
    return *this;
  }

private:

  /// Young's modulus
  real64 m_E = 0.0;

  /// Shear modulus
  real64 m_G = 0.0;

  /// Poisson's ratio
  real64 m_nu = -0.5;
};

/**
 * @class ShearModulus
 */
struct ShearModulus
{
public:

  /**
   * @brief Compute Shear modulus from other elastic parameters
   * @return Shear modulus
   */
  real64 getValue() const
  {
    if( m_E > 0 && m_K > 0 )
    {
      return 3.0 * m_K * m_E / ( 9.0 * m_K - m_E );
    }
    else if( m_E > 0 && m_nu > -0.5 && m_nu < 0.5 )
    {
      return m_E / ( 2.0 + 2.0 * m_nu );
    }
    else if( m_K > 0 && m_nu > -0.5 && m_nu < 0.5 )
    {
      return m_K * ( 3.0 - 6.0 * m_nu) / ( 2.0 + 2.0 * m_nu );
    }
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific pair of elastic constants is required: (E,K) or (E,nu) or (K,nu)" );
    }
  }

  ShearModulus setYoungModulus( real64 const E )
  {
    m_E = E;
    return *this;
  }

  ShearModulus setBulkModulus( real64 const K )
  {
    m_K = K;
    return *this;
  }

  ShearModulus setPoissonRatio( real64 const nu )
  {
    m_nu = nu;
    return *this;
  }

private:

  /// Young's modulus
  real64 m_E = 0.0;

  /// Bulk modulus
  real64 m_K = 0.0;

  /// Poisson's ratio
  real64 m_nu = -0.5;
};

/**
 * @class LameModulus
 */
struct LameModulus
{
public:

  /**
   * @brief Compute Lame modulus from other elastic parameters
   * @return Lame modulus
   */
  real64 getValue() const
  {
    if( m_K > 0 && m_G > 0 )
    {
      return m_K - 2.0 * m_G / 3.0;
    }
    if( m_K > 0 && m_E > 0 )
    {
      return m_K * ( 9.0 * m_K - 3.0 * m_E ) / ( 9.0 * m_K - m_E );
    }
    if( m_G > 0 && m_E > 0 )
    {
      return ( m_E - 2.0 * m_G ) * m_G / ( 3.0 * m_G - m_E );
    }
    else if( m_K > 0 && m_nu > -0.5 && m_nu < 0.5 )
    {
      return 3.0 * m_K * m_nu / ( 1.0 + m_nu );
    }
    else if( m_G > 0 && m_nu > -0.5 && m_nu < 0.5 )
    {
      return 2.0 * m_G * m_nu / ( 1.0 - 2.0 * m_nu );
    }
    else if( m_E > 0 && m_nu > -0.5 && m_nu < 0.5 )
    {
      return m_E * m_nu / ( 1.0 + m_nu ) / ( 1.0 - 2.0 * m_nu );
    }
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific pair of elastic constants is required: (K,G), (K,E), (G,E), (K,nu), (G,nu) or (E,nu)" );
    }
  }

  LameModulus setBulkModulus( real64 const K )
  {
    m_K = K;
    return *this;
  }

  LameModulus setShearModulus( real64 const G )
  {
    m_G = G;
    return *this;
  }

  LameModulus setYoungModulus( real64 const E )
  {
    m_E = E;
    return *this;
  }

  LameModulus setPoissonRatio( real64 const nu )
  {
    m_nu = nu;
    return *this;
  }

private:

  /// Bulk modulus
  real64 m_K = 0.0;

  /// Shear modulus
  real64 m_G = 0.0;

  /// Young's modulus
  real64 m_E = 0.0;

  /// Poisson's ratio
  real64 m_nu = -0.5;
};

/**
 * @class BiotCoefficient
 */
struct BiotCoefficient
{
public:

  /**
   * @brief Compute Biot's coefficient from other poroelastic parameters
   * @return Biot's coefficient
   */
  real64 getValue() const
  {
    if( m_Kd > 0 && m_Ks > 0 )
    {
      return 1.0 - m_Kd / m_Ks;
    }
    else if( m_Ku > 0 && m_Kd > 0 && m_M > 0 && m_Ku > m_Kd )
    {
      return sqrt( ( m_Ku - m_Kd ) / m_M );
    }
    else if( m_Ku > 0 && m_Kd > 0 && m_B > 0 && m_Ku > m_Kd )
    {
      return ( m_Ku - m_Kd ) / m_Ku / m_B;
    }
    else if( m_Ku > 0 && m_M > 0 && m_B > 0 )
    {
      return m_Ku * m_B / m_M;
    }
    else if( m_nud > -0.5 && m_nud < 0.5 && m_nuu > -0.5 && m_nuu < 0.5 && m_B > 0 && m_nuu > m_nud )
    {
      return 3.0 * ( m_nuu - m_nud ) / m_B / ( 1.0 - 2.0 * m_nud ) / ( 1.0 + m_nuu );
    }
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific combination of poroelastic constants is required" );
    }
  }

  BiotCoefficient setDrainedBulkModulus( real64 const Kd )
  {
    m_Kd = Kd;
    return *this;
  }

  BiotCoefficient setUndrainedBulkModulus( real64 const Ku )
  {
    m_Ku = Ku;
    return *this;
  }

  BiotCoefficient setDrainedPoissonRatio( real64 const nud )
  {
    m_nud = nud;
    return *this;
  }

  BiotCoefficient setUndrainedPoissonRatio( real64 const nuu )
  {
    m_nuu = nuu;
    return *this;
  }

  BiotCoefficient setSolidBulkModulus( real64 const Ks )
  {
    m_Ks = Ks;
    return *this;
  }

  BiotCoefficient setBiotModulus( real64 const M )
  {
    m_M = M;
    return *this;
  }

  BiotCoefficient setSkemptonCoefficient( real64 const B )
  {
    m_B = B;
    return *this;
  }

private:

  /// Drained Bulk modulus
  real64 m_Kd = 0.0;

  /// Undrained Bulk modulus
  real64 m_Ku = 0.0;

  /// Drained Poisson's ratio
  real64 m_nud = -0.5;

  /// Undrained Poisson's ratio
  real64 m_nuu = -0.5;

  /// Bulk modulus of the solid phase
  real64 m_Ks = 0.0;

  /// Biot's modulus
  real64 m_M = 0.0;

  /// Skempton's coefficient
  real64 m_B = 0.0;

};


/**
 * @class BiotModulus
 */
struct BiotModulus
{
public:

  /**
   * @brief Compute Biot's modulus from other poroelastic parameters
   * @return Biot's modulus
   */
  real64 getValue() const
  {
    if( m_phi > 0 && m_Kf > 0  && m_Ks > 0 && m_alpha > 0 )
    {
      return 1.0 / ( m_phi / m_Kf + ( m_alpha - m_phi ) / m_Ks );
    }
    else if( m_Ku > 0 && m_Kd > 0 && m_alpha > 0 && m_Ku > m_Kd )
    {
      return ( m_Ku - m_Kd ) / m_alpha / m_alpha;
    }
    else if( m_Ku > 0 && m_B > 0 && m_alpha > 0 )
    {
      return m_Ku * m_B / m_alpha;
    }
    else if( m_Kd > 0 && m_B > 0 && m_alpha > 0 && m_alpha * m_B < 1 )
    {
      return m_Kd * m_B / m_alpha / ( 1.0 - m_alpha * m_B );
    }
    else if( m_G > 0 && m_nuu > -0.5 && m_nuu < 0.5 && m_nud > -0.5 && m_nud < 0.5 && m_alpha > 0 && m_nuu > m_nud )
    {
      return 2.0 * m_G * ( m_nuu - m_nud ) / m_alpha / m_alpha / ( 1.0 - 2.0 * m_nuu ) / ( 1.0 - 2.0 * m_nud );
    }
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific combination of poroelastic constants is required" );
    }
  }


  BiotModulus setSolidBulkModulus( real64 const Ks )
  {
    m_Ks = Ks;
    return *this;
  }

  BiotModulus setFluidBulkModulus( real64 const Kf )
  {
    m_Kf = Kf;
    return *this;
  }

  BiotModulus setPorosity( real64 const phi )
  {
    m_phi = phi;
    return *this;
  }

  BiotModulus setBiotCoefficient( real64 const alpha )
  {
    m_alpha = alpha;
    return *this;
  }

  BiotModulus setDrainedBulkModulus( real64 const Kd )
  {
    m_Kd = Kd;
    return *this;
  }

  BiotModulus setUndrainedBulkModulus( real64 const Ku )
  {
    m_Ku = Ku;
    return *this;
  }

  BiotModulus setShearModulus( real64 const G )
  {
    m_G = G;
    return *this;
  }

  BiotModulus setDrainedPoissonRatio( real64 const nud )
  {
    m_nud = nud;
    return *this;
  }

  BiotModulus setUndrainedPoissonRatio( real64 const nuu )
  {
    m_nuu = nuu;
    return *this;
  }

  BiotModulus setSkemptonCoefficient( real64 const B )
  {
    m_B = B;
    return *this;
  }

private:

  /// Bulk modulus of the solid phase
  real64 m_Ks = 0.0;

  /// Bulk modulus of the fluid phase
  real64 m_Kf = 0.0;

  /// Porosity
  real64 m_phi = 0.0;

  /// Biot's coefficient
  real64 m_alpha = 0.0;

  /// Drained Bulk modulus
  real64 m_Kd = 0.0;

  /// Undrained Bulk modulus
  real64 m_Ku = 0.0;

  /// Shear modulus
  real64 m_G = 0.0;

  /// Drained Poisson's ratio
  real64 m_nud = -0.5;

  /// Undrained Poisson's ratio
  real64 m_nuu = -0.5;

  /// Skempton's coefficient
  real64 m_B = 0.0;

};

/**
 * @class SkemptonCoefficient
 */
struct SkemptonCoefficient
{
public:

  /**
   * @brief Compute Skempton's coefficient from other poroelastic parameters
   * @return Skempton's coefficient
   */
  real64 getValue() const
  {
    if( m_M > 0 && m_alpha > 0 && m_Ku > 0 )
    {
      return m_M * m_alpha / m_Ku;
    }
    else if( m_Ku > 0 && m_Kd > 0 && m_alpha > 0 && m_Ku > m_Kd )
    {
      return ( m_Ku - m_Kd ) / m_Ku / m_alpha;
    }
    else if( m_nud > -0.5 && m_nud < 0.5 && m_nuu > -0.5 && m_nuu < 0.5 && m_alpha > 0 && m_nuu > m_nud )
    {
      return 3.0 * ( m_nuu - m_nud ) / m_alpha / ( 1.0 - 2.0 * m_nud ) / ( 1.0 + m_nuu );
    }
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific combination of poroelastic constants is required" );
    }
  }

  SkemptonCoefficient setDrainedBulkModulus( real64 const Kd )
  {
    m_Kd = Kd;
    return *this;
  }

  SkemptonCoefficient setUndrainedBulkModulus( real64 const Ku )
  {
    m_Ku = Ku;
    return *this;
  }

  SkemptonCoefficient setDrainedPoissonRatio( real64 const nud )
  {
    m_nud = nud;
    return *this;
  }

  SkemptonCoefficient setUndrainedPoissonRatio( real64 const nuu )
  {
    m_nuu = nuu;
    return *this;
  }

  SkemptonCoefficient setBiotCoefficient( real64 const alpha )
  {
    m_alpha = alpha;
    return *this;
  }

  SkemptonCoefficient setBiotModulus( real64 const M )
  {
    m_M = M;
    return *this;
  }

private:

  /// Drained Bulk modulus
  real64 m_Kd = 0.0;

  /// Undrained Bulk modulus
  real64 m_Ku = 0.0;

  /// Drained Poisson's ratio
  real64 m_nud = -0.5;

  /// Undrained Poisson's ratio
  real64 m_nuu = -0.5;

  /// Biot's coefficient
  real64 m_alpha = 0.0;

  /// Biot's modulus
  real64 m_M = 0.0;

};


/**
 * @class UndrainedBulkModulus
 */
struct UndrainedBulkModulus
{
public:

  /**
   * @brief Compute undrained Bulk modulus from other poroelastic parameters
   * @return undrained Bulk modulus
   */
  real64 getValue() const
  {
    if( m_Kd > 0 && m_M > 0 && m_alpha > 0 )
    {
      return m_Kd + m_M * m_alpha * m_alpha;
    }
    else if( m_M > 0 && m_alpha > 0 && m_B > 0 )
    {
      return m_M * m_alpha / m_B;
    }
    else if( m_Kd > 0 && m_B > 0 && m_alpha > 0 && m_alpha * m_B < 1 )
    {
      return m_Kd / ( 1.0 -  m_alpha * m_B );
    }
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific combination of poroelastic constants is required" );
    }
  }

  UndrainedBulkModulus setBiotCoefficient( real64 const alpha )
  {
    m_alpha = alpha;
    return *this;
  }

  UndrainedBulkModulus setBiotModulus( real64 const M )
  {
    m_M = M;
    return *this;
  }

  UndrainedBulkModulus setDrainedBulkModulus( real64 const Kd )
  {
    m_Kd = Kd;
    return *this;
  }

  UndrainedBulkModulus setSkemptonCoefficient( real64 const B )
  {
    m_B = B;
    return *this;
  }

private:

  /// Biot's coefficient
  real64 m_alpha = 0.0;

  /// Biot's modulus
  real64 m_M = 0.0;

  /// Drained Bulk modulus
  real64 m_Kd = 0.0;

  /// Skempton's coefficient
  real64 m_B = 0.0;

};

/**
 * @class DrainedBulkModulus
 */
struct DrainedBulkModulus
{
public:

  /**
   * @brief Compute drained Bulk modulus from other poroelastic parameters
   * @return drained Bulk modulus
   */
  real64 getValue() const
  {
    if( m_Ku > 0 && m_M > 0 && m_alpha > 0 )
    {
      return m_Ku - m_M * m_alpha * m_alpha;
    }
    else if( m_Ku > 0 && m_B > 0 && m_alpha > 0 && m_alpha * m_B < 1 )
    {
      return m_Ku * ( 1.0 - m_alpha * m_B );
    }
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific combination of poroelastic constants is required" );
    }
  }

  DrainedBulkModulus setBiotCoefficient( real64 const alpha )
  {
    m_alpha = alpha;
    return *this;
  }

  DrainedBulkModulus setBiotModulus( real64 const M )
  {
    m_M = M;
    return *this;
  }

  DrainedBulkModulus setUndrainedBulkModulus( real64 const Ku )
  {
    m_Ku = Ku;
    return *this;
  }

  DrainedBulkModulus setSkemptonCoefficient( real64 const B )
  {
    m_B = B;
    return *this;
  }

private:

  /// Biot's coefficient
  real64 m_alpha = 0.0;

  /// Biot's modulus
  real64 m_M = 0.0;

  /// Undrained Bulk modulus
  real64 m_Ku = 0.0;

  /// Skempton's coefficient
  real64 m_B = 0.0;

};

/**
 * @class DrainedVolumetricTEC
 */
struct DrainedVolumetricTEC
{
public:

  /**
   * @brief Compute drained volumetric Thermal Expansion Coefficient (TEC) 
   * from other thermoporoelastic parameters
   * @return drained volumetric TEC
   */
  real64 getValue() const
  {
    if( m_TECs > 0 && m_TECp > 0 && m_phi > 0 )
    {
      return m_TECs + m_TECp / ( 1.0 - m_phi );
    }
    else if( m_TECs > 0 )
    {
      return m_TECs; ///< Ideal porous medium
    }    
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific combination of thermoporoelastic constants is required" );
    }
  }

  DrainedVolumetricTEC setPorosity( real64 const phi )
  {
    m_phi = phi;
    return *this;
  }

  DrainedVolumetricTEC setSolidVolumetricTEC( real64 const TECs )
  {
    m_TECs = TECs;
    return *this;
  }

  DrainedVolumetricTEC setPorosityTEC( real64 const TECp )
  {
    m_TECp = TECp;
    return *this;
  }
private:

  /// Porosity
  real64 m_phi = 0.0;

  /// Volumetric TEC of the solid phase
  real64 m_TECs = 0.0;

  /// TEC of the porosity (porosity change vs temperature change)
  real64 m_TECp = 0.0; 

};

/**
 * @class FreeStressFluidExchangeTEC, 
 * This parameter quantifies fluid exchange (expelled from the porous medium) 
 * due to temperature change in drained free stress condition
 */
struct FreeStressFluidExchangeTEC
{
public:

  /**
   * @brief Compute the fluid exchange Thermal Expansion Coefficient 
   * from other thermoporoelastic parameters
   * in drained free stress
   * @return fluid exchange TEC
   */
  real64 getValue() const
  {
    if( m_TECs > 0 && m_TECf > 0 && m_TECp > 0 && m_phi > 0 )
    {
      return m_phi *( m_TECf - m_TECs ) - m_TECp / ( 1.0 - m_phi );
    }
    else if( m_TECs > 0 && m_TECf > 0 )
    {
      return m_phi * ( m_TECf - m_TECs ); ///< Ideal porous medium
    }    
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific combination of thermoporoelastic constants is required" );
    }
  }

  FreeStressFluidExchangeTEC setPorosity( real64 const phi )
  {
    m_phi = phi;
    return *this;
  }

  FreeStressFluidExchangeTEC setSolidVolumetricTEC( real64 const TECs )
  {
    m_TECs = TECs;
    return *this;
  }

  FreeStressFluidExchangeTEC setFluidVolumetricTEC( real64 const TECf )
  {
    m_TECf = TECf;
    return *this;
  }

  FreeStressFluidExchangeTEC setPorosityTEC( real64 const TECp )
  {
    m_TECp = TECp;
    return *this;
  }
private:

  /// Porosity
  real64 m_phi = 0.0;

  /// Volumetric TEC of the solid phase
  real64 m_TECs = 0.0;

  /// Volumetric TEC of the fluid phase
  real64 m_TECf = 0.0;

  /// TEC of the porosity (porosity change vs temperature change)
  real64 m_TECp = 0.0; 

};

/**
 * @class IsochoreFluidExchangeTEC, 
 * This parameter quantifies fluid exchange (expelled from the porous medium) 
 * due to temperature change in drained isochore condition
 */
struct IsochoreFluidExchangeTEC
{
public:

  /**
   * @brief Compute the fluid exchange Thermal Expansion Coefficient 
   * from other thermoporoelastic parameters
   * in drained isochore condition
   * @return fluid exchange TEC
   */
  real64 getValue() const
  {
    if( m_TECd > 0 && m_TECv > 0 && m_alpha > 0 )
    {
      return m_alpha * m_TECd + m_TECv;
    }   
    else if( m_TECs > 0 && m_TECf > 0 && m_TECp > 0 && m_phi > 0 && m_alpha > 0 )
    {
      return ( m_alpha - m_phi ) * m_TECs + m_phi * m_TECf - ( 1.0 - m_alpha ) * m_TECp / ( 1.0 - m_phi );
    }
    else if( m_TECs > 0 && m_TECf > 0 && m_phi > 0 && m_alpha > 0 )
    {
      return ( m_alpha - m_phi ) * m_TECs + m_phi * m_TECf; ///< Ideal porous medium
    }
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific combination of thermoporoelastic constants is required" );
    }
  }
  IsochoreFluidExchangeTEC setDrainedVolumetricTEC( real64 const TECd )
  {
    m_TECd = TECd;
    return *this;
  }

  IsochoreFluidExchangeTEC setFluidExchangeTEC( real64 const TECv )
  {
    m_TECv = TECv;
    return *this;
  }

  IsochoreFluidExchangeTEC setBiotCoefficient( real64 const alpha )
  {
    m_alpha = alpha;
    return *this;
  }

  IsochoreFluidExchangeTEC setPorosity( real64 const phi )
  {
    m_phi = phi;
    return *this;
  }

  IsochoreFluidExchangeTEC setSolidVolumetricTEC( real64 const TECs )
  {
    m_TECs = TECs;
    return *this;
  }

  IsochoreFluidExchangeTEC setFluidVolumetricTEC( real64 const TECf )
  {
    m_TECf = TECf;
    return *this;
  }

  IsochoreFluidExchangeTEC setPorosityTEC( real64 const TECp )
  {
    m_TECp = TECp;
    return *this;
  }

private:

  /// Volumetric drained TEC
  real64 m_TECd = 0.0;

  /// Fluid exchange (expelled) TEC in drained free stress condition
  real64 m_TECv = 0.0; 

  /// Biot's coefficient
  real64 m_alpha = 0.0;

  /// Porosity
  real64 m_phi = 0.0;

  /// Volumetric TEC of the solid phase
  real64 m_TECs = 0.0;

  /// Volumetric TEC of the fluid phase
  real64 m_TECf = 0.0;

  /// TEC of the porosity (porosity change vs temperature change)
  real64 m_TECp = 0.0; 
};

/**
 * @class UndrainedVolumetricTEC
 */
struct UndrainedVolumetricTEC
{
public:

  /**
   * @brief Compute undrained volumetric Thermal Expansion Coefficient (TEC) 
   * from other thermoporoelastic parameters
   * @return undrained volumetric TEC
   */
  real64 getValue() const
  {
    if( m_TECd > 0 && m_B > 0 && m_TECv > 0 )
    {
      return m_TECd + m_B * m_TECv;
    }
    else if( m_TECs > 0 && m_TECf > 0 && m_B > 0 && m_phi > 0 && m_TECp > 0 )
    {
      return ( 1.0 - m_phi * m_B ) * m_TECs + m_phi * m_B * m_TECf + ( 1.0 - m_B ) * m_TECp / ( 1.0 - m_phi );
    } 
    else if( m_TECs > 0 && m_TECf > 0 && m_B > 0 && m_phi > 0 )
    {
      return ( 1.0 - m_phi * m_B ) * m_TECs + m_phi * m_B * m_TECf; ///< Ideal porous medium
    }   
    else if( m_TECs > 0 && m_TECf > 0 && m_phi > 0 )
    {
      return ( 1.0 - m_phi ) * m_TECs + m_phi * m_TECf; ///< Incompressible solid and fluid
    } 
    else
    {
      return 0.0;
      GEOSX_ERROR( "A specific combination of thermoporoelastic constants is required" );
    }
  }


  UndrainedVolumetricTEC setDrainedVolumetricTEC( real64 const TECd )
  {
    m_TECd = TECd;
    return *this;
  }

  UndrainedVolumetricTEC setFluidExchangeTEC( real64 const TECv )
  {
    m_TECv = TECv;
    return *this;
  }

  UndrainedVolumetricTEC setSkemptonCoefficient( real64 const B )
  {
    m_B = B;
    return *this;
  }

  UndrainedVolumetricTEC setPorosity( real64 const phi )
  {
    m_phi = phi;
    return *this;
  }

  UndrainedVolumetricTEC setSolidVolumetricTEC( real64 const TECs )
  {
    m_TECs = TECs;
    return *this;
  }

  UndrainedVolumetricTEC setFluidVolumetricTEC( real64 const TECf )
  {
    m_TECf = TECf;
    return *this;
  }

  UndrainedVolumetricTEC setPorosityTEC( real64 const TECp )
  {
    m_TECp = TECp;
    return *this;
  }

private:

  /// Volumetric drained TEC
  real64 m_TECd = 0.0;

  /// Fluid exchange (expelled) TEC in drained free stress condition
  real64 m_TECv = 0.0; 

  /// Skempton's coefficient
  real64 m_B = 0.0;

  /// Porosity
  real64 m_phi = 0.0;

  /// Volumetric TEC of the solid phase
  real64 m_TECs = 0.0;

  /// Volumetric TEC of the fluid phase
  real64 m_TECf = 0.0;

  /// TEC of the porosity (porosity change vs temperature change)
  real64 m_TECp = 0.0; 

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_PROPERTYCONVERSIONS_HPP_ */
