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

/// @namespace Namespace to collect common property conversion functions (elastic, poroelastic, etc.)
namespace conversions
{

/// @namespace Bulk modulus and Shear modulus as input
namespace BulkModAndShearMod
{

/**
 * @brief Compute Young's modulus
 * @param[in] K Bulk modulus
 * @param[in] G Shear modulus
 * @return Young's modulus
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toYoungsMod( real64 const & K, 
                    real64 const & G )
{
  return 9.0 * K * G / ( 3.0 * K + G );
}
*/
/**
 * @brief Compute Poisson's ratio
 * @param[in] K Bulk modulus
 * @param[in] G Shear modulus
 * @return Poisson's ratio
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toPoissonRatio( real64 const & K, 
                       real64 const & G )
{
  return ( 3.0 * K - 2.0 * G ) / ( 6.0 * K + 2.0 * G );
}
*/
/**
 * @brief Compute First Lamé parameter
 * @param[in] K Bulk modulus
 * @param[in] G Shear modulus
 * @return First Lamé parameter
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toFirstLame( real64 const & K, 
                    real64 const & G )
{
  return K - 2.0 * G / 3.0;
}
*/
} /* namespace BulkModeAndShearMod */

/// @namespace Young's modulus and Poisson's ratio as input
namespace YoungsModAndPoissonRatio
{

/**
 * @brief Compute Bulk modulus
 * @param[in] E Young's modulus
 * @param[in] nu Poisson's ratio
 * @return Bulk modulus
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toBulkMod( real64 const & E, 
                  real64 const & nu )
{
  return E / ( 3.0 - 6.0 * nu );
}
*/
/**
 * @brief Compute Shear modulus
 * @param[in] E Young's modulus
 * @param[in] nu Poisson's ratio
 * @return Shear modulus
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toShearMod( real64 const & E, 
                   real64 const & nu )
{
  return E / ( 2.0 + 2.0 * nu );
}
*/
/**
 * @brief Compute First Lamé parameter
 * @param[in] E Young's modulus
 * @param[in] nu Poisson's ratio
 * @return First Lamé parameter
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toFirstLame( real64 const & E, 
                    real64 const & nu )
{
  return E * nu / ( 1.0 + nu ) / ( 1.0 - 2.0 * nu );
}
*/
} /* namespace YoungsModAndPoissonRatio*/

/// @namespace Shear modulus and Poisson's ratio as input
namespace ShearModAndPoissonRatio
{

/**
 * @brief Compute Bulk modulus
 * @param[in] G Shear modulus
 * @param[in] nu Poisson's ratio
 * @return Bulk modulus
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toBulkMod( real64 const & G, 
                  real64 const & nu )
{
  return G * ( 2.0 + 2.0 * nu ) / ( 3.0 - 6.0 * nu );
}
*/
/**
 * @brief Compute Young's modulus
 * @param[in] G Shear modulus
 * @param[in] nu Poisson's ratio
 * @return Young's modulus
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toYoungsMod( real64 const & G, 
                    real64 const & nu )
{
  return 2.0 * G * ( 1.0 + nu );
}
*/
/**
 * @brief Compute First Lamé parameter
 * @param[in] G Shear modulus
 * @param[in] nu Poisson's ratio
 * @return First Lamé parameter
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toFirstLame( real64 const & G, 
                    real64 const & nu )
{
  return 2.0 * G * nu / ( 1.0 - 2.0 * nu );
}
*/
} /* namespace ShearModAndPoissonRatio*/

/// @namespace Bulk modulus and Poisson's ratio as input
namespace BulkModAndPoissonRatio
{

/**
 * @brief Compute Young's modulus
 * @param[in] K Bulk modulus
 * @param[in] nu Poisson's ratio
 * @return Young's modulus
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toYoungsMod( real64 const & K, 
                    real64 const & nu )
{
  return K * ( 3.0 - 6.0 * nu );
}
*/
/**
 * @brief Compute Shear modulus
 * @param[in] K Bulk modulus
 * @param[in] nu Poisson's ratio
 * @return Shear modulus
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toShearMod( real64 const & K, 
                   real64 const & nu )
{
  return K * ( 3.0 - 6.0 * nu) / ( 2.0 + 2.0 * nu );
}
*/
/**
 * @brief Compute First Lamé parameter
 * @param[in] K Bulk modulus
 * @param[in] nu Poisson's ratio
 * @return First Lamé parameter
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toFirstLame( real64 const & K, 
                    real64 const & nu )
{
  return 3.0 * K * nu / ( 1.0 + nu );
}
*/
} /* namespace BulkModAndPoissonRatio */

/// @namespace Bulk modulus and Young's modulus as input
namespace BulkModAndYoungsMod
{

/**
 * @brief Compute Shear modulus
 * @param[in] K Bulk modulus
 * @param[in] E Young's ratio
 * @return Shear modulus
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toShearMod( real64 const & K, 
                   real64 const & E )
{
  return 3.0 * K * E / ( 9.0 * K - E );
}
*/
/**
 * @brief Compute Poisson's ratio
 * @param[in] K Bulk modulus
 * @param[in] E Young's modulus
 * @return Poisson's ratio
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toPoissonRatio( real64 const & K, 
                       real64 const & E )
{
  return ( 3.0 * K - E ) / ( 6.0 * K); 
}
*/
/**
 * @brief Compute First Lamé parameter
 * @param[in] K Bulk modulus
 * @param[in] E Young's modulus
 * @return First Lamé parameter
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toFirstLame( real64 const & K, 
                    real64 const & E )
{
  return K * ( 9.0 * K - 3.0 * E ) / ( 9.0 * K - E );
}
*/
} /* namespace BulkModAndYoungsMod */

/// @namespace Shear modulus and Young's modulus as input
namespace ShearModAndYoungsMod
{
/**
 * @brief Compute Poisson's ratio
 * @param[in] G Shear modulus
 * @param[in] E Young's modulus
 * @return Poisson's ratio
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toPoissonRatio( real64 const & G, 
                       real64 const & E )
{
  return 0.5 * E / G - 1.0;
}
*/
/**
 * @brief Compute Bulk modulus
 * @param[in] G Shear modulus
 * @param[in] E Young's modulus
 * @return Bulk modulus
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toBulkMod( real64 const & G, 
                  real64 const & E )
{
  return E * G / ( 9.0 * G - 3.0 * E );
}
*/
/**
 * @brief Compute First Lamé parameter
 * @param[in] G Shear modulus
 * @param[in] E Young's modulus
 * @return First Lamé parameter
 */
/**
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toFirstLame( real64 const & G, 
                    real64 const & E )
{
  return ( E - 2.0 * G ) * G / ( 3.0 * G - E );
}
*/
} /* namespace ShearModAndYoungsMod*/


/// @namespace Compute Biot coefficient from different input combinations
namespace BiotCoefficient
{

/**
 * @brief Compute Biot's coefficient
 * @param[in] K drained Bulk modulus
 * @param[in] Ks Bulk modulus of the solid phase
 * @return Biot's coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKKs( real64 const & K, 
               real64 const & Ks )
{
  return 1.0 - K / Ks;
}

/**
 * @brief Compute Biot's coefficient
 * @param[in] Ku undrained Bulk modulus
 * @param[in] K drained Bulk modulus
 * @param[in] M Biot's modulus
 * @return Biot's coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKuKM( real64 const & Ku, 
                real64 const & K, 
                real64 const & M )
{
  return sqrt( ( Ku - K ) / M );
}


/**
 * @brief Compute Biot's coefficient
 * @param[in] Lu undrained Lamé modulus
 * @param[in] L drained Lamé modulus
 * @param[in] M Biot's modulus
 * @return Biot's coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useLuLM( real64 const & Lu, 
                real64 const & L, 
                real64 const & M )
{
  return sqrt( ( Lu - L ) / M );
}

/**
 * @brief Compute Biot's coefficient
 * @param[in] Ku undrained Bulk modulus
 * @param[in] K drained Bulk modulus
 * @param[in] B Skempton's coefficient
 * @return Biot's coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKuKB( real64 const & Ku, 
                real64 const & K, 
                real64 const & B )
{
  return ( Ku - K ) / Ku / B;
}

/**
 * @brief Compute Biot's coefficient
 * @param[in] Ku undrained Bulk modulus
 * @param[in] M Biot's modulus
 * @param[in] B Skempton's coefficient
 * @return Biot's coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKuMB( real64 const & Ku, 
                real64 const & M, 
                real64 const & B )
{
  return Ku * B / M;
}

/**
 * @brief Compute Biot's coefficient
 * @param[in] nuu undrained Poisson's ratio
 * @param[in] nu drained Poisson's ratio
 * @param[in] B Skempton's coefficient
 * @return Biot's coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useNuuNuB( real64 const & nuu, 
                  real64 const & nu, 
                  real64 const & B )
{
  return 3.0 * ( nuu - nu ) / B / ( 1.0 - 2.0 * nu ) / ( 1.0 + nuu );
}

} /* BiotCoefficient */

/// @namespace Compute Biot modulus from different input combinations
namespace BiotModulus
{
/**
 * @brief Compute Biot's modulus
 * @param[in] Ks Bulk modulus of the solid phase
 * @param[in] Kf Bulk modulus of the fluid phase
 * @param[in] alpha Biot's coefficient
 * @param[in] phi porosity
 * @return Biot's modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKsKfBiotCoeffPorosity( real64 const & Ks, 
                                 real64 const & Kf, 
                                 real64 const & alpha, 
                                 real64 const & phi )
{
  return 1.0 / ( phi / Kf + ( alpha - phi ) / Ks );
}

/**
 * @brief Compute Biot's modulus
 * @param[in] Ku undrained Bulk modulus
 * @param[in] K drained Bulk modulus
 * @param[in] alpha Biot's coefficient
 * @return Biot's modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKuKBiotCoeff( real64 const & Ku, 
                        real64 const & K, 
                        real64 const & alpha )
{
  return ( Ku - K ) / alpha / alpha;
}

/**
 * @brief Compute Biot's modulus
 * @param[in] Lu undrained Lamé modulus
 * @param[in] L drained Lamé modulus
 * @param[in] alpha Biot's coefficient
 * @return Biot's modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useLuLBiotCoeff( real64 const & Lu, 
                        real64 const & L, 
                        real64 const & alpha )
{
  return ( Lu - L ) / alpha / alpha;
}

/**
 * @brief Compute Biot's modulus
 * @param[in] Ku undrained Bulk modulus
 * @param[in] B Skempton's coefficient
 * @param[in] alpha Biot's coefficient
 * @return Biot's modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKuBBiotCoeff( real64 const & Ku, 
                        real64 const & B, 
                        real64 const & alpha )
{
  return Ku * B / alpha;
}

/**
 * @brief Compute Biot's modulus
 * @param[in] K drained Bulk modulus
 * @param[in] B Skempton's coefficient
 * @param[in] alpha Biot's coefficient
 * @return Biot's modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKBBiotCoeff( real64 const & K, 
                       real64 const & B, 
                       real64 const & alpha )
{
  return K * B / alpha / ( 1.0 - alpha * B );
}

/**
 * @brief Compute Biot's modulus
 * @param[in] G Shear modulus
 * @param[in] nuu undrained Poisson's ratio
 * @param[in] nu drained Poisson's ratio
 * @param[in] alpha Biot's coefficient
 * @return Biot's modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useGNuuNuBiotCoeff( real64 const & G, 
                           real64 const & nuu, 
                           real64 const & nu, 
                           real64 const & alpha )
{
  return 2.0 * G * ( nuu - nu ) / alpha / alpha / ( 1.0 - 2.0 * nuu ) / ( 1.0 - 2.0 * nu );
}

} /* BiotModulus */

/// @namespace Compute Skempton coefficient from different input combinations
namespace SkemptonCoefficient
{
/**
 * @brief Compute Skempton's coefficient
 * @param[in] Ku undrained Bulk modulus
 * @param[in] M Biot's modulus
 * @param[in] alpha Biot's coefficient
 * @return Skempton's coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKuMBiotCoeff( real64 const & Ku, 
                        real64 const & M, 
                        real64 const & alpha )
{
  return  M * alpha / Ku;
}

/**
 * @brief Compute Skempton's coefficient
 * @param[in] Ku undrained Bulk modulus
 * @param[in] K drained Bulk modulus
 * @param[in] alpha Biot's coefficient
 * @return Skempton's coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKuKBiotCoeff( real64 const & Ku, 
                        real64 const & K, 
                        real64 const & alpha )
{
  return ( Ku - K ) / Ku / alpha;
}

/**
 * @brief Compute Skempton's coefficient
 * @param[in] nuu undrained Poisson's ratio
 * @param[in] nu drained Poisson's ratio
 * @param[in] alpha Biot's coefficient
 * @return Skempton's coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useNuuNuBiotCoeff( real64 const & nuu, 
                          real64 const & nu, 
                          real64 const & alpha )
{
  return 3.0 * ( nuu - nu ) / alpha / ( 1.0 - 2.0 * nu ) / ( 1.0 + nuu );
}

} /* SkemptonCoefficient */

/// @namespace Compute undrained Bulk modulus from different input combinations
namespace UndrainedBulkMod
{
/**
 * @brief Compute undrained Bulk modulus
 * @param[in] K drained Bulk modulus
 * @param[in] M Biot's modulus
 * @param[in] alpha Biot's coefficient
 * @return Undrained Bulk modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKMBiotCoeff( real64 const & K, 
                       real64 const & M, 
                       real64 const & alpha )
{
  return K + M * alpha * alpha;
}

/**
 * @brief Compute undrained Bulk modulus
 * @param[in] B Skempton's coefficient
 * @param[in] M Biot's modulus
 * @param[in] alpha Biot's coefficient
 * @return Undrained Bulk modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useBMBiotCoeff( real64 const & B, 
                       real64 const & M, 
                       real64 const & alpha )
{
  return M * alpha / B;
}

/**
 * @brief Compute undrained Bulk modulus
 * @param[in] K drained Bulk modulus
 * @param[in] B Skempton's coefficient
 * @param[in] alpha Biot's coefficient
 * @return Undrained Bulk modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKBBiotCoeff( real64 const & K, 
                       real64 const & B, 
                       real64 const & alpha )
{
  return K / ( 1.0 -  alpha * B );
}

} /* UndrainedBulkMod */


/// @namespace Compute drained Bulk modulus from different input combinations
namespace DrainedBulkMod
{
/**
 * @brief Compute drained Bulk modulus
 * @param[in] Ku undrained Bulk modulus
 * @param[in] M Biot's modulus
 * @param[in] alpha Biot's coefficient
 * @return Drained Bulk modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKuMBiotCoeff( real64 const & Ku, 
                        real64 const & M, 
                        real64 const & alpha )
{
  return Ku - M * alpha * alpha;
}

/**
 * @brief Compute drained Bulk modulus
 * @param[in] Ku undrained Bulk modulus
 * @param[in] B Skempton's coefficient
 * @param[in] alpha Biot's coefficient
 * @return Drained Bulk modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKuBBiotCoeff( real64 const & Ku, 
                        real64 const & B, 
                        real64 const & alpha )
{
  return Ku * ( 1.0 -  alpha * B );
}

} /* DrainedBulkMod */

/// @namespace Compute diffusion coefficient from different input combinations
namespace DiffusionCoeff
{
/**
 * @brief Compute diffusion coefficient
 * @param[in] Ku undrained Bulk modulus
 * @param[in] K drained Bulk modulus
 * @param[in] G Shear modulus
 * @param[in] M Biot's modulus
 * @param[in] kappa ratio permeability/fluid dynamic viscosity
 * @return Diffusion coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useKuKGMKappa( real64 const & Ku, 
                      real64 const & K, 
                      real64 const & G, 
                      real64 const & M, 
                      real64 const & kappa )
{
  return kappa * M * ( 3.0 * K + 4.0 * G ) / ( 3.0 * Ku + 4.0 * G );
}

/**
 * @brief Compute diffusion coefficient
 * @param[in] nuu undrained Poisson's ratio
 * @param[in] nu drained Poisson's ratio
 * @param[in] G Shear modulus
 * @param[in] alpha Biot's coefficient
 * @param[in] kappa ratio permeability/fluid dynamic viscosity
 * @return Diffusion coefficient
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 useNuuNuGBiotCoeffKappa( real64 const & nuu, 
                                real64 const & nu, 
                                real64 const & G, 
                                real64 const & alpha, 
                                real64 const & kappa )
{
  return 2.0 * kappa * G * ( 1.0 - nu ) * ( nuu - nu ) / alpha / alpha / ( 1.0 - 2.0 * nu ) / ( 1.0 - 2.0 * nu ) / ( 1.0 - nuu );
}

} /* DiffusionCoeff */

} /* namespace conversions */

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_PROPERTYCONVERSIONS_HPP_ */
