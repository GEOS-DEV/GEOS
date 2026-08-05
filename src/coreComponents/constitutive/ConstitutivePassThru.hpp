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
 * @file ConstitutivePassThru.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_
#define GEOS_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_

#include "ConstitutivePassThruHandler.hpp"
#include "NullModel.hpp"
#include "solid/Damage.hpp"
#include "solid/DamageVolDev.hpp"
#include "solid/DamageSpectral.hpp"
#include "solid/DruckerPrager.hpp"
#include "solid/DruckerPragerExtended.hpp"
#include "solid/PerfectlyPlastic.hpp"
#include "solid/ModifiedCamClay.hpp"
#include "solid/DelftEgg.hpp"
#include "solid/DuvautLionsSolid.hpp"
#include "solid/ElasticIsotropic.hpp"
#include "solid/ElasticIsotropicPressureDependent.hpp"
#include "solid/ElasticTransverseIsotropic.hpp"
#include "solid/ElasticOrthotropic.hpp"
#include "solid/PorousSolid.hpp"
#include "solid/PorousDamageSolid.hpp"
#include "solid/CompressibleSolid.hpp"
#include "solid/ProppantSolid.hpp"
#include "solid/CeramicDamage.hpp"
#include "solid/ReactiveSolid.hpp"
#include "solid/porosity/PressurePorosity.hpp"
#include "solid/porosity/ProppantPorosity.hpp"
#include "solid/porosity/ReactivePorosity.hpp"
#include "surfaceArea/ConstantSurfaceArea.hpp"
#include "surfaceArea/PowerLawSurfaceArea.hpp"
#include "surfaceArea/SubstrateCoverageSurfaceArea.hpp"
#include "permeability/ConstantPermeability.hpp"
#include "permeability/CarmanKozenyPermeability.hpp"
#include "permeability/ExponentialDecayPermeability.hpp"
#include "permeability/ParallelPlatesPermeability.hpp"
#include "permeability/PowerLawPermeability.hpp"
#include "permeability/PressurePermeability.hpp"
#include "permeability/ProppantPermeability.hpp"
#include "permeability/SlipDependentPermeability.hpp"
#include "permeability/WillisRichardsPermeability.hpp"
#include "contact/CoulombFriction.hpp"
#include "contact/RateAndStateFriction.hpp"


namespace geos
{
namespace constitutive
{

/**
 * @struct ConstitutivePassThru
 * @brief Struct to facilitate launching of lambda functions with a compile
 *   time knowledge of what constitutive model is used.
 *
 * This struct works by implementing an if-else or switch-case block for a
 * specific constitutive base type, and executing the lambda passing it a
 * casted pointer to the constitutive relation.
 */
template< typename BASETYPE >
struct ConstitutivePassThru;

/**
 * Specialization for models that derive from ElasticIsotropic.
 */
template<>
struct ConstitutivePassThru< ElasticIsotropic >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< ElasticIsotropic >::execute( constitutiveRelation,
                                                              std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for models that derive from FrictionBase.
 */
template<>
struct ConstitutivePassThru< FrictionBase >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CoulombFriction,
                                 RateAndStateFriction< std::integral_constant< bool, true > >,
                                 RateAndStateFriction< std::integral_constant< bool, false > > >::execute( constitutiveRelation,
                                                                                                           std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for models that derive from CoulombFriction.
 */
template<>
struct ConstitutivePassThru< CoulombFriction >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CoulombFriction >::execute( constitutiveRelation,
                                                             std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for models that derive from CoulombFriction.
 */
template<>
struct ConstitutivePassThru< RateAndStateFrictionBase >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< RateAndStateFriction< std::integral_constant< bool, true > >,
                                 RateAndStateFriction< std::integral_constant< bool, false > > >::execute( constitutiveRelation,
                                                                                                           std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for models that derive from SolidBase.
 */
template<>
struct ConstitutivePassThru< SolidBase >
{

  // NOTE: The switch order here can be fragile if a model derives from another
  //       model, as the dynamic_cast will also cast to a base version.
  //       Models should be ordered such that children come before parents.
  //       For example, DruckerPrager before ElasticIsotropic, DamageVolDev before
  //       Damage, etc.

  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< DamageSpectral< ElasticIsotropic >,
                                 DamageVolDev< ElasticIsotropic >,
                                 Damage< ElasticIsotropic >,
                                 DuvautLionsSolid< DruckerPrager >,
                                 DuvautLionsSolid< DruckerPragerExtended >,
                                 DuvautLionsSolid< ModifiedCamClay >,
                                 DruckerPragerExtended,
                                 ModifiedCamClay,
                                 DelftEgg,
                                 DruckerPrager,
                                 ElasticIsotropic,
                                 ElasticTransverseIsotropic,
                                 ElasticIsotropicPressureDependent,
                                 ElasticOrthotropic >::execute( constitutiveRelation,
                                                                std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * @struct ConstitutivePassThruMPM
 */
template< typename BASETYPE >
struct ConstitutivePassThruMPM;

/**
 * Specialization for models that derive from SolidBase that are used by the MPM solver.
 * NOTE: this is only a temporary dispatch to reduce the compilation time.
 */
template<>
struct ConstitutivePassThruMPM< SolidBase >
{

  // NOTE: The switch order here can be fragile if a model derives from another
  //       model, as the dynamic_cast will also cast to a base version.
  //       Models should be ordered such that children come before parents.
  //       For example, DruckerPrager before ElasticIsotropic, DamageVolDev before
  //       Damage, etc.

  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CeramicDamage,
                                 PerfectlyPlastic,
                                 ElasticIsotropic >::execute( constitutiveRelation,
                                                              std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * @struct ConstitutivePassThruTriaxialDriver
 */
template< typename BASETYPE >
struct ConstitutivePassThruTriaxialDriver;

/**
 * Specialization for models that derive from SolidBase.
 * NOTE: this is only a temporary dispatch to reduce the compilation time.
 */
template<>
struct ConstitutivePassThruTriaxialDriver< SolidBase >
{

  // NOTE: The switch order here can be fragile if a model derives from another
  //       model, as the dynamic_cast will also cast to a base version.
  //       Models should be ordered such that children come before parents.
  //       For example, DruckerPrager before ElasticIsotropic, DamageVolDev before
  //       Damage, etc.

  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< DamageSpectral< ElasticIsotropic >,
                                 DamageVolDev< ElasticIsotropic >,
                                 Damage< ElasticIsotropic >,
                                 DuvautLionsSolid< DruckerPrager >,
                                 DuvautLionsSolid< DruckerPragerExtended >,
                                 DuvautLionsSolid< ModifiedCamClay >,
                                 DruckerPragerExtended,
                                 ModifiedCamClay,
                                 DelftEgg,
                                 DruckerPrager,
                                 ElasticIsotropic,
                                 ElasticTransverseIsotropic,
                                 ElasticIsotropicPressureDependent,
                                 ElasticOrthotropic >::execute( constitutiveRelation,
                                                                std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for the NullModel.
 */
template<>
struct ConstitutivePassThru< NullModel >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    if( auto * const ptr = dynamic_cast< NullModel * >( &constitutiveRelation ) )
    {
      lambda( *ptr );
    }
    else
    {
      GEOS_ERROR( GEOS_FMT( "ConstitutivePassThru< NullModel >::execute failed on constitutive relation {}",
                            LvArray::system::demangleType( constitutiveRelation ) ),
                  constitutiveRelation.getDataContext() );
    }
  }
};

/**
 * Specialization for the PorousSolid< ElasticIsotropic, ConstantPermeability > model.
 */
template<>
struct ConstitutivePassThru< PorousSolid< ElasticIsotropic, ConstantPermeability > >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    if( auto * const ptr = dynamic_cast< PorousSolid< ElasticIsotropic, ConstantPermeability > * >( &constitutiveRelation ) )
    {
      lambda( *ptr );
    }
    else
    {
      GEOS_ERROR( GEOS_FMT( "ConstitutivePassThru< PorousSolid< ElasticIsotropic, ConstantPermeability > >::execute "
                            "failed on constitutive relation {}",
                            LvArray::system::demangleType( constitutiveRelation ) ),
                  constitutiveRelation.getDataContext() );
    }
  }
};

/**
 * Specialization for the Damage models.
 */
template<>
struct ConstitutivePassThru< DamageBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation,
                       LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< DamageSpectral< ElasticIsotropic >,
                                 DamageVolDev< ElasticIsotropic >,
                                 Damage< ElasticIsotropic > >::execute( constitutiveRelation,
                                                                        std::forward< LAMBDA >( lambda ) );
  }
};



/**
 * Specialization for the PorousSolid models.
 */
template<>
struct ConstitutivePassThru< PorousSolidBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< PorousSolid< DruckerPragerExtended, ConstantPermeability >,
                                 PorousSolid< ModifiedCamClay, ConstantPermeability >,
                                 PorousSolid< DelftEgg, ConstantPermeability >,
                                 PorousSolid< DruckerPrager, ConstantPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPrager >, ConstantPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPragerExtended >, ConstantPermeability >,
                                 PorousSolid< DuvautLionsSolid< ModifiedCamClay >, ConstantPermeability >,
                                 PorousSolid< ElasticIsotropic, ConstantPermeability >,
                                 PorousSolid< ElasticTransverseIsotropic, ConstantPermeability >,
                                 PorousSolid< ElasticIsotropicPressureDependent, ConstantPermeability >,
                                 PorousSolid< ElasticOrthotropic, ConstantPermeability >,
                                 PorousSolid< DruckerPragerExtended, CarmanKozenyPermeability >,
                                 PorousSolid< ModifiedCamClay, CarmanKozenyPermeability >,
                                 PorousSolid< DelftEgg, CarmanKozenyPermeability >,
                                 PorousSolid< DruckerPrager, CarmanKozenyPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPrager >, CarmanKozenyPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPragerExtended >, CarmanKozenyPermeability >,
                                 PorousSolid< DuvautLionsSolid< ModifiedCamClay >, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticIsotropic, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticTransverseIsotropic, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticIsotropicPressureDependent, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticOrthotropic, CarmanKozenyPermeability > >::execute( constitutiveRelation,
                                                                                                         std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for the PorousDamageSolid models.
 */
template<>
struct ConstitutivePassThru< PorousDamageSolidBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< PorousDamageSolid< DamageSpectral< ElasticIsotropic > >,
                                 PorousDamageSolid< DamageVolDev< ElasticIsotropic > >,
                                 PorousDamageSolid< Damage< ElasticIsotropic > > >::execute( constitutiveRelation,
                                                                                             std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for the CompressibleSolid models.
 */
template<>
struct ConstitutivePassThru< CompressibleSolidBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CompressibleSolid< PressurePorosity, ConstantPermeability >,
                                 CompressibleSolid< PressurePorosity, CarmanKozenyPermeability >,
                                 CompressibleSolid< PressurePorosity, ExponentialDecayPermeability >,
                                 CompressibleSolid< PressurePorosity, ParallelPlatesPermeability >,
                                 CompressibleSolid< PressurePorosity, PressurePermeability >,
                                 CompressibleSolid< PressurePorosity, SlipDependentPermeability >,
                                 CompressibleSolid< PressurePorosity, WillisRichardsPermeability >
                                 >::execute( constitutiveRelation,
                                             std::forward< LAMBDA >( lambda ) );
  }

  template< typename LAMBDA >
  static void execute( ConstitutiveBase const & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CompressibleSolid< PressurePorosity, ConstantPermeability >,
                                 CompressibleSolid< PressurePorosity, CarmanKozenyPermeability >,
                                 CompressibleSolid< PressurePorosity, ExponentialDecayPermeability >,
                                 CompressibleSolid< PressurePorosity, ParallelPlatesPermeability >,
                                 CompressibleSolid< PressurePorosity, PressurePermeability >,
                                 CompressibleSolid< PressurePorosity, SlipDependentPermeability >,
                                 CompressibleSolid< PressurePorosity, WillisRichardsPermeability >
                                 >::execute( constitutiveRelation,
                                             std::forward< LAMBDA >( lambda ) );
  }
};

/**
 * Specialization for the ReactiveSolid models.
 */
template<>
struct ConstitutivePassThru< ReactiveSolidBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< ReactiveSolid< ReactivePorosity, ConstantPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, ConstantPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, ConstantPermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, SubstrateCoverageSurfaceArea >
                                 >::execute( constitutiveRelation,
                                             std::forward< LAMBDA >( lambda ) );
  }

  template< typename LAMBDA >
  static void execute( ConstitutiveBase const & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< ReactiveSolid< ReactivePorosity, ConstantPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, ConstantPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, ConstantPermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, SubstrateCoverageSurfaceArea >
                                 >::execute( constitutiveRelation,
                                             std::forward< LAMBDA >( lambda ) );
  }
};


/**
 * Specialization for the ProppantModel.
 */
template<>
struct ConstitutivePassThru< ProppantSolid< ProppantPorosity, ProppantPermeability > >
{
  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    if( auto * const ptr = dynamic_cast< ProppantSolid< ProppantPorosity, ProppantPermeability > * >( &constitutiveRelation ) )
    {
      lambda( *ptr );
    }
    else
    {
      GEOS_ERROR( GEOS_FMT( "ConstitutivePassThru< ProppantSolid >::execute failed on constitutive relation {}",
                            LvArray::system::demangleType( constitutiveRelation ) ),
                  constitutiveRelation.getDataContext() );
    }
  }
};


/**
 * Specialization for all CoupledSolid models.
 */
template<>
struct ConstitutivePassThru< CoupledSolidBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CompressibleSolid< PressurePorosity, ConstantPermeability >,
                                 CompressibleSolid< PressurePorosity, CarmanKozenyPermeability >,
                                 CompressibleSolid< PressurePorosity, ExponentialDecayPermeability >,
                                 CompressibleSolid< PressurePorosity, ParallelPlatesPermeability >,
                                 CompressibleSolid< PressurePorosity, PressurePermeability >,
                                 CompressibleSolid< PressurePorosity, SlipDependentPermeability >,
                                 CompressibleSolid< PressurePorosity, WillisRichardsPermeability >,
                                 PorousSolid< DruckerPragerExtended, ConstantPermeability >,
                                 PorousSolid< ModifiedCamClay, ConstantPermeability >,
                                 PorousSolid< DelftEgg, ConstantPermeability >,
                                 PorousSolid< DruckerPrager, ConstantPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPrager >, ConstantPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPragerExtended >, ConstantPermeability >,
                                 PorousSolid< DuvautLionsSolid< ModifiedCamClay >, ConstantPermeability >,
                                 PorousSolid< ElasticIsotropic, ConstantPermeability >,
                                 PorousSolid< ElasticTransverseIsotropic, ConstantPermeability >,
                                 PorousSolid< ElasticIsotropicPressureDependent, ConstantPermeability >,
                                 PorousSolid< ElasticOrthotropic, ConstantPermeability >,
                                 PorousSolid< DruckerPragerExtended, CarmanKozenyPermeability >,
                                 PorousSolid< ModifiedCamClay, CarmanKozenyPermeability >,
                                 PorousSolid< DelftEgg, CarmanKozenyPermeability >,
                                 PorousSolid< DruckerPrager, CarmanKozenyPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPrager >, CarmanKozenyPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPragerExtended >, CarmanKozenyPermeability >,
                                 PorousSolid< DuvautLionsSolid< ModifiedCamClay >, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticIsotropic, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticTransverseIsotropic, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticIsotropicPressureDependent, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticOrthotropic, CarmanKozenyPermeability >,
                                 PorousDamageSolid< DamageSpectral< ElasticIsotropic > >,
                                 PorousDamageSolid< DamageVolDev< ElasticIsotropic > >,
                                 PorousDamageSolid< Damage< ElasticIsotropic > >,
                                 ReactiveSolid< ReactivePorosity, ConstantPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, ConstantPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, ConstantPermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, SubstrateCoverageSurfaceArea > >::execute( constitutiveRelation,
                                                                                                                                   std::forward< LAMBDA >( lambda ) );
  }

  template< typename LAMBDA >
  static void execute( ConstitutiveBase const & constitutiveRelation, LAMBDA && lambda )
  {
    ConstitutivePassThruHandler< CompressibleSolid< PressurePorosity, ConstantPermeability >,
                                 CompressibleSolid< PressurePorosity, CarmanKozenyPermeability >,
                                 CompressibleSolid< PressurePorosity, ExponentialDecayPermeability >,
                                 CompressibleSolid< PressurePorosity, ParallelPlatesPermeability >,
                                 CompressibleSolid< PressurePorosity, PressurePermeability >,
                                 CompressibleSolid< PressurePorosity, SlipDependentPermeability >,
                                 CompressibleSolid< PressurePorosity, WillisRichardsPermeability >,
                                 PorousSolid< DruckerPragerExtended, ConstantPermeability >,
                                 PorousSolid< ModifiedCamClay, ConstantPermeability >,
                                 PorousSolid< DelftEgg, ConstantPermeability >,
                                 PorousSolid< DruckerPrager, ConstantPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPrager >, ConstantPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPragerExtended >, ConstantPermeability >,
                                 PorousSolid< DuvautLionsSolid< ModifiedCamClay >, ConstantPermeability >,
                                 PorousSolid< ElasticIsotropic, ConstantPermeability >,
                                 PorousSolid< ElasticTransverseIsotropic, ConstantPermeability >,
                                 PorousSolid< ElasticIsotropicPressureDependent, ConstantPermeability >,
                                 PorousSolid< ElasticOrthotropic, ConstantPermeability >,
                                 PorousSolid< DruckerPragerExtended, CarmanKozenyPermeability >,
                                 PorousSolid< ModifiedCamClay, CarmanKozenyPermeability >,
                                 PorousSolid< DelftEgg, CarmanKozenyPermeability >,
                                 PorousSolid< DruckerPrager, CarmanKozenyPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPrager >, CarmanKozenyPermeability >,
                                 PorousSolid< DuvautLionsSolid< DruckerPragerExtended >, CarmanKozenyPermeability >,
                                 PorousSolid< DuvautLionsSolid< ModifiedCamClay >, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticIsotropic, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticTransverseIsotropic, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticIsotropicPressureDependent, CarmanKozenyPermeability >,
                                 PorousSolid< ElasticOrthotropic, CarmanKozenyPermeability >,
                                 PorousDamageSolid< DamageSpectral< ElasticIsotropic > >,
                                 PorousDamageSolid< DamageVolDev< ElasticIsotropic > >,
                                 PorousDamageSolid< Damage< ElasticIsotropic > >,
                                 ReactiveSolid< ReactivePorosity, ConstantPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, ConstantSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, ConstantPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, PowerLawSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, ConstantPermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, CarmanKozenyPermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PressurePermeability, SubstrateCoverageSurfaceArea >,
                                 ReactiveSolid< ReactivePorosity, PowerLawPermeability, SubstrateCoverageSurfaceArea > >::execute( constitutiveRelation,
                                                                                                                                   std::forward< LAMBDA >( lambda ) );
  }
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_ */
