/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file RigidBodyMPMEvent.cpp
 */

#include "RigidBodyMPMEvent.hpp"

namespace geos
{

using namespace dataRepository;

RigidBodyMPMEvent::RigidBodyMPMEvent( const string & name,
                                      Group * const parent ):
  MPMEventBase( name, parent ),
  m_maxGridFields( 4 ),
  m_linearDamping( 0.0 ),
  m_angularDamping( 0.0 ),
  m_contactCFL( -1.0 ),
  m_maxTimeStep( -1.0 ),
  m_rigidBodyPenetrationPenaltyBeta( 0.0 ),
  m_stopKineticEnergy( -1.0 ),
  m_maxForce( -1.0 ),
  m_minActiveTime( 0.0 ),
  m_writeHistory( 1 ),
  m_historyWriteInterval( 0.0 )
{
  registerWrapper( viewKeyStruct::maxGridFieldsString(), &m_maxGridFields ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_maxGridFields ).
    setDescription( "Maximum number of local nodal velocity fields used for rigid-body color partitioning. "
                    "The first maxGridFields-1 colors are mapped to pure fields; additional colors are "
                    "lumped into the last field with weighted force attribution." );

  registerWrapper( viewKeyStruct::linearDampingString(), &m_linearDamping ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_linearDamping ).
    setDescription( "Rigid-body linear viscous damping rate applied after force integration." );

  registerWrapper( viewKeyStruct::angularDampingString(), &m_angularDamping ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_angularDamping ).
    setDescription( "Rigid-body angular viscous damping rate applied after torque integration." );

  registerWrapper( viewKeyStruct::contactCFLString(), &m_contactCFL ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_contactCFL ).
    setDescription( "Optional rigid-event-only contact CFL. When positive, the rigid event time step is limited by "
                    "contactCFL * h / max(2*maxParticleSpeed, maxParticleSpeed + maxBoundarySpeed). "
                    "Nonpositive values disable this additional limiter." );

  registerWrapper( viewKeyStruct::maxTimeStepString(), &m_maxTimeStep ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_maxTimeStep ).
    setDescription( "Optional rigid-event-only maximum time step. Nonpositive values disable this cap. "
                    "The continuum MPM time-step calculation is unchanged." );

  registerWrapper( viewKeyStruct::rigidBodyPenetrationPenaltyBetaString(), &m_rigidBodyPenetrationPenaltyBeta ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_rigidBodyPenetrationPenaltyBeta ).
    setDescription( "Dimensionless rigid-event-only penetration penalty coefficient. For material contact, the "
                    "normal penalty force is beta * m_eff * penetration / dt^2. For a moving contact boundary, "
                    "the equivalent nodal penalty uses the mapped nodal mass. Zero disables the penalty." );

  registerWrapper( viewKeyStruct::stopKineticEnergyString(), &m_stopKineticEnergy ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_stopKineticEnergy ).
    setDescription( "Optional event-completion kinetic-energy threshold. Negative disables this stop criterion." );

  registerWrapper( viewKeyStruct::maxForceString(), &m_maxForce ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_maxForce ).
    setDescription( "Optional event-completion maximum resultant body force. Negative disables this stop criterion." );

  registerWrapper( viewKeyStruct::minActiveTimeString(), &m_minActiveTime ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_minActiveTime ).
    setDescription( "Minimum active time before optional RigidBodyMPM force or kinetic-energy stop criteria are allowed to complete the event. "
                    "Use this to prevent a transient first wall contact from switching the run to continuum before the compaction ramp has completed." );

  registerWrapper( viewKeyStruct::writeHistoryString(), &m_writeHistory ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_writeHistory ).
    setDescription( "Flag enabling rigid-body CSV history output." );

  registerWrapper( viewKeyStruct::historyWriteIntervalString(), &m_historyWriteInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_historyWriteInterval ).
    setDescription( "Time interval for rigidBodyHistory.csv. Zero writes every solver step." );
}

RigidBodyMPMEvent::~RigidBodyMPMEvent()
{}

void RigidBodyMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();
  GEOS_ERROR_IF( m_maxGridFields < 2, "RigidBodyMPM maxGridFields must be at least 2." );
  GEOS_ERROR_IF( m_linearDamping < 0.0, "RigidBodyMPM linearDamping must be non-negative." );
  GEOS_ERROR_IF( m_angularDamping < 0.0, "RigidBodyMPM angularDamping must be non-negative." );
  GEOS_ERROR_IF( m_rigidBodyPenetrationPenaltyBeta < 0.0,
                 "RigidBodyMPM rigidBodyPenetrationPenaltyBeta must be non-negative." );
  GEOS_ERROR_IF( m_minActiveTime < 0.0, "RigidBodyMPM minActiveTime must be non-negative." );
}

REGISTER_CATALOG_ENTRY( MPMEventBase, RigidBodyMPMEvent, string const &, Group * const )

} /* namespace geos */
