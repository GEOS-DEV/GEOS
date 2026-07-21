/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file RigidBodyMPMEvent.hpp
 */

#ifndef GEOSX_RIGIDBODYMPM_MPMEVENT_HPP_
#define GEOSX_RIGIDBODYMPM_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class RigidBodyMPMEvent
 *
 * Activates color-partitioned rigid-body MPM between startTime and endTime.
 * One rigid body is formed from all particles sharing a particleColor value.
 */
class RigidBodyMPMEvent : public MPMEventBase
{
public:
  RigidBodyMPMEvent( const string & name,
                     Group * const parent );

  virtual ~RigidBodyMPMEvent() override;

  static string catalogName() { return "RigidBodyMPM"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct
  {
    static constexpr char const * maxGridFieldsString() { return "maxGridFields"; }
    static constexpr char const * linearDampingString() { return "linearDamping"; }
    static constexpr char const * angularDampingString() { return "angularDamping"; }
    static constexpr char const * contactCFLString() { return "contactCFL"; }
    static constexpr char const * maxTimeStepString() { return "maxTimeStep"; }
    static constexpr char const * rigidBodyPenetrationPenaltyBetaString() { return "rigidBodyPenetrationPenaltyBeta"; }
    static constexpr char const * stopKineticEnergyString() { return "stopKineticEnergy"; }
    static constexpr char const * maxForceString() { return "maxForce"; }
    static constexpr char const * minActiveTimeString() { return "minActiveTime"; }
    static constexpr char const * writeHistoryString() { return "writeHistory"; }
    static constexpr char const * historyWriteIntervalString() { return "historyWriteInterval"; }
  } RigidBodyMPMEventViewKeys;

  int getMaxGridFields() const { return m_maxGridFields; }
  real64 getLinearDamping() const { return m_linearDamping; }
  real64 getAngularDamping() const { return m_angularDamping; }
  real64 getContactCFL() const { return m_contactCFL; }
  real64 getMaxTimeStep() const { return m_maxTimeStep; }
  real64 getRigidBodyPenetrationPenaltyBeta() const { return m_rigidBodyPenetrationPenaltyBeta; }
  real64 getStopKineticEnergy() const { return m_stopKineticEnergy; }
  real64 getMaxForce() const { return m_maxForce; }
  real64 getMinActiveTime() const { return m_minActiveTime; }
  int getWriteHistory() const { return m_writeHistory; }
  real64 getHistoryWriteInterval() const { return m_historyWriteInterval; }

protected:
  virtual void postInputInitialization() override final;

  int m_maxGridFields;
  real64 m_linearDamping;
  real64 m_angularDamping;
  real64 m_contactCFL;
  real64 m_maxTimeStep;
  real64 m_rigidBodyPenetrationPenaltyBeta;
  real64 m_stopKineticEnergy;
  real64 m_maxForce;
  real64 m_minActiveTime;
  int m_writeHistory;
  real64 m_historyWriteInterval;
};

} /* namespace geos */

#endif /* GEOSX_RIGIDBODYMPM_MPMEVENT_HPP_ */
