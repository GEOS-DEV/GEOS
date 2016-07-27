
#include "LagrangeHelperFunctions.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"


namespace LagrangeHelperFunctions
{
  // *********************************************************************************************************************
  /**
   * @author R. Settgast
   * @param dt time increment of current time step
   *
   * This function is an aggregate operation that multiplies that performs the following
   * operations for each node:
   *
   * push nodal velocity to half step
   * \f[ v_{a}^{n+1/2} = v_{a}^{n} + a_a^{n} * (dt/2) \f]
   *
   * calculate incremental displacement over the step
   * \f[ uhat_{a}^{n+1/2} = v_{a}^{n+1/2} * dt \f]
   *
   * calculate total displacement at end of step
   * \f[ u_{a}^{n+1} = u_{a}^{n} + uhat_{a}^{n+1/2} \f]
   *
   * zero nodal forces
   * \f[ f_{a}^{n+1} = 0 \f]
   *
   *
   */
   void LinearPointUpdatePart1( ObjectDataStructureBaseT& objectManager,
                                      const realT& time,
                                      const realT& dt,
                                      const bool clearForces)
  {
    if( objectManager.DataLengths() > 0 )
    {
      Array1dT<R1Tensor>& velocity = objectManager.GetFieldData<FieldInfo::velocity> ();
      const Array1dT<R1Tensor>& acceleration = objectManager.GetFieldData<FieldInfo::acceleration> ();
      Array1dT<R1Tensor>& incrementalDisplacement = objectManager.GetFieldData<FieldInfo::incrementalDisplacement> ();
      Array1dT<R1Tensor>& displacement = objectManager.GetFieldData<FieldInfo::displacement> ();
      Array1dT<R1Tensor>& force = objectManager.GetFieldData<FieldInfo::force> ();
      Array1dT<R1Tensor>* const hgForce = objectManager.GetFieldDataPointer<FieldInfo::hgforce> ();
      Array1dT<R1Tensor>* const dampingForce = objectManager.GetFieldDataPointer<R1Tensor> ("dampingForce");

      rArray1d& work = objectManager.GetFieldData<realT>("work");


      const realT dtdiv2 = 0.5 * dt;

      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {
        // push nodal velocity forward to the half-step
        velocity[a].plus_cA(dtdiv2, acceleration[a]);
      }

      BoundaryConditionFunctions::ApplyDirichletBoundaryCondition<R1Tensor>(objectManager,
                                                                            Field<FieldInfo::velocity>::Name(),
                                                                            time+0.5*dt);


      BoundaryConditionFunctions::ApplyRigidWallBoundaryCondition( objectManager, dt );

      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {
        // calculate incremental displacements over the step
        incrementalDisplacement[a].cA(dt, velocity[a]);

        work[a] += 0.5 * Dot(force[a],incrementalDisplacement[a]);

        // add incremental displacements to total displacements to bring to end of step
        displacement[a] += incrementalDisplacement[a];

      }


  /*
      const Array1dT<lArray1d>& childIndices = objectManager.GetVariableOneToManyMap( "childIndices" );
      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {
        if( !(childIndices[a].empty()) )
        {
          incrementalDisplacement[a] = 0.5*( incrementalDisplacement[childIndices[a][0]] + incrementalDisplacement[childIndices[a][1]] );
          displacement[a] = 0.5*( displacement[childIndices[a][0]] + displacement[childIndices[a][1]] );
        }
      }*/

      // set forces to zero
      if(clearForces)
      {
        force = 0.0;
        if( hgForce!=nullptr )
        {
          *hgForce = 0.0;
        }
        if( dampingForce!=nullptr )
        {
          *dampingForce = 0.0;
        }
      }
    }
  }





  // *********************************************************************************************************************
  /**
   * @author R. Settgast
   * @param dt time increment of current time step
   *
   * This function is an aggregate operation that multiplies that performs the following
   * operations for each node:
   *
   * calculate nodal acceleration at end of step
   * \f[ a_{a}^{n+1} = f_{a}^{n+1} / m_a^{n+1} \f]
   *
   * push nodal velocity to end of step
   * \f[ v_{a}^{n+1} = v_{a}^{n+1/2} + a_a^{n+1} * (dt/2) \f]
   *
   */
   void LinearPointUpdatePart2( ObjectDataStructureBaseT& objectManager,
                                const realT& time,
                                const realT& dt,
                                const R1Tensor& gravity,
                                const realT& massDamping,
                                const realT& stiffDamping )
  {
     ///TODO why does this change answers???
    if( objectManager.DataLengths() > 0 ) //&& dt>0.0 )
    {
      Array1dT<R1Tensor>& velocity = objectManager.GetFieldData<FieldInfo::velocity> ();
      Array1dT<R1Tensor>& acceleration = objectManager.GetFieldData<FieldInfo::acceleration> ();
      Array1dT<R1Tensor>& force = objectManager.GetFieldData<FieldInfo::force> ();

      Array1dT<R1Tensor>* const dampingForce = objectManager.GetFieldDataPointer<R1Tensor>("dampingForce");
      Array1dT<realT>& mass = objectManager.GetFieldData<FieldInfo::mass> ();
      rArray1d& work = objectManager.GetFieldData<realT>("work");

      const Array1dT<R1Tensor>& incrementalDisplacement = objectManager.GetFieldData<FieldInfo::incrementalDisplacement> ();


      const realT dtdiv2 = 0.5 * dt;
      realT dampedMass;
      R1Tensor massDampingForce;
      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {
        // calculate acceleration

        force[a] += gravity * mass[a];
#if 0
        if( massDamping > 0.0 )
        {
          dampedMass = mass[a] * ( 1 + 0.5 * dt * massDamping );
          massDampingForce.cA( -massDamping * mass[a], velocity[a] );
          massDampingForce += force[a];
          acceleration[a].Adivc(dampedMass, massDampingForce );
        }
#else
        if( massDamping > 0.0 || stiffDamping > 0.0  )
        {
          dampedMass = mass[a] * ( 1 + 0.5 * dt * massDamping );
          massDampingForce.cA( -massDamping * mass[a], velocity[a] );

          if( stiffDamping > 0.0 && dampingForce != nullptr && dt>0.0 )
          {
            for( int i=0 ; i<3 ; ++i )
            {
              if( fabs((*dampingForce)[a][i]) > 0.0 )
              {
                realT h = ( force[a][i] + mass[a] * velocity[a][i] / ( 0.5 * dt ) ) / (-(*dampingForce)[a][i]);

                /*
                realT h = ( force[a][i] + 0.5 * mass[a] * velocity[a][i] * ( 1.0 / ( 0.5 * dt ) ) - massDamping ) / (-dampingForce[a][i]);
                */
                if( h<0 )
                {
                  h = 0.0;
                }
                else if( h>1.0 )
                {
                  h = 1.0;
                }
                (*dampingForce)[a][i] *= h;
              }
            }
          }

          (*dampingForce)[a] += massDampingForce;
          massDampingForce = (*dampingForce)[a];
          massDampingForce += force[a];
          acceleration[a].Adivc(dampedMass, massDampingForce );
        }
#endif
        else
        {
          acceleration[a].Adivc(mass[a], force[a]);
        }

        //acceleration[a] += gravity * mass[a] / dampedMass;

        work[a] += 0.5 * Dot(force[a],incrementalDisplacement[a]);
      }

      BoundaryConditionFunctions::ApplyDirichletBoundaryCondition<R1Tensor>(objectManager,
                                                                            Field<FieldInfo::acceleration>::Name(),time+dt);

      for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
      {
        // push velocity forward to end of step
        velocity[a].plus_cA(dtdiv2, acceleration[a]);
      }

      BoundaryConditionFunctions::ApplyDirichletBoundaryCondition<R1Tensor>(objectManager,
                                                                            Field<FieldInfo::velocity>::Name(),time+dt);
    }
  }





  /**
   * @brief First part of the rotation update (also updates current positions of de's and nodes)
   * @author Scott Johnson
   *
   * This function is an aggregate operation that multiplies that performs the following
   * operations for each node (after Omelyan, 1998). Note that this interpretation of the
   * quaternion can be parameterized using the rotation, \f[\theta\f], about a unit axis, n:
   * \f[ q[0] = cos(0.5 * \theta) \f]
   * \f[ q[1..3] = sin(0.5 * \theta) * n\f]
   *
   * calculate quaternion time derivative at the step
   * \f[ \dot{q}_{a}^{n} = \frac{1}{2}Q\left(q_{a}^{n}\right)\omega_{a}^{n} \f]
   * \f[ Q\left(q_{a}^{n}\right)=\left[\begin{array}{ccc}-q_{a}^{n}[1] & -q_{a}^{n}[2] & -q_{a}^{n}[3] \\ q_{a}^{n}[0] & -q_{a}^{n}[3] & q_{a}^{n}[2] \\ q_{a}^{n}[3] & q_{a}^{n}[0] & -q_{a}^{n}[1] \\ -q_{a}^{n}[2] & q_{a}^{n}[1] & q_{a}^{n}[0]\end{array}\right]\f]
   *
   * push quaternion to the half step
   * \f[ q_{a}^{n+1/2} = q_{a}^{n} + \dot{q}_{a}^{n} * (dt/2) \f]
   *
   * push body frame nodal rotational velocity to half step
   * \f[ \omega_{a}^{n+1/2} = \omega_{a}^{n} + \alpha_a^{n} * (dt/2) \f]
   *
   * calculate quaternion time derivative at the half step
   * \f[ \dot{q}_{a}^{n+1/2} = \frac{1}{2}Q\left(q_{a}^{n+1/2}\right)\omega_{a}^{n+1/2} \f]
   * \f[ Q\left(q_{a}^{n+1/2}\right)=\left[\begin{array}{ccc}-q_{a}^{n+1/2}[1] & -q_{a}^{n+1/2}[2] & -q_{a}^{n+1/2}[3] \\ q_{a}^{n+1/2}[0] & -q_{a}^{n+1/2}[3] & q_{a}^{n+1/2}[2] \\ q_{a}^{n+1/2}[3] & q_{a}^{n+1/2}[0] & -q_{a}^{n+1/2}[1] \\ -q_{a}^{n+1/2}[2] & q_{a}^{n+1/2}[1] & q_{a}^{n+1/2}[0]\end{array}\right]\f]
   *
   * calculate incremental quaternion over the step
   * \f[ {dq}_{a}^{n+1/2} = \dot{q}_{a}^{n+1/2} * dt \f]
   *
   * calculate quaternion at end of step
   * \f[ q_{a}^{n+1} = q_{a}^{n} + {dq}_{a}^{n+1/2} \f]
   *
   * renormalize quaternion
   * \f[ |q| = 1 \f]
   *
   * @param discreteElementManager manager of the discrete elements to update
   * @param time current simulation time
   * @param dt current timestep length
   */
   void RotationalPointUpdatePart1(DiscreteElementManagerBaseT& discreteElementManager,
                                         const realT& time, const realT& dt)
  {

    if( discreteElementManager.DataLengths() > 0 )
    {

    Array1dT<R1Tensor>& rotationalVelocity        = discreteElementManager.GetFieldData<FieldInfo::rotationalVelocity> ();
    Array1dT<R1Tensor>& rotationalAcceleration    = discreteElementManager.GetFieldData<FieldInfo::rotationalAcceleration> ();
    Array1dT<R1Tensor>& rotationalAxisIncrement   = discreteElementManager.GetFieldData<FieldInfo::rotationalAxisIncrement> ();
    Array1dT<realT>& rotationalMagnitudeIncrement = discreteElementManager.GetFieldData<FieldInfo::rotationalMagnitudeIncrement> ();
    Array1dT<R1Tensor>& rotationAxis              = discreteElementManager.GetFieldData<FieldInfo::rotationAxis> ();
    Array1dT<realT>& rotationMagnitude            = discreteElementManager.GetFieldData<FieldInfo::rotationMagnitude> ();

    const realT dtdiv2 = 0.5 * dt;

    //note: all quantities are in body frame
    for (localIndex a = 0; a < discreteElementManager.DataLengths(); ++a)
    {
      //push the quaternion to the half step
      realT qh[] = {rotationMagnitude[a], rotationAxis[a][0], rotationAxis[a][1], rotationAxis[a][2]};
      realT dqdt[nsdof + 1];
      discreteElementManager.Calculate_dqdt(a, qh, dqdt);
      qh[0] = rotationMagnitude[a] + dtdiv2 * dqdt[0];
      for (unsigned int i = 0; i < nsdof; i++)
        qh[i + 1] = rotationAxis[a][i] + dtdiv2 * dqdt[i + 1];

      // push rotational velocity forward to the half-step
      rotationalVelocity[a].plus_cA(dtdiv2, rotationalAcceleration[a]);

      // calculate dq based on the half step
      discreteElementManager.Calculate_dqdt(a, qh, dqdt);
      rotationalMagnitudeIncrement[a] = dqdt[0] * dt;
      rotationalAxisIncrement[a][0] = dqdt[1] * dt;
      rotationalAxisIncrement[a][1] = dqdt[2] * dt;
      rotationalAxisIncrement[a][2] = dqdt[3] * dt;

      // add incremental rotational displacements to total displacements to bring to end of step
      rotationMagnitude[a] += rotationalMagnitudeIncrement[a];
      rotationAxis[a] += rotationalAxisIncrement[a];

      // renormalize quaternion if necessary
      {
        realT fct = rotationMagnitude[a] * rotationMagnitude[a];
        for (unsigned int i = 0; i < nsdof; i++)
          fct += rotationAxis[a][i] * rotationAxis[a][i];
        if (fct > 0 && !isEqual(fct,1.0) )
        {
          fct = 1. / sqrt(fct);
          rotationMagnitude[a] *= fct;
          for (unsigned int i = 0; i < nsdof; i++)
            rotationAxis[a][i] *= fct;
        }
      }
    }

    }
  }

   void RotationalPointUpdatePart1b(DiscreteElementManagerT& discreteElementManager)
  {
    if( discreteElementManager.DataLengths() > 0 )
    {
      //discrete element positions
      Array1dT<R1Tensor>& deCurrentPosition         = discreteElementManager.GetFieldData<FieldInfo::currentPosition> ();
      Array1dT<R1Tensor>& deReferencePosition       = discreteElementManager.GetFieldData<FieldInfo::referencePosition> ();
      Array1dT<R1Tensor>& deDisplacement            = discreteElementManager.GetFieldData<FieldInfo::displacement> ();

      //nodal positions
      Array1dT<R1Tensor>& nodeRelativePosition      = discreteElementManager.m_nodeManager->GetFieldData<FieldInfo::relativePosition> ();
      Array1dT<R1Tensor>& nodeCurrentPosition       = discreteElementManager.m_nodeManager->GetFieldData<FieldInfo::currentPosition> ();

      //update nodal and discrete element current positions
      //note: all quantities are in body frame
      for (localIndex a = 0; a < discreteElementManager.DataLengths(); ++a)
      {
        R2Tensor rotation;
        discreteElementManager.RotationTensor(a, rotation);

        //--APPLY THE ROTATION TENSOR--
        deCurrentPosition[a] = deReferencePosition[a];
        deCurrentPosition[a] += deDisplacement[a];
        for(localIndex b=0; b<discreteElementManager.m_discreteElementToExternalNodesMap[a].size(); b++)
        {
          int nn = discreteElementManager.m_discreteElementToExternalNodesMap[a][b];
          //****local to global direction transform (see DiscreteElementManagerBaseT.h)
          //    R2Tensor Rt;
          //    RotationTensorTranspose(a, Rt);
          //    global.AijBj(Rt, local);
          nodeCurrentPosition[nn].AijBi(rotation, nodeRelativePosition[nn]);
          nodeCurrentPosition[nn] += deCurrentPosition[a];
        }
      }
    }
  }

  // *********************************************************************************************************************
  /**
   * @author Scott Johnson
   * @param objectManager manager of objects to operate on
   * @param time current simulation time
   * @param dt time increment of current time step
   *
   * This function is an aggregate operation that multiplies that performs the following
   * operations for each node:
   *
   * calculate body frame nodal rotational acceleration at end of step
   * \f[ \alpha_{a}^{n+1} = \tau_{a}^{n+1} / I_a^{n+1} \f]
   *
   * push body frame nodal rotational velocity to end of step
   * \f[ \omega_{a}^{n+1} = \omega_{a}^{n+1/2} + \alpha_a^{n+1} * (dt/2) \f]
   *
   * zero body frame moments
   * \f[ m_{a}^{n+1} = 0 \f]
   *
   */
   void RotationalPointUpdatePart2(ObjectDataStructureBaseT& objectManager,
                                         const realT& time, const realT& dt)
  {
    if( objectManager.DataLengths() > 0 )
    {
    Array1dT<R1Tensor>& rotationalVelocity = objectManager.GetFieldData< FieldInfo::rotationalVelocity> ();
    Array1dT<R1Tensor>& rotationalAcceleration = objectManager.GetFieldData< FieldInfo::rotationalAcceleration> ();
    Array1dT<R1Tensor>& moment = objectManager.GetFieldData<FieldInfo::moment> ();
    Array1dT<R1Tensor>& rotationalInertia = objectManager.GetFieldData<FieldInfo::rotationalInertia> ();

    const realT dtdiv2 = 0.5 * dt;
    for (localIndex a = 0; a < objectManager.DataLengths(); ++a)
    {
      // calculate rotational acceleration
      rotationalAcceleration[a] = moment[a];
      rotationalAcceleration[a] /= rotationalInertia[a];

      // push rotational velocity forward to end of step
      rotationalVelocity[a].plus_cA(dtdiv2, rotationalAcceleration[a]);
    }

    // set moments to zero
    moment = 0.0;
    }
  }


   realT CalculateMaxStableExplicitTimestep( const realT& density,
                                             const realT& stiffness,
                                             const realT& BB )
   {
     realT rval =  sqrt( density / ( stiffness *BB ) );

     return rval;
   }

}
