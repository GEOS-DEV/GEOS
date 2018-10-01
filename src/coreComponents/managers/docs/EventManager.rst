###############################################################################
Event Management
###############################################################################

The goal of the GEOSX event manager is to be flexible with regards to event type, application order, and method of triggering.  The event manager is configured via the ``Event`` block in an input .xml file, i.e.:

.. code-block:: xml
  <Events maxTime="1.0e-2">
    <PeriodicEvent name="event_a"
                   target="/path/to/event" />
  </Events>


Event Manager Configuration
=====================================

The Event block includes two attributes (by default, they are set to their max values):
* ``maxTime`` - Sets the maximum time for the global event loop (real64, optional)
* ``maxCycle`` - Sets the maximum number of cycles for the global event loop (integer, optional)

Event candidates are indicated by appending children to the Event block.  These children must point to an object of the type ``EventBase``.  The common attributes for all events are:
* ``name`` - A unique identifier for the event (string)
* ``target`` - A unix-style path to the object that should be executed if the event criteria are met.  The path can be either absolute (i.e.: '/Solvers/solver_a') or relative (i.e.: '../../Solvers/solver_a') (string, optional)
* ``beginTime`` - This requires that (time >= beginTime) to execute. (real64, optional)
* ``endTime`` - This requires that (time < endTime) to execute. (real64, optional)
* ``forceDt`` - This will override the timestep requests from its target (real64, optional)
* ``allowSuperstep`` - This will override the time-stepping behavior for its targets, and is explained further below (integer, optional)
* ``allowSubstep`` - This will override the time-stepping behavior for its targets, and is explained further below (integer, optional)
* ``substepFactor`` - This sets the substepping behavior for the target (integer, optional)

The primary type of event used in GEOSX is of type ``PeriodicEvent``.  As its name suggests, it will execute periodically during a simulation.  It can be triggered based upon a cycleFrequency, timeFrequency, a time-function, or a function applied to an object.  The unique attributes for this event are:
* ``cycleFrequency`` - This will instruct the event to execute every N cycles.  A value of "1" (default) will cause the event to trigger every cycle, a value of "2" will trigger every other cycle, and so on. (integer, optional)
* ``timeFrequency`` - This will instruct the event to execute every X seconds.  If this parameter is set, it will supersede the cycle-driven behavior. (real64, optional)
* ``targetExactTimestep`` - If this is set, will allow the event to limit its timestep requests in an attempt to execute on integer multiples of timeFrequency. (bool, optional)
* ``function`` - If this is set, the event will evaluate a function to test if its target should execute.  Because some functions may be time-consuming to compute, the function is only evaluated after the cycle/time criteria are met.  The function can be a function of time or can be applied to an object. (string, optional) 
* ``threshold`` - If the optional function control is used, the event will execute if f(inputs) > threshold.  The default value is 0.  (real64, optional)
* ``object`` - If this value is set, the function will be applied to an object, and the min, mean, or max value of the function will be compared to the threshold. (string, optional)
* ``set`` - If the target of a function is an object, then this may be used to limit the sets within the object to apply the function to.  Otherwise, it will be applied to the entire object. (string, optional)
* ``stat`` - If the target of a function is an object, then this will select which value to compare to the threshold. 0 = min, 1 = mean, 2 = max.  (integer, optional)

The other event type used in GEOSX is of the type ``HaltEvent``.  This event will track the wall clock, and if it is executed it will set a flag that instructs the manager to exit.  The unique attribute for this object is:
* ``maxRunTime`` - The event will trigger once (wallTime > maxRunTime) (real64)


Basic Event Execution Rules
=====================================

During a simulation, the event manager will loop through the list of the events **in the order they are defined in the xml**.  The simulation *cycle* denotes the number of times this loop has completed, and *dt* denotes the timestep.  During each loop, each event will do the following:

1. Calculate a *forecast*, which is defined as the expected number of cycles until the event is expected to execute.
2. ``if (forecast == 1)`` the event will signal its target to prepare to execute.  This is useful for preparing time-consuming I/O operations.
3. ``if (forecast <= 0)`` the event will execute its target
4. ``if (forecast <= 1)`` the event will obtain a timestep request from its target for the next cycle
5. Check to see if the main loop execution flag has been set

To initialize the simulation, the value of *dt* for the first cycle is set to 0.  At the end of each loop, the *dt* for the next cycle will be set to the smallest timestep requested by the events.  The event manager loop will continue until the maximum time, maximum number of cycles, and/or the exit flag is set.  After exiting the main loop, the event manager will call the ``Cleanup`` method for each of its children (to produce final plots, etc.).


Event Progress Indicator
=====================================
Because the event manager allows the user to specify the order of events, it could introduce ambiguity into the timestamps of output files.  To resolve this, we pass the *progress*, which is defined as the percent completion of the main loop, to the event targets.  Currently, this value is included in the headers of plot files.

The event manager will also test to see if a given target is expected to execute **after all** calls to objects of type ``SolverBase``.  If this is the case, then the event will be executed with ``time = time + dt``.  Otherwise, the event will be executed with ``time = time``.  This is useful for automatically aligning the timestamps for output files.



Event Superstepping and Substepping Behavior
=====================================

If the ``allowSuperstep`` attribute of an event is set, when its criteria are met, it will execute its target with ``time = lastTime`` and ``dt = dt + time - lastTime`` instead of their typical values.

If the ``allowSubstep`` attribute of an event is set, when its criteria are met, it will execute its target ``N = substepFactor`` times with ``dt = dt_o / N`` and an the appropriate timestamp.


Event Forecast Calculation
=====================================
Again, the *forecast* is defined as the expected number of cycles until the event will execute.  If ``time < beginTime`` or ``time >= endTime``, this value will be equal to ``std::numeric_limits<integer>::max()``.  Otherwise, it is calculated by the specific event types:
* cycle-driven ``PeriodicEvent`` - ``forecast = cycleFrequency - (cycle - lastCycle)``
* time-driven ``PeriodicEvent`` - if (dt > 0), ``forecast = (timeFrequency - (time - lastTime)) / dt``, otherwise forecast is set to the max value.
For 
* ``HaltEvent`` - ``forecast = (maxRuntime - (currentTime - startTime)) / realDt``



Nested Events
=====================================
The event manager allows its child events to be nested.  If this feature is used, then the manager follows the basic execution rules, with the following exception:  When its criteria are met, an event will first execute its (optional) target.  It will then estimate the forecast for its own sub-events, and execute them following the same rules as in the main loop.  For example:

.. code-block:: xml
  <Events maxTime="1.0e-2">
    <PeriodicEvent name="event_a"
                   target="/path/to/target_a" />

    <PeriodicEvent name="event_b"
                   timeFrequency="100"

      <PeriodicEvent name="subevent_b_a"
                     target="/path/to/target_c"/>

      <PeriodicEvent name="subevent_b_b"
                     target="/path/to/target_d"/>
    <PeriodicEvent/>
  </Events>

In this example, event_a will trigger during every cycle and call the Execute method on the object located at /path/to/target_a.  Because it is time-driven, event_b will execute every 100 s.  When this occurs, it will execute subevent_b_a and subevent_b_b in order. (Note: these are both cycle-driven events which, by default would occur every cycle.  However, they will not execute until each of their parents, grandparents, etc. execution criteria are met as well.)

