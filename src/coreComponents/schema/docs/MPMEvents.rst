

======== ======= ======= =========== 
Name     Type    Default Description 
======== ======= ======= =========== 
logLevel integer 0       Log level   
======== ======= ======= =========== 



Uniform confining pressure
-------------------------

``ConfiningPressure`` accepts ``uniformPressure=1`` (default: 0), without
``confiningPressureBoxMin``, ``confiningPressureBoxMax`` or
``confiningPressureSphereRadius``. It applies the interpolated pressure as a
uniform background stress. ``startPressure`` and ``endPressure`` remain
required; ``startTime``, ``endTime`` and ``interpType`` retain their usual event
semantics. Axisymmetric MPM requires this uniform mode. Box/sphere selection
remains available in Cartesian formulations. Pressure load state is retained
after the event and in restart files.

See the PFW ``examples/axisymmetricJacketedSample`` README for a runnable case.
