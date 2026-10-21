.. Generated from current source wrappers by extract_mpm_metadata.py.
.. Full inherited XML schema generation requires a GEOS build.

Current directly registered input attributes are listed below. Defaults
are resolved from source where possible; remaining C++ expressions are
shown verbatim. Inherited attributes are not repeated.

.. list-table:: SolidMechanics_MPM source wrapper reference
   :header-rows: 1
   :widths: 24 16 18 42

   * - Attribute
     - Input / type
     - Default
     - Description
   * - FSubcycles
     - OPTIONAL / integer
     - 1
     - Deprecated compatibility input. The fixed-endpoint matrix exponential uses internal scaling and squaring; external F subcycles are not used.
   * - LBar
     - OPTIONAL / int
     - 0
     - Used to choose the Lbar method. 0 = do nothing, 1 = run trD through the SPH kernel, 2 = straight-up volume average
   * - LBarScale
     - OPTIONAL / real64
     - 0.0
     - Scales the contribution of volume-averaged trD to particle L, alleviates checkerboarding
   * - LRtolerance
     - OPTIONAL / real64
     - 0.01
     - Tolerance for logistic regression calculations
   * - applyVelocityGradientWhenSplittingParticles
     - OPTIONAL / int
     - 1
     - Flag for assigning split-particle velocities with the parent velocity gradient: v_child = v_parent + L_parent * dx_child. Mass-weighted centroid offsets conserve linear momentum; kinetic energy is recomputed and can increase. Zero copies parent velocity.
   * - areaIntegrationMethod
     - OPTIONAL / string (enum)
     - BruteForce
     - Method for performing nodal area integration
   * - batchSize
     - OPTIONAL / localIndex
     - std::numeric_limits< localIndex >::max()
     - Batch sizing for constitutive model update kernels on device
   * - bcTable
     - OPTIONAL / array2d< real64 >
     - not explicitly set
     - Array that stores time-dependent bc types on x-, x+, y-, y+, z- and z+ faces.
   * - binSizeMultiplier
     - OPTIONAL / int
     - 1
     - Multiplier for setting bin size, used to speed up particle neighbor sorting
   * - bodyForce
     - OPTIONAL / array1d< real64 >
     - not explicitly set
     - Array that stores uniform body force
   * - boundaryConditionTypes
     - OPTIONAL / array1d< int >
     - not explicitly set
     - Boundary conditions on x-, x+, y-, y+, z- and z+ faces. Options are:; * ; *
   * - boundaryFaceCoefficientsOfRestitution
     - OPTIONAL / array1d< real64 >
     - not explicitly set
     - Boundary face coefficients of restitution for BC type = 3 on x-, x+, y-, y+, z- and z+ faces.
   * - boundaryFaceFrictionCoefficients
     - OPTIONAL / array1d< real64 >
     - not explicitly set
     - Boundary face fricition coefficients for BC type = 3 on x-, x+, y-, y+, z- and z+ faces.
   * - boxAverageHistory
     - OPTIONAL / int
     - 0
     - Flag for whether to output box average history
   * - boxAverageMax
     - OPTIONAL / array1d< real64 >
     - not explicitly set
     - Maximum corner position of box average
   * - boxAverageMin
     - OPTIONAL / array1d< real64 >
     - not explicitly set
     - Minimum corner position of box average
   * - boxAverageResizeWithDomain
     - OPTIONAL / int
     - 0
     - Flag for whether the box-average integration domain follows the current diagonal domain stretch.
   * - boxAverageWriteInterval
     - OPTIONAL / real64
     - 0.0
     - Interval between writing box averages to files
   * - computeCZInterfacesFromDamage
     - OPTIONAL / int
     - 0
     - Flag to compute CZ interfaces from damage
   * - computeInternalEnergyAndTemperature
     - OPTIONAL / int
     - 0
     - Flag to enable computing changes to internal energy and temperature
   * - computeNodalArea
     - OPTIONAL / int
     - 0
     - Option to compute nodal area from grid surface normals and positions (results may be inaccurate if implicit surface normals and positiosn are used)
   * - computeParticleSurfaceNormalsAndPositions
     - OPTIONAL / int
     - 0
     - Flag to compute particle surface normals and positions
   * - computeSPHJacobian
     - OPTIONAL / int
     - 0
     - Overlap correction flag (0=off, 1=increase normal force to remove gap. (not fully developed), 2=compute the SPH Jacobian and use it to scale particle density to improve the overlap. )
   * - contactGapCorrection
     - OPTIONAL / string (enum)
     - Implicit
     - Flag for mitigating contact gaps
   * - contactNormalExponent
     - OPTIONAL / real64
     - 1.0
     - Exponent value for surface normal weights in contact normal type aligned
   * - contactNormalType
     - OPTIONAL / string (enum)
     - MassWeighted
     - Flag for contact normal type
   * - cpdiDomainScaling
     - OPTIONAL / int
     - 0
     - Option for CPDI domain scaling
   * - cpdiDomainScalingType
     - OPTIONAL / string (enum)
     - homel
     - CPDI domain-scaling algorithm. Options are:; * ; *
   * - crackTipDetectionThreshold
     - OPTIONAL / real64
     - 0.5
     - Set threshold for crack-tip detection
   * - damageFieldPartitioning
     - OPTIONAL / int
     - 0
     - Flag for using the gradient of the particle damage field to partition material into separate velocity fields
   * - damageHessianSurfaceThreshold
     - OPTIONAL / real64
     - 1e16
     - Threshold to update surface flag of damaged material
   * - defGradResetMaxParticleDomainToGridCellRatio
     - OPTIONAL / real64
     - -1.0
     - Maximum full particle-domain diagonal divided by active grid-cell diagonal before optional oversized-particle treatment. Deformation-gradient reset fallback is intended for the same fluid, melted, damaged, or hypo-elastic material classes as resetDefGradForOversizedParticles.
   * - directionalOverlapCorrection
     - OPTIONAL / int
     - 0
     - Flag for mitigating pile-up of particles at contact interfaces
   * - disableSurfaceNormalsAndPositionsOnCPDIScaling
     - OPTIONAL / int
     - 0
     - Option for disabling explicit surface normals and positions when a particle has severely deformed
   * - disableSurfaceNormalsAndPositionsOnDamage
     - OPTIONAL / int
     - 1
     - Option for disabling explicit surface normals and positions when a particle has been severely damaged
   * - disableSurfaceNormalsAndPositionsOnMelt
     - OPTIONAL / int
     - 1
     - Option for disabling explicit surface normals and positions when a particle has melted
   * - disableSurfaceNormalsAndPositionsOnOversizedDomainReset
     - OPTIONAL / int
     - 0
     - Option for disabling explicit surface normals and positions when an oversized particle domain triggers deformation-gradient reset
   * - disableSurfaceTractionsOnDamage
     - OPTIONAL / int
     - 1
     - Option for disabling explicit surface tractions when a particle has been severely damaged
   * - disableSurfaceTractionsOnMelt
     - OPTIONAL / int
     - 1
     - Option for disabling explicit surface tractions when a particle has melted
   * - domainResetType
     - OPTIONAL / string (enum)
     - isotropicPolar
     - Deformation-gradient reset algorithm for fluid, melted, damaged, explicitly requested scaled-surface particles, and oversized particles treated with deformation-gradient reset. Options are:; * ; *
   * - enableContact
     - OPTIONAL / int
     - 1
     - Flag to enable contact between velocity fields
   * - enablePrescribedBoundaryTransverseVelocities
     - OPTIONAL / array1d< int >
     - not explicitly set
     - Per-face flags enabling prescribed transverse velocity components. If a face is not flagged, moving-boundary enforcement leaves tangential components free.
   * - enableSurfaceTension
     - OPTIONAL / int
     - 0
     - Flag to enable surface tension forces on particles
   * - eventReporting
     - OPTIONAL / int
     - 0
     - Flag to print event status information
   * - exactJIntegration
     - OPTIONAL / int
     - 0
     - Deprecated compatibility input. The endpoint exponential update always integrates particle volume exactly from exp(dt*tr(L)).
   * - explicitSurfaceNormalInfluence
     - OPTIONAL / real64
     - 0.0
     - Dimensionless relative weight for explicit particle surface normals; the solver divides this value by the active grid-cell diagonal before scattering to grid normals and by neighborRadius when forming the DFG partition-reference vector
   * - fTable
     - OPTIONAL / array2d< real64 >
     - not explicitly set
     - Array of time-dependent grid-aligned diagonal stretches interpreted as the prescribed domain F.
   * - fTableInterpType
     - OPTIONAL / string (enum)
     - Linear
     - The type of F table interpolation. Options are 0 (linear), 1 (cosine), 2 (quintic polynomial).
   * - flagParticlesWithBadMappingArraysForDeletion
     - OPTIONAL / int
     - 1
     - Flag particles for deletion if particle-grid mapping construction detects an out-of-bounds particle/stencil location; compact active-ordinal mapping rows before using them.
   * - frictionCoefficient
     - OPTIONAL / real64
     - -1.0
     - Coefficient of friction, currently assumed to be the same everywhere
   * - frictionCoefficientTable
     - OPTIONAL / array2d< real64 >
     - not explicitly set
     - Friction coefficient table for different groups
   * - generalizedVortexMMS
     - OPTIONAL / int
     - 0
     - Option for CPDI domain scaling
   * - geometryReset
     - OPTIONAL
     - FOnly
     - FOnly preserves current material geometry; FAndGeometry resets it to F_new G0. Applies to automatic and healing resets.
   * - gpuScheme
     - OPTIONAL / string (enum)
     - Atomics
     - GPU Scheme (Atomics, No atomics, Random shuffling, Minimal atomics)
   * - gridToParticleMapping
     - OPTIONAL / string (enum)
     - Precomputed
     - Mapping evaluation used by FLIP, PIC, XPIC, and FMPM grid-to-particle updates. Options are:; * ; * . Precomputed reuses the raw arrays built once per step; OnTheFly recomputes the raw node connectivity, shape values, gradients, and field assignment inside each G2P particle kernel. RigidBodyMPM events require Precomputed because their color-field assignments are stored in the field map.
   * - logMomentum
     - OPTIONAL / int
     - not explicitly set
     - Flag for writing end-of-step momentum balance diagnostics to mpm_momentumHistory.csv.
   * - logStartCycle
     - OPTIONAL / int
     - not explicitly set
     - Cycle to start writing mpm_momentumHistory.csv rows.
   * - maxLRIterations
     - OPTIONAL / integer
     - 30
     - Maximum number of iterations for logistic regression calculations
   * - maxParticleJacobian
     - OPTIONAL / real64
     - 10.0
     - Jacobian value above which particles are deleted
   * - maxParticleVelocity
     - OPTIONAL / real64
     - 1e6
     - Velocity above which particles are deleted
   * - maxSingleFieldStateFractionForSeparability
     - OPTIONAL / real64
     - -1.0
     - Maximum volume-weighted mapped fraction of fully damaged, melted, or domain-reset particles allowed for separability. A negative value disables this restriction; active values are in [0,1].
   * - minParticleJacobian
     - OPTIONAL / real64
     - 0.1
     - Jacobian value below which particles are deleted
   * - needsNeighborList
     - OPTIONAL / int
     - 0
     - Flag for whether to construct neighbor list
   * - needsNodalNeighborList
     - OPTIONAL / int
     - 0
     - Flag for whether to construct nodal neighbor list
   * - neighborRadius
     - OPTIONAL / real64
     - -1.0
     - Neighbor radius for SPH-type calculations
   * - normalAndPositionMethod
     - OPTIONAL / string (enum)
     - LogisticRegression
     - Method for computing particle surface normals and positions
   * - numSurfaceIntegrationPoints
     - OPTIONAL / integer
     - 200
     - Number of surface integration points for cohesive grid nodes
   * - overlapCorrection
     - OPTIONAL / string (enum)
     - Off
     - Overlap correction flag (0=off, 1=increase normal force to remove gap. (not fully developed), 2=compute the SPH Jacobian and use it to scale particle density to improve the overlap. )
   * - overlapThreshold1
     - OPTIONAL / real64
     - 1.00
     - Overlap correction threshold 1
   * - overlapThreshold2
     - OPTIONAL / real64
     - 1.10
     - Overlap correction threshold 2
   * - oversizedParticleResetRankParticleCountThreshold
     - OPTIONAL / integer
     - -1
     - Rank-local active master-particle count threshold for oversizedParticleTreatment="split". Oversized particles are split while the start-of-step rank count is less than or equal to this value, and deformation-gradient reset is used only on ranks whose count exceeds it. A negative value disables the reset fallback.
   * - oversizedParticleTreatment
     - OPTIONAL / string (enum)
     - none
     - Action for particles whose current domain diagonal exceeds defGradResetMaxParticleDomainToGridCellRatio. "none" preserves the legacy resetDefGradForOversizedParticles flag, "resetDeformationGradient" resets F, and "split" uses the particle subdivision path unless the rank particle-count threshold requests reset fallback. Allowed values are "none", "resetDeformationGradient", and "split".
   * - overwriteExistingNormalsAndPositions
     - OPTIONAL / int
     - 0
     - Flag to that determines if existing particle normals and positions are overriden
   * - particleDataWriteInterval
     - OPTIONAL / real64
     - 0.0
     - Interval to write particle data to file
   * - planeStrain
     - OPTIONAL / int
     - 0
     - Flag for performing plane strain calculations
   * - plotGridFields
     - OPTIONAL / int
     - 0
     - Flag to enable plotting of grid fields in vtk output
   * - plottableFields
     - OPTIONAL / string_array
     - not explicitly set
     - If specified, determines which fields are written to vtk. Helpful for reducing the vtk output file sizes
   * - prescribedBcTable
     - OPTIONAL / int
     - 0
     - Flag for whether to have time-dependent boundary condition types
   * - prescribedBoundaryFTable
     - OPTIONAL / int
     - not explicitly set
     - Enable moving-boundary deformation from fTable. Boundary velocities use the end-of-step Fdot*Finv rate, while grid geometry uses the exact finite ratio F(t_{n+1})/F(t_n).
   * - prescribedBoundaryTransverseVelocities
     - OPTIONAL / array2d< real64 >
     - not explicitly set
     - Prescribed transverse velocity components for faces whose enablePrescribedBoundaryTransverseVelocities entry is set.
   * - prescribedFTable
     - OPTIONAL / int
     - not explicitly set
     - Flag for whether to have time-dependent superimposed velocity gradient for triply periodic simulations
   * - preventCZInterpenetration
     - OPTIONAL / int
     - 0
     - Flag to use contact normal forces instead of CZ forces in compression
   * - profileCycleInterval
     - OPTIONAL / int
     - 0
     - Cycle interval for MPM profile output. Use either profileCycleInterval or profileWriteInterval. A value of 0 means no cycle-based interval.
   * - profileDirection
     - OPTIONAL / string
     - "x"
     - Direction for MPM profile output. Valid values are x, y, and z.
   * - profileHistory
     - OPTIONAL / int
     - 0
     - Flag to enable MPM profile CSV output. One CSV file is written per requested profile variable.
   * - profileNumSlices
     - OPTIONAL / int
     - 0
     - Number of profile slices. A value of 0 uses the background-grid resolution in profileDirection.
   * - profileVariables
     - OPTIONAL / string_array
     - not explicitly set
     - Variables to write to MPM profile CSV files. Supported names include density, damage, temperature, internalEnergy, supplementalPressure, kineticEnergy, velocityX, velocityY, velocityZ, volumeFraction, plasticStrainMagnitude, stressXX, stressYY, stressZZ, stressYZ, stressXZ, stressXY, and plasticStrainXX/YY/ZZ/YZ/XZ/XY.
   * - profileWriteInterval
     - OPTIONAL / real64
     - 0.0
     - Time interval for MPM profile output. Use either profileWriteInterval or profileCycleInterval. A value of 0 means no time-based interval.
   * - reactionHistory
     - OPTIONAL / int
     - 0
     - Flag for whether to output face reaction history
   * - reactionWriteInterval
     - OPTIONAL / real64
     - 0.0
     - Interval between writing reactions to files
   * - resetDefGradForFullyDamagedParticles
     - OPTIONAL / int
     - 0
     - Flag for resetting deformation gradient of fully damaged particles
   * - resetDefGradForMeltedParticles
     - OPTIONAL / int
     - 0
     - Flag for resetting deformation gradient of melted particles to rotation times isotropic stretch
   * - resetDefGradForOversizedParticles
     - OPTIONAL / int
     - 0
     - Compatibility flag for resetting deformation gradient of particles whose current domain diagonal exceeds the configured grid-cell ratio. Equivalent to oversizedParticleTreatment="resetDeformationGradient" when oversizedParticleTreatment is left at "none". Use only for materials where the stored deviatoric deformation gradient is not constitutive state, such as fluids, melted material, fully damaged material, or hypo-elastic deviatoric models; do not enable for finite-strain hyperelastic, crystal-plasticity, fiber, or other orientation/stretch-history-dependent solids unless the material model explicitly supports this reset.
   * - resetDefGradForScaledSurfaceParticles
     - OPTIONAL / int
     - 0
     - Flag for resetting deformation gradient of CPDI scaled surface particles
   * - rigidBodyMaxGridFields
     - OPTIONAL / int
     - 0
     - Optional solver-level reservation for color-partitioned rigid-body MPM grid fields. RigidBodyMPM events also reserve their maxGridFields during initialization.
   * - separabilityMinDamage
     - OPTIONAL / real64
     - 0.5
     - Damage threshold for field separability
   * - setDomainTemperature
     - OPTIONAL / int
     - not explicitly set
     - Flag that activates domain temperature interpolation from table.
   * - setDomainTemperatureRate
     - OPTIONAL / int
     - not explicitly set
     - Flag that activates domain temperature rateinterpolation from table.
   * - shapeFunctionDiagnostics
     - OPTIONAL / int
     - not explicitly set
     - Flag to perform diagnostics on particle shape functions.
   * - smallMass
     - OPTIONAL / real64
     - DBL_MAX
     - Optional small mass threshold for ignoring extremely low-mass nodes. If omitted, an automatic threshold is computed from the global minimum particle mass.
   * - solverProfiling
     - OPTIONAL / int
     - not explicitly set
     - Flag for timing subroutines in the solver
   * - stressControl
     - OPTIONAL / array1d< int >
     - not explicitly set
     - Flag for whether stress control using box averages is enabled
   * - stressControlAdaptGain
     - OPTIONAL / real64
     - 5.0e-2
     - Normalized rank-one update gain for the adaptive stress-control tangent.
   * - stressControlAdaptStrainWindow
     - OPTIONAL / real64
     - 5.0e-4
     - Accumulated macroscopic strain window for each adaptive tangent update.
   * - stressControlAdaptiveSeek
     - OPTIONAL / int
     - 0
     - Use the adaptive-seek stress controller for normal stress-control components. This opt-in controller learns a normal stress tangent, filters box stress over a strain scale, and enters a bounded compression seek mode when a compressive target has no load path.
   * - stressControlCommandFilterStrain
     - OPTIONAL / real64
     - 0.0
     - Optional macroscopic strain interval used to smooth the applied adaptive stress-control strain rate. A value of zero disables command filtering.
   * - stressControlFilterStrain
     - OPTIONAL / real64
     - 1.0e-3
     - Macroscopic strain interval used to low-pass filter stress for adaptive stress control.
   * - stressControlJammingPressureRatio
     - OPTIONAL / real64
     - 2.0e-2
     - Fraction of target compressive pressure indicating that a load path has appeared during seek mode.
   * - stressControlKd
     - OPTIONAL / real64
     - not explicitly set
     - Derivative gain of stress PID controller
   * - stressControlKi
     - OPTIONAL / real64
     - not explicitly set
     - Integral gain of stress PID controller
   * - stressControlKp
     - OPTIONAL / real64
     - not explicitly set
     - Proportional gain of stress PID controller
   * - stressControlMaxCouplingRatio
     - OPTIONAL / real64
     - 10.0
     - Maximum absolute learned off-diagonal normal tangent divided by the corresponding diagonal tangent.
   * - stressControlMaxFeedbackRateRatio
     - OPTIONAL / real64
     - 2.0
     - Maximum feedback strain-rate norm divided by the prescribed/reference strain-rate norm.
   * - stressControlMaxRateChangeRatio
     - OPTIONAL / real64
     - 0.0
     - Optional per-step cap on the change of an adaptive stress-control strain rate, divided by the reference strain rate. A value of zero disables this cap.
   * - stressControlMaxRateRatio
     - OPTIONAL / real64
     - 5.0
     - Maximum active stress-control strain-rate norm divided by the prescribed/reference strain-rate norm.
   * - stressControlMaxSeekRateRatio
     - OPTIONAL / real64
     - 2.0
     - Maximum seek strain-rate norm divided by the prescribed/reference strain-rate norm.
   * - stressControlMaxSeekStrain
     - OPTIONAL / real64
     - 5.0e-2
     - Maximum accumulated seek strain before a compressive target is declared unreachable.
   * - stressControlMaxSeekStrainIncrement
     - OPTIONAL / real64
     - 1.0e-4
     - Maximum absolute stress-control strain increment in one step during adaptive/seek control.
   * - stressControlMinStrainRate
     - OPTIONAL / real64
     - 1.0e-30
     - Minimum reference strain rate for adaptive stress control when no prescribed deformation rate is active.
   * - stressControlReengageRampStrain
     - OPTIONAL / real64
     - 5.0e-3
     - Macroscopic strain interval over which seek mode is blended back to adaptive active control.
   * - stressControlReengageTangentRatio
     - OPTIONAL / real64
     - 5.0e-3
     - Fraction of initial tangent scale needed to re-engage active control after seek mode.
   * - stressControlResponseStrain
     - OPTIONAL / real64
     - 1.0e-2
     - Macroscopic strain interval over which adaptive stress control removes filtered stress error.
   * - stressControlSeekRateRatio
     - OPTIONAL / real64
     - 0.5
     - Bounded seek strain-rate ratio used while a compressive target has no load path.
   * - stressControlSolverDampingRatio
     - OPTIONAL / real64
     - 5.0e-2
     - Damping used in the small damped least-squares solve for active adaptive stress control.
   * - stressControlStressFloor
     - OPTIONAL / real64
     - 0.0
     - Stress floor used in adaptive stress-control ratios. If zero, a target-based floor is used.
   * - stressControlTangentCeilingRatio
     - OPTIONAL / real64
     - 1.0e2
     - Upper bound for learned tangent diagonals as a fraction of the initial tangent scale.
   * - stressControlTangentFloorRatio
     - OPTIONAL / real64
     - 1.0e-3
     - Lower bound for learned tangent diagonals as a fraction of the initial tangent scale.
   * - stressControlWaveCutbackRatio
     - OPTIONAL / real64
     - 1.0e-1
     - Raw-minus-filtered stress ratio above which tangent adaptation remains frozen during re-engagement.
   * - stressControlWaveWarnRatio
     - OPTIONAL / real64
     - 5.0e-2
     - Raw-minus-filtered stress ratio that warns about inertial/wave dominated stress control.
   * - stressTable
     - OPTIONAL / array2d< real64 >
     - not explicitly set
     - Array that stores the time-depended grid aligned stresses
   * - stressTableInterpType
     - OPTIONAL / string (enum)
     - Linear
     - The type of stress table interpolation. Options are 0 (linear), 1 (cosine), 2 (quintic polynomial).
   * - subdivideParticles
     - OPTIONAL / int
     - -1
     - Particle subdivision policy: -1 selects the automatic default, 0 disables subdivision, and 1 enables it. Automatic subdivision is enabled only for vpHencky CPDI domain scaling with linear CPDI shape functions; it remains disabled for B-spline particles and deformation-gradient reset alone. Explicit 1 also splits native tets/wedges by current diameter using 0.99998 times the minimum active grid spacing. Native splits preserve shape and are bounded to four children in plane strain or eight in 3D.
   * - surfaceDetection
     - OPTIONAL / int
     - 0
     - Flag for automatic surface detection on the 1st cycle
   * - surfaceNormalAndPositionDamageThreshold
     - OPTIONAL / real64
     - 0.9999
     - Sets the damage threshold at which explicit particle surface normals, positions, and tractions are disabled
   * - surfaceQualityThreshold
     - OPTIONAL / real64
     - 0.0
     - Minimum mapping-normal tensor quality required for same-group DFG separability. Valid range is [0,1]; 0 disables the quality restriction and 1 requires a fully coherent, fully covered normal.
   * - surfaceTensionForceSign
     - OPTIONAL / real64
     - 1.0
     - Sign multiplier for the grid continuum-surface-stress force. Leave at +1; set to -1 only when validating against an opposite normal/color convention.
   * - surfaceTensionGradientTolerance
     - OPTIONAL / real64
     - 0.0
     - Lower bound for |grad color| used to regularize surface-tension normals. Units are inverse length. A non-positive value selects an automatic value.
   * - surfaceTensionPairs
     - OPTIONAL / array2d< real64 >
     - not explicitly set
     - Rows of { phaseA, phaseB, gamma }. Phases are particleGroup/contact-group ids. Use phaseB=-1 for a phase-void/free-surface pair. gamma has units of force/length.
   * - thinFeatureDFGThreshold
     - OPTIONAL / real64
     - DBL_MAX
     - Suppresses DFG separability when the squared distance between field centers is less than (thinFeatureDFGThreshold * neighborRadius)^2.
   * - timeIntegrationOption
     - OPTIONAL / string (enum)
     - ExplicitDynamic
     - Time integration method. Options are:; * ; *
   * - totalBinderVolume
     - OPTIONAL / real64
     - 0.0
     - Total volume of binder
   * - tracerCoordinates
     - OPTIONAL / array2d< real64 >
     - not explicitly set
     - Coordinates used to initialize MPM lagrangian tracers. Each row is x,y,z; each tracer follows the nearest particle selected during initialization.
   * - tracerCycleInterval
     - OPTIONAL / int
     - 0
     - Cycle interval for MPM lagrangian tracer CSV output. Use either tracerCycleInterval or tracerWriteInterval. A value of 0 means no cycle-based interval.
   * - tracerHistory
     - OPTIONAL / int
     - 0
     - Flag to enable MPM lagrangian tracer CSV output. One CSV file is written per tracer point.
   * - tracerOutputPrefix
     - OPTIONAL / string
     - "tracerPoint"
     - Filename prefix for MPM lagrangian tracer CSV files.
   * - tracerVariables
     - OPTIONAL / string_array
     - not explicitly set
     - Particle variables to write to each MPM lagrangian tracer CSV file. Supported names include particleID, mass, volume, materialType, porosity, density, damage, temperature, internalEnergy, kineticEnergy, speed, velocityX, velocityY, velocityZ, plasticStrainMagnitude, stressXX, stressYY, stressZZ, stressYZ, stressXZ, stressXY, and plasticStrainXX/YY/ZZ/YZ/XZ/XY.
   * - tracerWriteInterval
     - OPTIONAL / real64
     - 0.0
     - Time interval for MPM lagrangian tracer CSV output. Use either tracerWriteInterval or tracerCycleInterval. A value of 0 means no time-based interval.
   * - updateMethod
     - OPTIONAL / string (enum)
     - FLIP
     - Update method. Options are:; * ; *
   * - updateOrder
     - OPTIONAL / integer
     - 2
     - Order for update method, only applies to XPIC and FMPM
   * - useCrackTipDetection
     - OPTIONAL / int
     - 0
     - Set to activate crack-tip detection and distance calc
   * - useEvents
     - OPTIONAL / int
     - 0
     - Enable events
   * - useInternalForceAsFaceReaction
     - OPTIONAL / int
     - 0
     - Will use internal force component at boundary node as reaction
   * - useNodePosForArea
     - OPTIONAL / int
     - 0
     - Option to compute nodal area from grid surface normals and positions (results may be inaccurate if implicit surface normals and positiosn are used)
   * - useReferenceVectorsForParticleUpdate
     - OPTIONAL / int
     - 0
     - Flag for determining the integration scheme of particle vector attributes
   * - useSurfacePositionForContact
     - OPTIONAL / int
     - 0
     - Flag for to use particle mapped surface positions in contact calculations
   * - writeParticleData
     - OPTIONAL / int
     - 0
     - Flag to enable writing particle data to file
