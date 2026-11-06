.. _Capabilities:

###############################################################################
Capabilities Overview
###############################################################################

The table below lists the physics GEOS can simulate, with the XML tag of each solver
(placed in the ``<Solvers>`` block). Each tag links to its complete list of attributes;
the last column points to the user guide page or example that explains it.
Coupled solvers combine single-physics solvers, which they reference by name
(e.g. ``flowSolverName``, ``solidSolverName``, ``wellSolverName``).

.. list-table::
   :header-rows: 1
   :widths: 18 52 30

   * - Physics
     - Solver XML tags
     - Documentation
   * - Single-phase flow
     - :ref:`SinglePhaseFVM <XML_SinglePhaseFVM>`, :ref:`SinglePhaseHybridFVM <XML_SinglePhaseHybridFVM>`
     - :ref:`SinglePhaseFlow`, :ref:`TutorialSinglePhaseFlowWithInternalMesh`
   * - Multiphase / compositional flow
     - :ref:`CompositionalMultiphaseFVM <XML_CompositionalMultiphaseFVM>`,
       :ref:`CompositionalMultiphaseHybridFVM <XML_CompositionalMultiphaseHybridFVM>`,
       :ref:`ImmiscibleMultiphaseFlow <XML_ImmiscibleMultiphaseFlow>`,
       :ref:`ReactiveCompositionalMultiphaseOBL <XML_ReactiveCompositionalMultiphaseOBL>`
     - :ref:`CompositionalMultiphaseFlow`, :ref:`ImmiscibleMultiphaseFlow`, :ref:`TutorialDeadOilBottomLayersSPE10`
   * - Thermal flow
     - Flow and poromechanics solvers with ``isThermal="1"``
     - :ref:`ValidationStudiesThermoPoromechanics`
   * - Wells
     - :ref:`WellManager <XML_WellManager>` (well definitions and controls), coupled to the reservoir by
       :ref:`SinglePhaseReservoir <XML_SinglePhaseReservoir>`,
       :ref:`CompositionalMultiphaseReservoir <XML_CompositionalMultiphaseReservoir>`
     - :ref:`CompositionalMultiphaseWell`, :ref:`TutorialDeadOilEgg`, :ref:`TutorialCO2FieldCase`
   * - Solid mechanics
     - :ref:`SolidMechanicsLagrangianFEM <XML_SolidMechanicsLagrangianFEM>`,
       :ref:`SolidMechanics_MPM <XML_SolidMechanics_MPM>` (material point method)
     - :ref:`SolidMechanicsLagrangianFEM`, :ref:`TutorialElasticity`, :ref:`ValidationViscoplasticModels`
   * - Contact and fractures
     - :ref:`SolidMechanicsLagrangeContact <XML_SolidMechanicsLagrangeContact>`,
       :ref:`SolidMechanicsLagrangeContactBubbleStab <XML_SolidMechanicsLagrangeContactBubbleStab>`,
       :ref:`SolidMechanicsAugmentedLagrangianContact <XML_SolidMechanicsAugmentedLagrangianContact>`,
       :ref:`SolidMechanicsEmbeddedFractures <XML_SolidMechanicsEmbeddedFractures>`,
       :ref:`SurfaceGenerator <XML_SurfaceGenerator>`,
       :ref:`EmbeddedSurfaceGenerator <XML_EmbeddedSurfaceGenerator>`
     - :ref:`ContactMechanicsSolver`, :ref:`ValidationStudiesFaultMechanics`
   * - Poromechanics
     - :ref:`SinglePhasePoromechanics <XML_SinglePhasePoromechanics>`,
       :ref:`MultiphasePoromechanics <XML_MultiphasePoromechanics>`,
       with wells: :ref:`SinglePhasePoromechanicsReservoir <XML_SinglePhasePoromechanicsReservoir>`,
       :ref:`SinglePhaseReservoirPoromechanics <XML_SinglePhaseReservoirPoromechanics>`,
       :ref:`MultiphasePoromechanicsReservoir <XML_MultiphasePoromechanicsReservoir>`,
       :ref:`CompositionalMultiphaseReservoirPoromechanics <XML_CompositionalMultiphaseReservoirPoromechanics>`
     - :ref:`PoroelasticSolver`, :ref:`TutorialPoroelasticity`, :ref:`ValidationStudiesPoromechanics`
   * - Poromechanics with fractures
     - :ref:`SinglePhasePoromechanicsConformingFractures <XML_SinglePhasePoromechanicsConformingFractures>`,
       :ref:`SinglePhasePoromechanicsConformingFracturesALM <XML_SinglePhasePoromechanicsConformingFracturesALM>`,
       :ref:`SinglePhasePoromechanicsConformingFracturesReservoir <XML_SinglePhasePoromechanicsConformingFracturesReservoir>`,
       :ref:`SinglePhaseReservoirPoromechanicsConformingFractures <XML_SinglePhaseReservoirPoromechanicsConformingFractures>`,
       :ref:`SinglePhaseReservoirPoromechanicsConformingFracturesALM <XML_SinglePhaseReservoirPoromechanicsConformingFracturesALM>`,
       :ref:`MultiphasePoromechanicsConformingFractures <XML_MultiphasePoromechanicsConformingFractures>`,
       :ref:`MultiphasePoromechanicsConformingFracturesALM <XML_MultiphasePoromechanicsConformingFracturesALM>`,
       :ref:`CompositionalMultiphaseReservoirPoromechanicsConformingFractures <XML_CompositionalMultiphaseReservoirPoromechanicsConformingFractures>`,
       :ref:`CompositionalMultiphaseReservoirPoromechanicsConformingFracturesALM <XML_CompositionalMultiphaseReservoirPoromechanicsConformingFracturesALM>`,
       :ref:`SinglePhasePoromechanicsEmbeddedFractures <XML_SinglePhasePoromechanicsEmbeddedFractures>`,
       :ref:`OneWayCoupledFractureFlowContactMechanics <XML_OneWayCoupledFractureFlowContactMechanics>`
     - :ref:`ValidationStudiesFaultMechanics`
   * - Hydraulic fracturing and proppant
     - :ref:`Hydrofracture <XML_Hydrofracture>`,
       :ref:`ProppantTransport <XML_ProppantTransport>`,
       :ref:`SinglePhaseProppantFVM <XML_SinglePhaseProppantFVM>`,
       :ref:`FlowProppantTransport <XML_FlowProppantTransport>`
     - :ref:`TutorialHydraulicFractureWithAdvancedXML`, :ref:`ValidationStudiesHydraulicFracture`, :ref:`ProppantTransport`
   * - Phase-field fracture
     - :ref:`PhaseFieldDamageFEM <XML_PhaseFieldDamageFEM>`,
       :ref:`PhaseFieldFracture <XML_PhaseFieldFracture>`,
       :ref:`PhaseFieldPoromechanics <XML_PhaseFieldPoromechanics>`
     - Datastructure pages only
   * - Wave propagation (seismic)
     - :ref:`AcousticSEM <XML_AcousticSEM>`, :ref:`AcousticFirstOrderSEM <XML_AcousticFirstOrderSEM>`,
       :ref:`AcousticVTISEM <XML_AcousticVTISEM>`, :ref:`AcousticDG <XML_AcousticDG>`,
       :ref:`AcousticElasticSEM <XML_AcousticElasticSEM>`, :ref:`ElasticSEM <XML_ElasticSEM>`,
       :ref:`ElasticFirstOrderSEM <XML_ElasticFirstOrderSEM>`
     - Datastructure pages only
   * - Induced seismicity
     - :ref:`SeismicityRate <XML_SeismicityRate>`,
       :ref:`ExplicitSpringSlider <XML_ExplicitSpringSlider>`, :ref:`ImplicitSpringSlider <XML_ImplicitSpringSlider>`,
       :ref:`ExplicitQuasiDynamicEQ <XML_ExplicitQuasiDynamicEQ>`, :ref:`ImplicitQuasiDynamicEQ <XML_ImplicitQuasiDynamicEQ>`
     - Datastructure pages only
   * - Reactive transport
     - :ref:`SinglePhaseReactiveTransport <XML_SinglePhaseReactiveTransport>`
     - Datastructure pages only
   * - Other
     - :ref:`LaplaceFEM <XML_LaplaceFEM>` (Laplace equation, used as a simple example and in tests)
     - :ref:`AddingNewSolver`

Supporting features used across these solvers include structured and unstructured meshes
(:ref:`Meshes`), finite-volume (TPFA, hybrid mimetic) and finite-element discretizations
(:ref:`NumericalMethodsManager`), a library of constitutive models (:ref:`Constitutive`),
standalone material-model drivers for calibration (:ref:`ConstitutiveDrivers`),
and a Python interface (:ref:`pygeosxInterface`).
