.. _FieldSpecificationBase:

####################################################
FieldSpecification Options
####################################################

Overview
======================

`FieldSpecification` entries live under the `FieldSpecifications` block and are used to define initial and boundary conditions for fields registered on mesh managers (nodes, faces, or element subregions).
This page documents options added for per-component values and region-based targeting.

Region-based targetting
===============================

Use `regionNames` to target one or more element regions without specifying an explicit `objectPath`.
When `regionNames` is provided, GEOS builds the `objectPath` internally and applies the specification to all regions.
You can also specify subregions in `regionNames` by using a syntax similar to the objectPath (cf. examples).

.. note::
    `regionNames` and `objectPath` are mutually exclusive. Use one or the other.

Non-scalar values
======================

The `scale` and `functionName` attributes accept either a single value or a list of values surrounded by braces.

- `scale="1.0"` and `scale="{ 1.0 }"` are equivalent.
- `scale={ a, b, c }` applies one scale factor per component.
- `functionName` follows the same convention. For convenience, if a single function name is set it will be applied to every component.
- When `scale` has more than one entry, `component` must not be set.

Examples
======================

The snippet below shows a specification using scalar values in `scale` and `functionName`.

.. code-block:: xml

    <FieldSpecifications>
        <FieldSpecification
            name="perm"
            initialCondition="1"
            objectPath="ElementRegions/region1/block1"
            setNames="{ all }"
            fieldName="permeability"
            scale="1.0e-22"
            functionName="func"
        />
    </FieldSpecifications>
The snippet below shows a specification using non-scalar values in `scale` and `functionName`.

.. code-block:: xml

    <FieldSpecifications>
        <FieldSpecification
            name="perm"
            initialCondition="1"
            objectPath="ElementRegions/region1/block1"
            setNames="{ all }"
            fieldName="permeability"
            scale="{ 1.0e-22, 2.0e-22, 3.0e-22 }"
            functionName="{ func1, func2, func3 }"
        />
    </FieldSpecifications>

The snippet below shows a specification using a single `functionName` applied to every components.

.. code-block:: xml

    <FieldSpecifications>
        <FieldSpecification
            name="perm"
            initialCondition="1"
            objectPath="ElementRegions/region1/block1"
            setNames="{ all }"
            fieldName="permeability"
            scale="{ 1.0e-22, 2.0e-22, 3.0e-22 }"
            functionName="func"
        />
    </FieldSpecifications>

The snippet below shows a specification using `regionNames` on a single region.

.. code-block:: xml

    <FieldSpecifications>
        <FieldSpecification
            name="perm"
            initialCondition="1"
            regionNames="region1"
            setNames="{ all }"
            fieldName="permeability"
            scale="{ 1.0e-22, 2.0e-22, 3.0e-22 }"
        />
    </FieldSpecifications>

The snippet below shows a specification using `regionNames` on a multiple subregions.

.. code-block:: xml

    <FieldSpecifications>
        <FieldSpecification
            name="perm"
            initialCondition="1"
            regionNames="{ region1/block1, region1/block2 }"
            setNames="{ all }"
            fieldName="permeability"
            scale="{ 1.0e-22, 2.0e-22, 3.0e-22 }"
        />
    </FieldSpecifications>