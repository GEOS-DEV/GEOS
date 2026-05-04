.. _FieldSpecificationBase:

####################################################
FieldSpecification Options
####################################################

Overview
======================

`FieldSpecification` entries live under the `FieldSpecifications` block and are used to define initial and boundary conditions for fields registered on mesh managers (nodes, faces, or element subregions).
This page documents options added for per-component values and region-based targeting.

Region-based targetting
======================

Use `regionNames` to target one or more element regions without specifying an explicit `objectPath`.
When `regionNames` is provided, GEOS builds the `objectPath` internally and applies the specification to all regions.
You can also specify subregions in `regionNames` by using a syntax similar to the objectPath (cf. examples).

.. note::
    `regionNames` and `objectPath` are mutually exclusive. Use one or the other.

Non-scalar values
======================

Use `scales` and (optionally) `functionNames` to apply per-component values.

- `scales` replaces `scale` when you want one value per component.
- `functionNames` is optional and must be empty or have the same size as `scales`.
- When `scales` is used, do not set `component` or `functionName`.

For scalar fields, you can use `scale`, or `scales` with a single value.

Examples
======================

The snippet below shows a specification using `scales` and `functionNames`.

.. code-block:: xml

    <FieldSpecifications>
        <FieldSpecification
            name="perm"
            initialCondition="1"
            objectPath="ElementRegions/region1/block1"
            setNames="{ all }"
            fieldName="permeability"
            scales="{ 1.0e-22, 2.0e-22, 3.0e-22 }"
            functionNames="{ func1, func2, func3 }"
        />
    </FieldSpecifications>

The snippet below shows a specification using `regionNames` on a single region.

.. code-block:: xml

    <FieldSpecifications>
        <FieldSpecification
            name="perm"
            initialCondition="1"
            regionNames="{ region1 }"
            setNames="{ all }"
            fieldName="permeability"
            scales="{ 1.0e-22, 2.0e-22, 3.0e-22 }"
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
            scales="{ 1.0e-22, 2.0e-22, 3.0e-22 }"
        />
    </FieldSpecifications>