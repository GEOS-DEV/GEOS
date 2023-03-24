.. _InternalNodesetNames:

####################################################
List of internal nodeset names
####################################################

Nodeset name is required for defining initial and boundary conditions. When using the ``InternalMesh`` generator, following nodesets are automatically generated:

* ``all``: this nodeset contains all nodes of the generated mesh.
* ``xneg``: this nodeset contains all nodes on the face that corresponds to the smallest x-coordinate.
* ``xpos``: this nodeset contains all nodes on the face that corresponds to the highest x-coordinate.
* ``yneg``: this nodeset contains all nodes on the face that corresponds to the smallest y-coordinate.
* ``ypos``: this nodeset contains all nodes on the face that corresponds to the highest y-coordinate.
* ``zneg``: this nodeset contains all nodes on the face that corresponds to the smallest z-coordinate.
* ``zpos``: this nodeset contains all nodes on the face that corresponds to the highest z-coordinate.

When using the ``InternalWellbore`` generator, following additional nodesets are automatically generated:

* ``rneg``: this nodeset contains all nodes on the face that corresponds to the smallest radial coordinate.
* ``rpos``: this nodeset contains all nodes on the face that corresponds to the highest radial coordinate.
* ``tneg``: this nodeset contains all nodes on the face that corresponds to the smallest tangent coordinate.
* ``tpos``: this nodeset contains all nodes on the face that corresponds to the highest tangent coordinate.
* ``rCasingCementInterface``: this nodeset contains all nodes on the interface between the casing and the cement layers of a cased wellbore.
* ``rCementRockInterface``: this nodeset contains all nodes on the interface between the cement layer and rock formation around a cased wellbore.

An example of using pre-defined nodeset for declaring the boundary condition is given below:

.. code-block:: xml

   <FieldSpecifications>		
	     
      <FieldSpecification
      name="itConstraint_x"
      objectPath="nodeManager"
      fieldName="totalDisplacement"
      component="0"
      scale="0.0"
      setNames="{ rCementRockInterface }"/>

   </FieldSpecifications>


