.. _CodeRules:

##########
Code Rules
##########

.. _CodeRules_Logging:

*******
Logging
*******

Rules
=====

Console Output
--------------

**Never use** ``std::cout`` for logging. Always use the GEOS logging macros (``GEOS_LOG*``).

.. code-block:: cpp

   // Bad
   std::cout << "Processing mesh..." << std::endl;
   printf("Value: %d\n", value);

   // Good
   GEOS_LOG_RANK_0("Processing mesh...");
   GEOS_LOG_RANK_0("Value: " << value);

Log Levels
----------

Use appropriate log levels to avoid spam. Choose between ``logLevel`` for verbose information and ``logInfo`` for important messages.

.. code-block:: cpp

   GEOS_LOG_LEVEL(logInfo::logOutput, "Simulation started");
   GEOS_LOG_LEVEL(logLevel, "Iteration " << iter << " completed");

GPU Logging
-----------

.. warning::
   Logging from GPU kernels should **only be used for debugging purposes** and must **never be committed to the develop branch**.

Structured Data Output
======================

For structured data output (tables, CSV files), use the appropriate formatters:

- ``TableFormatter`` for console table output
- ``writeCSV`` for CSV file generation

.. code-block:: cpp

   TableFormatter table;
   table.addRow("Column1", "Column2", "Column3");
   // ...

.. _CodeRules_ErrorManagement:

*****************
Error Management
*****************

Errors
======

When to Use
-----------

Use ``GEOS_ERROR`` for unrecoverable errors that should terminate execution immediately.

.. code-block:: cpp

   if (invalidCondition)
   {
     GEOS_ERROR("Invalid configuration detected");
   }

Error Flag (``-e``)
-------------------

The ``-e`` flag controls error handling behavior in GEOS. [TODO: Add details about -e flag usage]

Kernel Restrictions
-------------------

.. warning::
   **Never use** ``GEOS_ERROR`` inside GPU kernels. Errors in kernels will cause undefined behavior.

GPU Usage
---------

.. warning::
   Error checking in GPU code should **only be used for debugging purposes** and must **never be committed to the develop branch**.

Exceptions
==========

When to Throw
-------------

Use exceptions (``throw``) for exceptional conditions that may be recoverable at a higher level.

.. code-block:: cpp

   if (fileNotFound)
   {
     throw std::runtime_error("Configuration file not found");
   }

Exception Handling Policy
-------------------------

.. important::
   Exceptions should **never be caught and suppressed**. If you catch an exception, either:
   
   - Re-throw it after cleanup
   - Convert it to a GEOS_ERROR if truly unrecoverable

.. code-block:: cpp


Warnings
========

When to Use
-----------

Use ``GEOS_WARNING`` for problematic conditions that don't prevent execution but should be brought to the user's attention.

.. code-block:: cpp

   if (suboptimalCondition)
   {
     GEOS_WARNING("Performance may be degraded due to...");
   }

GPU Warnings
------------

.. warning::
   Warnings from GPU code should **only be used for debugging purposes** and must **never be committed to the develop branch**.

Contextualization
=================

DataContext
-----------

``DataContext`` provides contextual information for error messages, making debugging easier.

When to Add Context
^^^^^^^^^^^^^^^^^^^

Add ``DataContext`` when:

- Processing multiple objects in a loop
- Dealing with hierarchical data structures
- Error location may not be obvious from the call stack

Syntax
^^^^^^

.. code-block:: cpp

   DataContext context;
   context.addData("meshName", mesh.getName());
   context.addData("elementIndex", elemIdx);
   
   GEOS_ERROR_WITH_CONTEXT("Invalid element configuration", context);

.. _CodeRules_Precision:

Z
***********
Performance
***********

General Rules
=============

Profiling Requirements
----------------------

.. important::
   When modifying computational kernels, **always profile** your changes using Caliper on the following sample cases:
   
   - [TODO: List specific sample cases]
   - [TODO: Add profiling procedure]

Performance profiling ensures that optimizations actually improve performance and don't introduce regressions.

.. TODO: Add additional performance rules

.. _CodeRules_Contribution:

************
Contribution
************

Documentation
=============

Doxygen Wrappers
----------------

All public classes and functions must have Doxygen documentation.

.. code-block:: cpp

   /**
    * @brief Brief description of the function
    * @param[in] input Input parameter description
    * @param[out] output Output parameter description
    * @return Description of return value
    */
   real64 computeValue(real64 input, real64 & output);

RST Documentation
-----------------

User-facing features must be documented in RST format in the appropriate documentation section.

Testing
=======

Unit Tests
----------

.. important::