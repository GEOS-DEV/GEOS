.. _DeviceKernelDataMembers:


###############################################################################
Device Kernel Data Members
###############################################################################

Store scalar parameters and LvArray views **by value** in objects copied into
device kernels. Reference data members can retain host addresses that a GPU
cannot access, even when CPU tests pass. Constructor reference parameters are
allowed if their values are copied into members:

.. code-block:: cpp

   class KernelWrapper
   {
     // BAD: copying this wrapper preserves a reference to the host parameter.
     real64 const & m_hostParameter;

     // GOOD: the parameter is copied into the wrapper.
     real64 const m_parameter;

   public:
     KernelWrapper( real64 const & parameter ) :
       m_hostParameter( parameter ),
       m_parameter( parameter )
     {}
   };

The ``check_code_rules`` CI job rejects reference members in ``KernelWrapper``,
``*Update`` and ``*Updates`` types, types named by ``KernelWrapper`` aliases,
and their known base classes and members stored by value.

Run the check from the repository root; no GEOS build or GPU is required:

.. code-block:: bash

   python3 -m pip install -r scripts/code_checks_requirements.txt
   python3 scripts/check_device_kernel_members.py

Pass file paths to check selected sources. Violations report file and line
numbers and cause a nonzero exit status.
