.. _LogLevelDocumentation:

Log levels documentation
========================

Add a log level
---------------

To add a log level, you must respect the following structure and add it to the appropriate ``LogLevelsInfos.hpp`` :

.. code-block:: c++

    struct MyMessage
    {
        static constexpr int getMinLogLevel() { return 2; }
        static constexpr std::string_view getDescription() { return msg; }
    };

If there is no ``LogLevelsInfos.hpp`` in the corresponding folder, you can create a ``LogLevelsInfos.hpp``.
As the system filters any doubles, we encourage anyone adding a logInfo to add it in the constructor(s) of the calling class(es), 
while ignoring any polymorphism concern (avoid adding it in base class so there will not be undesired documentation entry for other final classes).

Example of usage
----------------

To log a message with a log level, 4 macros are defined in LogLevelsInfo.hpp located in dataRepository:

* ``GEOS_LOG_LEVEL_INFO( logInfoStruct, msg )``: Output messages based on current ``Group``'s log level.
* ``GEOS_LOG_LEVEL_INFO_RANK_0( logInfoStruct, msg )``: Output messages (only on rank 0) based on current ``Group``'s log level.
* ``GEOS_LOG_LEVEL_INFO_BY_RANK( logInfoStruct, msg )``: Output messages (with one line per rank) based on current ``Group``'s log level.
* ``GEOS_LOG_LEVEL_INFO_RANK_0_NLR( logInfoStruct, msg )``: Output messages (only on rank 0) based on current ``Group``'s log level without the line return.

If you want to add a log level with a target log level belonging to another group, 
4 others macro with the same mecanism are defined in LogLevelsInfo.hpp located in dataRepository:

* ``GEOS_LOG_LEVEL_INFO_ON_GROUP( logInfoStruct, msg, group )``: Output messages based on targeted ``group``'s log level.
* ``GEOS_LOG_LEVEL_RANK_0_ON_GROUP( logInfoStruct, msg, group )``: Output messages (only on rank 0) based on targeted ``group``'s log level.
* ``GEOS_LOG_LEVEL_BY_RANK_ON_GROUP( logInfoStruct, msg, group )``: Output messages (with one line per rank) based on targeted ``group``'s log level.
* ``GEOS_LOG_LEVEL_RANK_0_NLR_ON_GROUP( logInfoStruct, msg, group )``: Output messages (only on rank 0) based on targeted ``group``'s log level without the line return.

An exemple of adding and using a log level in a ``group``:

.. code-block:: c++

    MyGroup::MyGroup( string const & name );
    {
        // To use a log level, make sure it is declared in the constructor of the group
        addLogLevel< logInfo::MyMessage >();
    }

    void MyGroup::outputMessage()
    {
        // effectively output the message taking into account the log level of the group instance
        GEOS_LOG_LEVEL_INFO( logInfo::MyMessage, "output some message" );
    }
