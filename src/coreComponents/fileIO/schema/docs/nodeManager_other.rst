

=========================== ====================================================================================================================== ================================================= ============= 
Name                        Type                                                                                                                   Description                                       Registered By 
=========================== ====================================================================================================================== ================================================= ============= 
localToGlobalMap            globalIndex_array                                                                                                      (no description available)                                      
globalToLocalMap            map< long long, long, less< long long >, allocator< pair< long long const, long > > >                                  (no description available)                                      
isExternal                  integer_array                                                                                                          (no description available)                                      
ghostRank                   integer_array                                                                                                          (no description available)                                      
domainBoundaryIndicator     integer_array                                                                                                          (no description available)                                      
ReferencePosition           r1_array                                                                                                               (no description available)                                      
edgeList                    InterObjectRelation< Array< SortedArray< long, long >, 1, long, ChaiVector< SortedArray< long, long > > > >            (no description available)                                      
faceList                    InterObjectRelation< Array< SortedArray< long, long >, 1, long, ChaiVector< SortedArray< long, long > > > >            (no description available)                                      
elemRegionList              Array< Array< long, 1, long, ChaiVector< long > >, 1, long, ChaiVector< Array< long, 1, long, ChaiVector< long > > > > (no description available)                                      
elemSubRegionList           Array< Array< long, 1, long, ChaiVector< long > >, 1, long, ChaiVector< Array< long, 1, long, ChaiVector< long > > > > (no description available)                                      
elemList                    Array< Array< long, 1, long, ChaiVector< long > >, 1, long, ChaiVector< Array< long, 1, long, ChaiVector< long > > > > (no description available)                                      
velocityTilde               r1_array                                                                                                               (no description available)                                      
uhatTilde                   r1_array                                                                                                               (no description available)                                      
TotalDisplacement           r1_array                                                                                                               The total displacement vector.                    [':ref:`DATASTRUCTURE_SolidMechanicsLagrangianSSLEblah`', ':ref:`DATASTRUCTURE_SolidMechanics_LagrangianFEM`'] 
IncrementalDisplacement     r1_array                                                                                                               (no description available)                                      
Velocity                    r1_array                                                                                                               (no description available)                                      
Acceleration                r1_array                                                                                                               (no description available)                                      
Mass                        real64_array                                                                                                           (no description available)                                      
trilinosIndex               globalIndex_array                                                                                                      (no description available)                                      
                            real64_array                                                                                                           Primary field variable                                          
blockLocalDofNumber_Laplace globalIndex_array                                                                                                      Global DOF numbers for the primary field variable               
parentIndex                 localIndex_array                                                                                                       Parent index of node.                                           
degreeFromCrack             integer_array                                                                                                          connectivity distance from crack.                               
sets                        node                                                                                                                   :ref:`DATASTRUCTURE_sets`                                       
neighborData                node                                                                                                                   :ref:`DATASTRUCTURE_neighborData`                               
=========================== ====================================================================================================================== ================================================= ============= 


