

================== ======= =========== ============================================== 
Name               Type    Default     Description                                    
================== ======= =========== ============================================== 
amgCoarseSolver    string  direct      AMG coarsest level solver/smoother type        
amgNumSweeps       integer 2           AMG smoother sweeps                            
amgSmootherType    string  gaussSeidel AMG smoother type                              
amgThreshold       real64  0           AMG strength-of-connection threshold           
iluFill            integer 0           ILU(K) fill factor                             
iluThreshold       real64  0           ILU(T) threshold factor                        
krylovAdaptiveTol  integer 0           Use Eisenstat-Walker adaptive linear tolerance 
krylovMaxIter      integer 200         Maximum iterations allowed                     
krylovTol          real64  1e-06       Relative convergence tolerance                 
krylovWeakestTol   real64  0.001       Weakest-allowed tolerance for adaptive method  
logLevel           integer 0           Log level                                      
preconditionerType string  iluk        Preconditioner type                            
solverType         string  direct      Linear solver type                             
================== ======= =========== ============================================== 


