

=================== ======= ======= ================================================================================================================== 
Name                Type    Default Description                                                                                                        
=================== ======= ======= ================================================================================================================== 
allowNonConverged   integer 0       Allow non-converged solution to be accepted                                                                        
<<<<<<< HEAD
doLineSearch        integer 1       Line search option                                                                                                 
=======
>>>>>>> develop
dtCutIterLimit      real64  0.7     Fraction of the Max Newton iterations above which the solver asks for the time-step to be cut for the next dt.     
dtIncIterLimit      real64  0.4     Fraction of the Max Newton iterations below which the solver asks for the time-step to be doubled for the next dt. 
ilut_drop           real64  0       (no description available)                                                                                         
ilut_fill           real64  3       (no description available)                                                                                         
krylovTol           real64  1e-06   Allowable tolerance for krylov solve                                                                               
kspace              integer 0       (no description available)                                                                                         
lineSearchCutFactor real64  0.5     Line search cut factor                                                                                             
logLevel            integer 0       Log level                                                                                                          
maxIterNewton       integer 5       Maximum number of Newton iterations                                                                                
maxLineSearchCuts   integer 4       Max number of line search cuts                                                                                     
maxSubSteps         integer 10      Maximum number of time sub-steps allowed for the solver                                                            
maxTimeStepCuts     integer 2       Max number of time step cuts                                                                                       
minIterNewton       integer 1       Minimum number of Newton iterations.                                                                               
newtonTol           real64  1e-06   (no description available)                                                                                         
numKrylovIter       integer 100     Maximum number of Krylov Iterations                                                                                
scalingOption       integer 0       (no description available)                                                                                         
solverType          string          (no description available)                                                                                         
timestepCutFactor   real64  0.5     Time step cut factor                                                                                               
useBicgstab         integer 0       (no description available)                                                                                         
useDirectSolver     integer 0       (no description available)                                                                                         
useInnerSolver      integer 0       (no description available)                                                                                         
useMLPrecond        integer 0       (no description available)                                                                                         
useNewtonSolve      integer 0       (no description available)                                                                                         
=================== ======= ======= ================================================================================================================== 


