

================== ======= =================== =============================================================================================== 
Name               Type    Default             Description                                                                                     
================== ======= =================== =============================================================================================== 
slaveDirectory     string                      slave directory path                                                                            
parallelThreads    integer 1                   Number of plot files.                                                                           
outputPath         string  required            Path at which the geosx to chombo file will be written.                                         
beginCycle         real64  required            Cycle at which the coupling will commence.                                                      
inputPath          string  /INVALID_INPUT_PATH Path at which the chombo to geosx file will be written.                                         
waitForInput       integer required            True iff geosx should wait for chombo to write out a file. When true the inputPath must be set. 
useChomboPressures integer 0                   True iff geosx should use the pressures chombo writes out.                                      
name               string  required            A name is required for any non-unique nodes                                                     
================== ======= =================== =============================================================================================== 


