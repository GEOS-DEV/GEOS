

============= ============ ======== ============================================================================ 
Name          Type         Default  Description                                                                  
============= ============ ======== ============================================================================ 
expression    string       required Symbolic math expression                                                     
inputVarNames string_array {}       Name of fields are input to function.                                        
inputVarScale real64_array {1}      Scaling applied to function inputs before function evaluation.               
name          string       required A name is required for any non-unique nodes                                  
variableNames string_array required List of variables in expression.  The order must match the evaluate argument 
============= ============ ======== ============================================================================ 


