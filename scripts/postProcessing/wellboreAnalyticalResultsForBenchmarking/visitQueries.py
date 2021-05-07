import numpy as np


# Rotate a vector in local coodinate of an inclined borehole to the global coordinate
def vectorRotation(x,y,z,phi_x,phi_z):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x),0.],[-np.sin(phi_x), np.cos(phi_x),0.],[0.,0.,1.]])
    rotz = np.array([[np.cos(phi_z),0., np.sin(phi_z)],[0.,1.,0.],[-np.sin(phi_z),0., np.cos(phi_z)]])
    
    localCoord = np.array([x,y,z])
    return np.dot( rotz, np.dot( rotx,localCoord ) )

start = vectorRotation( 0.0, 0, 0, -90./180*3.1416 , 70./180*3.1416 )
end   = vectorRotation( 1.0, 0, 0, -90./180*3.1416 , 70./180*3.1416 )

print(start)
print(end)



OpenDatabase("localhost:/scratchrd/SCR/GEOS/GEOSX/SyTuan/dev/deviatedWellborePoroelastic/siloFiles/plot_* database", 0)
AddPlot("Pseudocolor", "Omega_ElementFields/pressure", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=( end[0], end[1], end[2] ), num_samples=50, start_point=( start[0], start[1], start[2] ), use_sampling=0, vars=("Omega_ElementFields/pressure"))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)


SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "pressure"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()

SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()

AddPlot("Pseudocolor", "Omega_NodalFields/TotalDisplacement_1", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=( end[0], end[1], end[2] ), num_samples=50, start_point=( start[0], start[1], start[2] ), use_sampling=0, vars=("Omega_NodalFields/TotalDisplacement_1"))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)


SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "radialDisplacement"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()


SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()

AddPlot("Pseudocolor", "Omega_Solid_MaterialFields/stress_22", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=( end[0], end[1], end[2] ), num_samples=50, start_point=( start[0], start[1], start[2] ), use_sampling=0, vars=("Omega_Solid_MaterialFields/stress_22"))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)


SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "radialStress"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()

SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()

AddPlot("Pseudocolor", "Omega_Solid_MaterialFields/stress_11", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
Query("Lineout", end_point=( end[0], end[1], end[2] ), num_samples=50, start_point=( start[0], start[1], start[2] ), use_sampling=0, vars=("Omega_Solid_MaterialFields/stress_11"))
SetActiveWindow(2)
SetTimeSliderState(TimeSliderGetNStates()-1)


SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "tangentStress"
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.CURVE  
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
DeleteWindow()


SetActiveWindow(1)
SetActivePlots(0)
DeleteActivePlots()

