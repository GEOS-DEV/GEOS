
def SaveWindowPlot( name, TYPE ):
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 1
    SaveWindowAtts.outputDirectory = "."
    SaveWindowAtts.fileName = name
    SaveWindowAtts.family = 0
    if TYPE==0:
        SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
    if TYPE==1:
        SaveWindowAtts.format = SaveWindowAtts.PLY
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()


def SaveTimeHistory( name, array, header ):
    with open(name, 'w') as filehandle:
        filehandle.write("%s\n" % header )
        for row in array:
            for entry in row:
                filehandle.write("%10.4g" % entry)
            filehandle.write("\n")




def RadiusScatterPlot( field ):
    SetTimeSliderState(TimeSliderGetNStates()-1)
    AddPlot("Scatter", "face_radius", 1, 0)
    ScatterAtts = ScatterAttributes()
    ScatterAtts.var1 = "face_radius"
    ScatterAtts.var1Role = ScatterAtts.Coordinate0  # Coordinate0, Coordinate1, Coordinate2, Color, None
    ScatterAtts.var1MinFlag = 0
    ScatterAtts.var1MaxFlag = 0
    ScatterAtts.var1Min = 0
    ScatterAtts.var1Max = 1
    ScatterAtts.var1Scaling = ScatterAtts.Linear  # Linear, Log, Skew
    ScatterAtts.var1SkewFactor = 1
    ScatterAtts.var2Role = ScatterAtts.Coordinate1  # Coordinate0, Coordinate1, Coordinate2, Color, None
    ScatterAtts.var2 = field
    ScatterAtts.var2MinFlag = 0
    ScatterAtts.var2MaxFlag = 0
    ScatterAtts.var2Min = 0
    ScatterAtts.var2Max = 1
    ScatterAtts.var2Scaling = ScatterAtts.Linear  # Linear, Log, Skew
    ScatterAtts.var2SkewFactor = 1
    ScatterAtts.pointSize = 0.05
    ScatterAtts.pointSizePixels = 1
    ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
    ScatterAtts.scaleCube = 0
    ScatterAtts.colorType = ScatterAtts.ColorByForegroundColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
    ScatterAtts.singleColor = (255, 0, 0, 255)
    ScatterAtts.colorTableName = "Default"
    ScatterAtts.invertColorTable = 0
    ScatterAtts.legendFlag = 1
    SetPlotOptions(ScatterAtts)
    AddOperator("Threshold", 0)
    ThresholdAtts = ThresholdAttributes()
    ThresholdAtts.outputMeshType = 0
    ThresholdAtts.listedVarNames = ("Fracture_ElementFields/elementAperture", "Fracture_ElementFields/ghostRank")
    ThresholdAtts.zonePortions = (1, 1)
    ThresholdAtts.lowerBounds = (1.0e-6, -1e+37)
    ThresholdAtts.upperBounds = (1e+37, -1)
    ThresholdAtts.defaultVarName = "face_radius"
    ThresholdAtts.defaultVarIsScalar = 1
    SetOperatorOptions(ThresholdAtts, 0)
    DrawPlots()


def TimehistQuery( state, timehist ):
    print "processing state ", state
    SetTimeSliderState(state)
    SetQueryFloatFormat("%g")
    time = Query("Time")[:-1].split(' ')[-1]
    injectionPressure = ZonePick(coord=(0.5, 0.5, 0.0), vars=("Fracture_ElementFields/pressure"))['Fracture_ElementFields/pressure']
    injectionAperture = ZonePick(coord=(0.5, 0.5, 0.0), vars=("Fracture_ElementFields/elementAperture"))['Fracture_ElementFields/elementAperture']
    fractureArea = Query("Variable Sum").split(' ')[-1]
    timehist.append([float(time), float(injectionPressure), float(injectionAperture), 4*float(fractureArea)])


    

if len(sys.argv) < 4:
    sys.exit('Usage: %s path/To/database database outputroot' % sys.argv[0])

if not os.path.exists(sys.argv[1]):
    sys.exit('ERROR: path %s was not found!' % sys.argv[1])

database = "localhost:" + sys.argv[1] + "/" + sys.argv[2]
outputroot = sys.argv[3]

print "using " + database

#try:
#    AddArgument("-np8")
#    AddArgument("-ppdebug")    
#    AddArgument("-nowin")   # no windows
#    Launch()
#except VisItException:
#    print 'Visit opened, no need to relaunch

OpenDatabase(database)
print "opened database"


#             1234567890    1234567890    1234567890    1234567890
header   = [['      time', '  pressure', '  aperture', '      area']]
timehist = []

AddPlot("Pseudocolor", "Fracture_ElementFields/elementArea", 1, 1)

AddOperator("Threshold", 1)
ThresholdAtts = ThresholdAttributes()
ThresholdAtts.outputMeshType = 0
ThresholdAtts.boundsInputType = 0
ThresholdAtts.listedVarNames = ("Fracture_ElementFields/ghostRank")
ThresholdAtts.zonePortions = (1)
ThresholdAtts.lowerBounds = (-1e+37)
ThresholdAtts.upperBounds = (-1)
ThresholdAtts.defaultVarName = "Fracture_ElementFields/elementArea"
ThresholdAtts.defaultVarIsScalar = 1
ThresholdAtts.boundsRange = ("-1e+37:-1")
SetOperatorOptions(ThresholdAtts, 1)
DrawPlots()

TimehistQuery( 1, timehist )

ThresholdAtts = ThresholdAttributes()
ThresholdAtts.outputMeshType = 0
ThresholdAtts.boundsInputType = 0
ThresholdAtts.listedVarNames = ("Fracture_ElementFields/ghostRank", "Fracture_ElementFields/elementAperture")
ThresholdAtts.zonePortions = (1, 1)
ThresholdAtts.lowerBounds = (-1e+37, 0.06e-3)
ThresholdAtts.upperBounds = (-1, 1e+37)
ThresholdAtts.defaultVarName = "Fracture_ElementFields/elementArea"
ThresholdAtts.defaultVarIsScalar = 1
ThresholdAtts.boundsRange = ("-1e+37:-1", "0.06e-3:1e+37")
SetOperatorOptions(ThresholdAtts, 1)
DrawPlots()



for state in range(2,TimeSliderGetNStates()):
    TimehistQuery( state, timehist )
#    print "processing state ", state
#    SetTimeSliderState(state)
#    SetQueryFloatFormat("%g")
#    time = Query("Time")[:-1].split(' ')[-1]
#    injectionPressure = ZonePick(coord=(0.5, 0.5, 0.0), vars=("Fracture_ElementFields/pressure"))['Fracture_ElementFields/pressure']
#    injectionAperture = ZonePick(coord=(0.5, 0.5, 0.0), vars=("Fracture_ElementFields/elementAperture"))['Fracture_ElementFields/elementAperture']
#    fractureArea = Query("Variable Sum").split(' ')[-1]
#    timehist.append([float(time), float(injectionPressure), float(injectionAperture), float(fractureArea)])

SaveTimeHistory(outputroot + '_timehist.txt', timehist, header)
print timehist


DefineScalarExpression("node_radius", "cylindrical_radius(Fracture_Fluid)")
DefineScalarExpression("face_radius", "recenter(node_radius, \"zonal\")")

DeleteActivePlots()
RadiusScatterPlot("Fracture_ElementFields/elementAperture")
SaveWindowPlot( outputroot+"_aperture", 1 )
DeleteActivePlots()
RadiusScatterPlot("Fracture_ElementFields/pressure")
SaveWindowPlot( outputroot+"_pressure", 1 )
DeleteActivePlots()

quit()
