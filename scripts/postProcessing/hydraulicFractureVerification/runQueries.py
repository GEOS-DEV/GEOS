
import os
import sys
sys.path.append('/usr/gapps/visit/2.7.2/linux-x86_64/lib/site-packages/')
#import visit
from visit import *



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







if len(sys.argv) < 4:
    sys.exit('Usage: %s path/To/database database outputroot' % sys.argv[0])

if not os.path.exists(sys.argv[1]):
    sys.exit('ERROR: path %s was not found!' % sys.argv[1])

database = "localhost:" + sys.argv[1] + "/" + sys.argv[2]
outputroot = sys.argv[3]
print "using " + database

try:
#    AddArgument("-np8")
#    AddArgument("-ppdebug")    
    AddArgument("-nowin")   # no windows
    Launch()
except VisItException:
    print 'Visit opened, no need to relaunch'
    


OpenDatabase(database)
print "opened database"

AddPlot("Mesh", "face_mesh", 1, 1)
AddPlot("Pseudocolor", "FaceFields/Aperture", 1, 1)
SetActivePlots((0, 1))
AddOperator("Threshold", 1)
ThresholdAtts = ThresholdAttributes()
ThresholdAtts.outputMeshType = 0
ThresholdAtts.listedVarNames = ("FaceFields/flowFaceType", "FaceFields/ghostRank")
ThresholdAtts.zonePortions = (1, 1)
ThresholdAtts.lowerBounds = (1, -1e+37)
ThresholdAtts.upperBounds = (1e+37, -1)
ThresholdAtts.defaultVarName = "default"
ThresholdAtts.defaultVarIsScalar = 0
SetOperatorOptions(ThresholdAtts, 1)
DrawPlots()


SetActivePlots(0)
QueryOverTimeAtts = GetQueryOverTimeAttributes()
QueryOverTimeAtts.timeType = 1
SetQueryOverTimeAttributes(QueryOverTimeAtts)

print "running area query"
SetQueryFloatFormat("%g")
#QueryOverTime("3D surface area", end_time=200, start_time=0, stride=1)
QueryOverTime("3D surface area")
SetActiveWindow(2)
SaveWindowPlot( outputroot+"_area", 0 )
DeleteWindow()




print "running aperture at wellbore query"
PickByZone(curve_plot_type=0, do_time=1, domain=2, element=200, preserve_coord=0, vars=("FaceFields/Aperture",))
SetActiveWindow(2)
SaveWindowPlot( outputroot+"_aper0", 0 )
DeleteWindow()

SetActivePlots((0, 1))
DeleteActivePlots()

DefineScalarExpression("node_radius", "cylindrical_radius(face_mesh)")
DefineScalarExpression("face_radius", "recenter(node_radius, \"zonal\")")



SetTimeSliderState(TimeSliderGetNStates()-1)
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
ScatterAtts.var2 = "FaceFields/Aperture"
ScatterAtts.var2MinFlag = 0
ScatterAtts.var2MaxFlag = 0
ScatterAtts.var2Min = 0
ScatterAtts.var2Max = 1
ScatterAtts.var2Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var2SkewFactor = 1
ScatterAtts.pointSize = 0.05
ScatterAtts.pointSizePixels = 1
ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
ScatterAtts.scaleCube = 1
ScatterAtts.colorType = ScatterAtts.ColorByForegroundColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (255, 0, 0, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 1
SetDefaultPlotOptions(ScatterAtts)
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
ScatterAtts.var2 = "FaceFields/Aperture"
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
ThresholdAtts.listedVarNames = ("FaceFields/flowFaceType", "FaceFields/ghostRank")
ThresholdAtts.zonePortions = (1, 1)
ThresholdAtts.lowerBounds = (1, -1e+37)
ThresholdAtts.upperBounds = (1e+37, -1)
ThresholdAtts.defaultVarName = "face_radius"
ThresholdAtts.defaultVarIsScalar = 1
SetOperatorOptions(ThresholdAtts, 0)
DrawPlots()


SaveWindowPlot( outputroot+"_aperture", 1 )
DeleteActivePlots()

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
ScatterAtts.var2 = "FaceFields/Pressure"
ScatterAtts.var2MinFlag = 0
ScatterAtts.var2MaxFlag = 0
ScatterAtts.var2Min = 0
ScatterAtts.var2Max = 1
ScatterAtts.var2Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var2SkewFactor = 1
ScatterAtts.pointSize = 0.05
ScatterAtts.pointSizePixels = 1
ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
ScatterAtts.scaleCube = 1
ScatterAtts.colorType = ScatterAtts.ColorByForegroundColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (255, 0, 0, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 1
SetDefaultPlotOptions(ScatterAtts)
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
ScatterAtts.var2 = "FaceFields/Aperture"
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
ThresholdAtts.listedVarNames = ("FaceFields/flowFaceType", "FaceFields/ghostRank")
ThresholdAtts.zonePortions = (1, 1)
ThresholdAtts.lowerBounds = (1, -1e+37)
ThresholdAtts.upperBounds = (1e+37, -1)
ThresholdAtts.defaultVarName = "face_radius"
ThresholdAtts.defaultVarIsScalar = 1
SetOperatorOptions(ThresholdAtts, 0)
DrawPlots()


SaveWindowPlot( outputroot+"_pressure", 1 )
DeleteActivePlots()
