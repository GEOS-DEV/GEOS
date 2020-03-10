
import numpy as np


def WriteVTK(previewName, fractureCenter, fractureLength, fractureHeight, fractureAngle, fracturePlaneID, numIntersections, header='dfn'):
  """
  Write a VTK file for previewing the generated DFN

  @arg previewName Name of the origianl XML file
  @arg fractureCenter Nx3 List of nodeset centers
  @arg fractureLength List of nodeset lengths
  @arg fractureHeight List of nodeset heights
  @arg fractureAngle List of nodeset angles
  @arg numIntersections Number of fracture intersections
  @arg header Additional header information
  """

  print ('Writing preview file')
  myfile = open(previewName, 'w')
  myfile.write("# vtk DataFile Version 2.0" + '\n')
  myfile.write(header + '\n')
  myfile.write("ASCII" + '\n')
  myfile.write("DATASET UNSTRUCTURED_GRID" + '\n')
  myfile.write("POINTS" + " " + str(len(fractureCenter) * 4) + " float" + '\n')

  for ii in range(0, len(fractureCenter)):
    myfile.write(str(fractureCenter[ii][0] - np.cos(fractureAngle[ii]) * fractureLength[ii] * 0.5) + ' ')
    myfile.write(str(fractureCenter[ii][1] - np.sin(fractureAngle[ii]) * fractureLength[ii] * 0.5) + ' ')
    myfile.write(str(fractureCenter[ii][2] - fractureHeight[ii] * 0.5) + '\n')

    myfile.write(str(fractureCenter[ii][0] + np.cos(fractureAngle[ii]) * fractureLength[ii] * 0.5) + ' ')
    myfile.write(str(fractureCenter[ii][1] + np.sin(fractureAngle[ii]) * fractureLength[ii] * 0.5) + ' ')
    myfile.write(str(fractureCenter[ii][2] - fractureHeight[ii] * 0.5) + '\n')

    myfile.write(str(fractureCenter[ii][0] + np.cos(fractureAngle[ii]) * fractureLength[ii] * 0.5) + ' ')
    myfile.write(str(fractureCenter[ii][1] + np.sin(fractureAngle[ii]) * fractureLength[ii] * 0.5) + ' ')
    myfile.write(str(fractureCenter[ii][2] + fractureHeight[ii] * 0.5) + '\n')

    myfile.write(str(fractureCenter[ii][0] - np.cos(fractureAngle[ii]) * fractureLength[ii] * 0.5) + ' ')
    myfile.write(str(fractureCenter[ii][1] - np.sin(fractureAngle[ii]) * fractureLength[ii] * 0.5) + ' ')
    myfile.write(str(fractureCenter[ii][2] + fractureHeight[ii] * 0.5) + '\n')

  myfile.write("CELLS" + " " + str(len(fractureCenter)) + " " + str(len(fractureCenter) * 5) + '\n')

  for i in range(0, len(fractureCenter)):
    myfile.write(str(4) + " " + str(i * 4) + " " + str(i * 4 + 1) + " " + str(i * 4 + 2) + " " + str(i * 4 + 3) + '\n')

  myfile.write("CELL_TYPES" + " " + str(len(fractureCenter)) + '\n')
  for i in range(0, len(fractureCenter)):
    myfile.write('9 ')
    if (i % 10 == 9):
      myfile.write('\n')

  myfile.write('\n' + "CELL_DATA" + " " + str(len(fractureLength)) + '\n')

  myfile.write("SCALARS Length float 1" + '\n')
  myfile.write("LOOKUP_TABLE default" + '\n')
  for i in range(0, len(fractureLength)):
    myfile.write(str(fractureLength[ii]) + ' ')
    if (i % 10 == 9):
      myfile.write('\n')

  myfile.write('\n' + "SCALARS planeID int 1" + '\n')
  myfile.write("LOOKUP_TABLE default" + '\n')
  for i in range(0, len(fracturePlaneID)):
    myfile.write(str(fracturePlaneID[ii]) + ' ')
    if (i % 10 == 9):
      myfile.write('\n')

  myfile.write('\n' + "SCALARS numIntersect int 1" + '\n')
  myfile.write("LOOKUP_TABLE default" + '\n')
  for i in range(0, len(numIntersections)):
    myfile.write(str(numIntersections[ii]) + ' ')
    if (i % 10 == 9):
      myfile.write('\n')


