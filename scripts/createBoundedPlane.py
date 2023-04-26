import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from geosx_xml_tools import xml_formatter as xmlFormatter


def arrayToXMLArray( array ):
    stringArray = '{%s}' % ', '.join(['%1.10f' % (x) for x in array])
    return stringArray

class BoundedPlane:
    
    def __init__(self, name, strikeAngle, dipAngle, origin) -> None:
        self.name = name
        self.normal = np.array([0.0, 0.0, 0.0], dtype=float )
        self.origin = np.array( origin, dtype=float )
        self.lengthVector = np.array([0.0, 0.0, 0.0], dtype=float )
        self.widthVector = np.array([0.0, 0.0, 0.0], dtype=float )
        self.dimensions = np.array([0.0, 0.0], dtype=float )
        self.strikeAngle = strikeAngle
        self.dipAngle = dipAngle

        self.compute_normal()
        self.compute_tangentVectors()

    def compute_tangentVectors( self ) -> None:   
        self.lengthVector = np.array( np.cross( self.normal, np.array([0.0,1.0,0.0]) ), dtype=float )
        self.lengthVector = self.lengthVector / np.linalg.norm(self.lengthVector) 
        self.widthVector = np.array( np.cross( self.normal, self.lengthVector ), dtype=float )  
        self.widthVector = self.widthVector / np.linalg.norm(self.widthVector)

        assert np.allclose([1,1,1], [self.normal.dot(self.normal), self.lengthVector.dot(self.lengthVector), self.widthVector.dot(self.widthVector)])
        assert np.allclose([0,0,0], [self.normal.dot(self.lengthVector), self.lengthVector.dot(self.widthVector), self.normal.dot(self.widthVector)])

    def compute_normal(self):
        '''
        Function to compute the plane unit normal vector given the strike and dip angles
        (in degrees). The strike should lie between 0 and 360 (negative ok) and the
        dip is restricted to lie between 0 and 90 measured in the direction such that 
        when you look in the strike direction, the fault dips to your right.
        '''
        deg_to_rad = np.pi/180
        strike = self.strikeAngle * deg_to_rad
        dip = self.dipAngle * deg_to_rad
        self.normal[0] = -np.sin(dip)*np.sin(strike) # x component
        self.normal[1] =  np.sin(dip)*np.cos(strike) # y component
        self.normal[2] = -np.cos(dip)                # z component


    def print_XML_block(self):
        
        root = ET.Element("geos")
        element = ET.SubElement(root, "BoundedPlane")
        element.set( "name" , self.name )
        element.set( "normal", arrayToXMLArray( self.normal ) )
        element.set( "origin", arrayToXMLArray( self.origin ) )
        element.set( "lengthVector", arrayToXMLArray( self.lengthVector ) )
        element.set( "widthVector", arrayToXMLArray( self.widthVector ) )
        element.set( "dimensions", arrayToXMLArray( self.dimensions ) )

        from xml.dom import minidom

        xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(indent="   ")
        fileName = self.name + "_" + str(self.strikeAngle)+ "_" + str(self.dipAngle) + ".xml"
        f = open( fileName, "w" )
        f.write( xmlstr )
        f.close()
        xmlFormatter.format_file( fileName )

def plotPlane(xx, yy, plane):
    # a plane is a*x+b*y+c*z+d=0
    # [a,b,c] is the normal. Thus, we have to calculate
    # d and we're set
    d = -plane.origin.dot(plane.normal)

    # calculate corresponding z
    z = (-plane.normal[0] * xx - plane.normal[1] * yy - d) * 1. /plane.normal[2]

    # plot the surface
    title = str(plane.strikeAngle) + " -  " + str(plane.dipAngle)
    fig = plt.figure( title )
    plt3d = fig.add_subplot(111, projection='3d')
    plt3d.plot_surface(xx, yy, z)
    plt3d.set_zlim3d(0, 140)

def plot(planes):

     # create x,y
    xx, yy = np.meshgrid(range(146), range(75))

    for plane in planes:
        plotPlane(xx, yy, plane)
    
                    
def main():
    origin = [40, 37.5, 10]
    planes = []
    planes.append(BoundedPlane("FracturePlane", -20, 60, origin))
    plot(planes)
    for plane in planes:
        plane.print_XML_block()
    plt.show()
    
    
if __name__ == "__main__":
    main()

    