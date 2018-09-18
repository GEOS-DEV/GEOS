
from pylab import *
import geos_inputs


# Setup grid
x = linspace(-10, 10, 10)
y = linspace(-10, 10, 10)
z = linspace(-10, 10, 10)
X, Y, Z = meshgrid(x, y, z, indexing='ij')


properties = {'fx': X,
              'fy': Y,
              'fz': Z,
              'r': sqrt(X**2 + Y**2 + Z**2),
              'poly': X + Y**2 + Z**3}

geos_inputs.writeGEOSTable((x, y, z), properties)

