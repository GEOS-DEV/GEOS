'''
Created on 10/02/2021

@author: macpro
'''
import numpy as np

class Source:
    """A class representing a point source
    
    Attributes
    ----------
    coords :
        Coordinates of the source
    f :
        Source function (Ricker)
    
    domainFlag :
        0 the source is not in the domain, 
        1 the source is in the domain
    """
    def __init__(self, meshBox, pos, wavelet):
        """Constructor for the point source
        
        Parameters
        ----------
        meshBox :
            Numpy array containing min/max boundary coordinates of the domain
        
        pos :
            Coordinates for the point source
        
        wavelet :
            Source function (Ricker)
        """
        
        self.coords = np.array([pos[0], pos[1], pos[2]])
        self.f = wavelet #Ricker
            
        self.setDomainFlag(meshBox) #Flag wether or not the source is in the domain
    
    def __repr__(self):
        return '('+str(self.coords[0])+','+str(self.coords[1])+','+str(self.coords[2])+')'
    
    def getCoord(self):
        return self.coords
    
    def getDomainFlag(self):
        return self.domainFlag
    
    def getFunction(self):
        return self.f
    
    def setDomainFlag(self, meshBox):
        if meshBox[0][0] > self.coords[0] or self.coords[0] > meshBox[0][1] \
            or meshBox[1][0] > self.coords[1] or self.coords[1] > meshBox[1][1] \
            or meshBox[2][0] > self.coords[2] or self.coords[2] > meshBox[2][1]:
                self.domainFlag=0
        else:
            self.domainFlag = 1
            
    def x(self):
        return self.coords[0]
        
    def y(self):
        return self.coords[1]
    
    def z(self):
        return self.coords[2]
            
    
class SourceSet:
    """A class representing a point source
    
    Attributes
    ----------
    source_list :
        List of Source object
    n :
        Number of Source
    """
    
    def __init__(self, source_list=None):
       if source_list is None:
           self.source_list = [] 
           self.n = 0 
       else:
           self.source_list = source_list 
           self.n = len(source_list) #Number of sources
                
                
    def append(self, source):
       self.source_list.append(source)
        

        
