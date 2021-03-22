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
        
        self.x = pos[0]
        self.y = pos[1]
        self.z = pos[2]
        self.coords = np.array([self.x, self.y, self.z])
        self.f = wavelet #Ricker
            
        self.setDomainFlag(meshBox) #Flag wether or not the source is in the domain
    
    def __repr__(self):
        return '('+str(self.x)+','+str(self.y)+','+str(self.z)+')'
    
    def getCoord(self):
        return self.coords
    
    def getDomainFlag(self):
        return self.domainFlag
    
    def getFunction(self):
        return self.f
    
    def setDomainFlag(self, meshBox):
        if meshBox[0][0]>self.x or self.x>meshBox[0][1] \
            or meshBox[1][0]>self.y or self.y>meshBox[1][1] \
            or meshBox[2][0]>self.z or self.z>meshBox[2][1]:
                self.domainFlag=0
        else:
            self.domainFlag = 1
            
    
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
        

        
