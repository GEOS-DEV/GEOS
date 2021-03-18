'''
Created on 18/02/2021

@author: macpro
'''
import numpy as np
import pygeosx
import sys
from mpi4py import MPI
import numpy as np
from functions import *
    

def initialize_pyMesh(problem):

    node_manager = problem.get_group("domain/MeshBodies/mesh/Level0/nodeManager")

    number_of_node = len(node_manager.get_wrapper("ReferencePosition").value().to_numpy())
    
    boundary_node_list = []    
    for i in range(number_of_node):
    	if node_manager.get_wrapper("domainBoundaryIndicator").value().to_numpy()[i]==1:
    	     boundary_node_list.append(i)
    	     

    elem_manager = problem.get_group("domain/MeshBodies/mesh/Level0/ElementRegions/elementRegionsGroup/Region/elementSubRegions/cb")
  
    elem_speed  = problem.get_wrapper("FieldSpecifications/cellVelocity/scale").value()[0]
    
    number_of_elem = len(elem_manager.get_wrapper("nodeList").value().to_numpy())
    
    nb_node_per_elem = elem_manager.get_wrapper("numNodesPerElement").value()[0]
    if nb_node_per_elem==8:
    	disc_ord=1
    elif nb_node_per_elem==64:
    	disc_ord=3
    elif nb_node_per_elem==216:
    	disc_ord=5
    
    
    Node_list = [Node(node_manager.get_wrapper("ReferencePosition").value().to_numpy()[i], i) for i in range(number_of_node)]
    node_array=np.array(Node_list)
    Elem_list = [Element(node_array[elem_manager.get_wrapper("nodeList").value().to_numpy()[i]].tolist(), 
                         elem_speed, 
                         elem_manager.get_wrapper("elementCenter").value().to_numpy()[i], 
                         elem_manager.get_wrapper("elementVolume").value().to_numpy()[i], 
                         i) 
                         for i in range(number_of_elem)]
   
    mesh = Mesh(Elem_list, 
                Node_list, 
                node_array[boundary_node_list].tolist(), 
                disc_ord)
    
    return mesh
    
    
    
    
    
class Node:
    def __init__(self, nodecoord, nodenumber):
        self.x = nodecoord[0]                          #OK
        self.y = nodecoord[1]                          #OK
        self.z = nodecoord[2]                          #OK
        self.coords = np.array([self.x, self.y, self.z])         #OK
        self.n = nodenumber                            #OK
        
    def getCoords(self):
        return self.coords
    
    def x(self):
    	return self.x
    
    def y(self):
    	return self.y
    
    def z(self):
    	return self.z
    
    def getNodeNumber(self):
    	return self.nodenumber
        

class Element:
    def __init__(self, node_list, elemspeed, coordcenter, elemvolume, elemnumber):
        self.node_list = node_list   #OK
        self.speed  = elemspeed      #Acoustic solver
        self.center = np.array(coordcenter)    #OK
        self.volume = elemvolume     #OK
        self.n = elemnumber          #OK
    
    def getNode_List(self):
        return self.node_list
        
    def getCenter(self):
        return self.center
    
    def getSpeed(self):
        return self.speed
    
    def getVolume(self):
        return self.volume
    
    def getNbNode(self):
        return self.nbnode
    
    
        
class Mesh:
    def __init__(self, elem_list, node_list, boundary_nodes, order):
        self.elem_list = elem_list                      #OK
        self.node_list = node_list                      #OK
        self.nbelem = len(elem_list)                    #OK
        self.nbnode = len(node_list)                    #OK
        self.ord = order                                #OK
        self.b_node_list = boundary_nodes               #OK
        
    
    def __repr__(self):
    	return 'Number of elements : ' + str(self.nbelem) + '\n' + 'Number of nodes : ' + str(self.nbnode) + '\n'
        
    def getElem_List(self):
        return self.elem_list    
    
    def getNode_List(self):
        return self.node_list
    
    def getNbElem(self):
        return self.nbelem
    
    def getNbNode(self):
        return self.nbnode
        
    def getOrd(self):
    	return self.ord
    	
    def getBoundaryNode(self):
    	return self.b_node_list
    
    def getMinMaxBoundary(self):
    	xmin = self.b_node_list[0].x
    	xmax = xmin
    	ymin = self.b_node_list[0].y
    	ymax = ymin
    	zmin = self.b_node_list[0].z
    	zmax = zmin
    	
    	for b_node in self.b_node_list:
    	    if b_node.x < xmin:
    	    	xmin = b_node.x
    	    elif b_node.x > xmin:
    	    	xmax = b_node.x
    	    
    	    if b_node.y < ymin:
    	    	ymin = b_node.y
    	    elif b_node.y > ymin:
    	    	ymax = b_node.y
    	    	
    	    if b_node.z < zmin:
    	    	zmin = b_node.z
    	    elif b_node.z > zmin:
    	    	zmax = b_node.z
    	    
    	box = [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
    	
    	return box
    	    
    	
    def getMeshSection(self, sourceCoords):
        maxDistance = 0.05
        
        new_elem_list = []
        node_ind = []
        for elem in self.getElem_List():
            if abs(sourceCoords[0] - elem.getCenter()[0]) <= maxDistance and \
            abs(sourceCoords[1] - elem.getCenter()[1]) <= maxDistance and \
            abs(sourceCoords[2] - elem.getCenter()[2]) <= maxDistance:
               new_elem_list.append(elem)
        node_ind = [node.n for elem in new_elem_list for node in elem.getNode_List()]
        node_ind = list(set(node_ind))
        node_array = np.array(self.node_list)
        
        meshSection = Mesh(new_elem_list, node_array[node_ind].tolist(), self.b_node_list, self.ord)
        
        return meshSection
    
    
        
        
