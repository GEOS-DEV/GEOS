'''
Created on 18/02/2021

@author: macpro
'''
import numpy as np
import pygeosx
    

def initialize_pyMesh(problem):
    """Get mesh informations from GEOSX mesh for the seismic acquisition
    
    Parameters
    ----------
    problem :
        A pygeosx.initialize object
    
    Return
    ------
    mesh :
        A Mesh object
    """
    
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
    
    node_position = node_manager.get_wrapper("ReferencePosition").value().to_numpy()
    
    Node_list = [Node(node_position[i], i) for i in range(number_of_node)]
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
    """A class representing a node of the mesh
    
    Attributes
    ----------
    coords :
        Node coordinate
    
    n :
        The global number of the node
    """
        
    def __init__(self, nodecoord, nodenumber):
        """Constructor of Node
        
        Parameters
        ----------
        nodecoord :
            Coordinates of the node
        
        nodenumber :
            Global number of the node in the mesh
        """
                                
        self.coords = np.array([nodecoord[0], nodecoord[1], nodecoord[2])         
        self.n = nodenumber                            
        
    def getCoords(self):
        return self.coords
    
    def x(self):
    	return self.coords[0]
    
    def y(self):
    	return self.coords[1]
    
    def z(self):
    	return self.coords[2]
    
    def getNodeNumber(self):
    	return self.nodenumber
        

class Element:
    """A class representing an element of the mesh
    
    Attributes
    ----------
    node_list :
        List of nodes in that element
        
    speed :
        Wave speed in the element
    
    center :
        Element center coordinates
    
    volume :
        Volume of the element
    
    n :
        The global number of the element
    """
    def __init__(self, node_list, elemspeed, coordcenter, elemvolume, elemnumber):
        """Constructor of Element
        
        Parameters
        ----------
        node_list :
            List of nodes in that element
        
        elemspeed :
            Wave speed in the element
    
        coordcenter :
            Element center coordinates
    
        elemvolume :
            Volume of the element
    
        elemnumber :
            The global number of the element
        """
        
        self.node_list = node_list   
        self.speed  = elemspeed      
        self.center = np.array(coordcenter)   
        self.volume = elemvolume     
        self.n = elemnumber         
    
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
    """A class representing the mesh
    
    Attributes
    ----------
    elem_list :
        List of elements
        
    node_list :
        List of nodes
        
    nbelem :
        Number of elements
    
    nbnode :
        Number of nodes
    
    b_node_list :
        List of boundary nodes
    """
    def __init__(self, elem_list, node_list, boundary_nodes, order):
        """Constructor of Mesh
        
        Parameters
        ----------
        elem_list :
            List of Element objects
        
        node_list :
            List of Node objects
        
        boundary_nodes :
            Number of the Node objects on the boundary of the domain
        
        order :
            Space discretization order
        """
        
        self.elem_list = elem_list                      
        self.node_list = node_list                      
        self.nbelem = len(elem_list)                    
        self.nbnode = len(node_list)                    
        self.ord = order                                
        self.b_node_list = boundary_nodes               
        
    
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
        
        xmin = self.b_node_list[0].coords[0]
        xmax = xmin
        ymin = self.b_node_list[0].coords[1]
        ymax = ymin
        zmin = self.b_node_list[0].coords[2]
        zmax = zmin
    	
        for b_node in self.b_node_list:
            if b_node.x < xmin:
                xmin = b_node.coords[0]
            elif b_node.x > xmin:
                xmax = b_node.coords[0]
            if b_node.y < ymin:
                ymin = b_node.coords[1]
            elif b_node.y > ymin:
                ymax = b_node.coords[1]
            if b_node.z < zmin:
                zmin = b_node.coords[2]
            elif b_node.z > zmin:
                zmax = b_node.coords[2]
        
        box = [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
        
        return box
    	    
    	
    def getMeshSection(self, sourceCoords):
        #Hard coded value, isn't used for now
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
    
    
        
        
