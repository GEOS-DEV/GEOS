# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 09:00:00 2017
@author: homel1
Geometry object functions for the particle file writer.
"""
import numpy as np                      # math stuff
from sklearn.neighbors import KDTree    # neighbor search with KDTree
import math
from scipy.optimize import minimize_scalar
from scipy.spatial import Voronoi, ConvexHull
np.random.seed(1)

# ===========================================
# DEFAULTS
# ===========================================
defaultSurfaceFlag = 2

# ===========================================
# UTILITY FUNCTIONS
# ===========================================

log_file = "job_file_log"
def log2file(msg):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")

    msg = current_time + " Rank " + str(g_rank) + ": " + msg + "\n"
    mode = MPI.MODE_WRONLY | MPI.MODE_CREATE | MPI.MODE_APPEND

    file_handle = MPI.File.Open(MPI.COMM_SELF, log_file + "_" + str(g_rank), mode)
    file_handle.Set_atomicity(True)

    b = bytearray()
    b.extend(map(ord, msg))

    file_handle.Write(b)
    file_handle.Sync()
    file_handle.Close()


def print2file(file_name, text):
    with open(file_name, 'a') as f:
        f.write(text + "\n")


def countFileLines(fname):
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b: break
            yield b

    with open(fname, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count

def fOffsetLineNum(fname, line_number):
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b: break
            yield b

    if line_number == 0:
        return 0

    count = 0
    with open(fname, "rb") as f:
        for buf in _make_gen(f.raw.read):
            count += buf.count(b"\n")
            if count >= line_number:
                chunk_pos = f.tell()
                break

        # Back track to line_number
        for ci in range(len(buf)-1, -1, -1):
            c = buf[ci]
            # print("Character ",c, )
            chunk_pos -= 1
            if c == 10: #b"\n":
                if count == line_number:
                    # Add one char because first char of new line immediately 
                    #     follows new line char
                    return chunk_pos + 1 
                count -= 1


#Might be a faster way to do this
def fileOffsetFromLineNum(fname, line_number):
    char_pos = 0
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            if i == line_number:
                offset = len(line)
                return char_pos
            else:
                char_pos += len(line)


# used for periodicity
def mapToRange(x, xmin, xmax):
    # Make sure all inputs are np arrays
    xmin = np.array(xmin)
    xmax = np.array(xmax)
    dx = xmax-xmin
    x = np.array(x) - xmin

    x = x * (x < dx) + (x >= dx) * (np.fmod(x, dx))
    x = x * (x >= 0) + (x < 0) * (dx - np.fmod(np.absolute(x), dx))

    return x + xmin


# Smooth sigmoid interpolation between rf and rf0
def smoothStep3(r, rf, rf0):
  if r >= rf0:
    return 0

  if r <= rf:
    return 1

  r = (rf0-r)/(rf0-rf)

  return r*r*(3-2*r)


# Checks if point x lies outside of box with minimum corner at x0 and 
#     dimensions dx
def outside(x, x0, dx, periodic):
    x = periodic * (np.fmod(x - x0, dx) + x0) + np.logical_not(periodic) * x
    return np.any(np.logical_or(x < x0, x > x0 + dx)), x
    

# Grid object used for efficiently performing random sequential absorption and 
#     poisson disc sampling
class Grid:
    def __init__(self, cell_size, **kwargs):
        self.cell_size = cell_size
        self.dim = kwargs.get('dim', 3)
        self.x0 = np.array(kwargs.get('x0', np.zeros(self.dim)))
        self.dx = np.array(kwargs.get('dx', np.ones(self.dim)))
        self.periodic = kwargs.get('periodic', [True for k in range(self.dim)])
        self.num_cells = np.ceil(self.dx / self.cell_size).astype(int)
        self.has_fractional_cells = np.multiply(self.num_cells, self.dx) > self.dx
        # Cell value of -1 indicates empty cell
        self.cells = -np.ones(tuple(self.num_cells)).astype(int)  

    def outside(self, x):
        return outside(x, self.x0, self.dx, self.periodic)

    def pos2cell(self, x):
        return np.floor((x - self.x0) / self.cell_size).astype(int)

    def place(self, cell, x):
        self.cells[tuple(cell)] = x
        return self

    def at(self, cell):
        return self.cells[tuple(cell)]

    def occupied(self, cell):
        return self.cells[tuple(cell)] >= 0

    def cells_in_range(self, cell, window):
        search_min = cell - window
        search_max = cell + window

        # If direction doesn't have periodic boundaries then cap search window 
        #     to min and max cells
        # If grid has fractional cells and candidate cell is within range of 
        #     them add 1 to search direction
        search_min = ( (search_min - np.ones(self.dim).astype(int) 
            * np.logical_and(self.has_fractional_cells, cell < 1)) 
            * self.periodic + np.logical_not(self.periodic) * np.maximum(
            np.zeros(self.dim).astype(int), search_min) )

        search_max = (search_max + np.ones(self.dim).astype(int) * np.logical_and(self.has_fractional_cells, cell >= self.num_cells - 3)) * self.periodic + np.logical_not(self.periodic) * np.minimum(
            self.num_cells-1, search_max)

        search_range = search_max-search_min
        larger_than_cells = search_range > ( self.num_cells - 1 )

        search_min = np.zeros(self.dim).astype(int) * larger_than_cells + np.logical_not(larger_than_cells) * search_min
        search_max = (self.num_cells - 1) * larger_than_cells + np.logical_not(larger_than_cells) * search_max

        # Create list of all cells indices in search area
        cells_to_search_ranges = [np.array([]) for k in range(self.dim)]
        for k in range(self.dim):
            cells_to_search_ranges[k] = np.array(list(range(search_min[k], search_max[k] + 1))).astype(int)
        cells_to_search = combvec(cells_to_search_ranges).astype(int)


        cell_vals = []
        for c in cells_to_search:
            # Check if c index needs to be wrapped in periodic direction
            c = wrap(c, np.zeros(self.dim), self.num_cells - 1)

            val = self.cells[tuple(c)]
            if val >= 0:
                cell_vals.append(val)

        return cell_vals


# wraps integer index i to range of min and max (both inclusive)
# accepts nd arrays
def wrap(i, lower, upper):
    # Make sure all inputs are np arrays
    i = np.array(i).astype(int)
    lower = np.array(lower).astype(int)
    upper = np.array(upper).astype(int)

    r = upper - lower + 1
    i = i - lower

    i = i * (i < r) + (i >= r) * (lower + np.mod(i, r))
    i = i * (i >= 0) + (i < 0) * (upper - np.mod(np.absolute(i + 1), r))

    return i.astype(int)


# Generates all possible combinations of elements from list of offsets for each direction
def combvec(offsets):
    dim = len(offsets)

    combs = np.reshape(offsets[0], (len(offsets[0]), 1))
    for i in range(1, dim):
        new_combs = np.empty((0, i+1))
        for j in offsets[i]:
            for k in combs:
                new_combs = np.append(new_combs, np.reshape(np.append(k, j), (1, len(k)+1)), axis=0)

        combs = new_combs

    return combs


# Randomly samples points on surface of unit n-sphere
def random_direction(dim=3):
    direction = np.random.normal(size=dim)
    return direction / math.sqrt(np.sum(direction**2))


# Adds images of pores for nearest neighbors
def add_images(pores, **kwargs):
    # Setup parameters
    dim = kwargs.get('dim', 3)
    x0 = np.array(kwargs.get('x0', np.zeros(dim)))
    dx = np.array(kwargs.get('dx', np.ones(dim)))
    periodic = kwargs.get('periodic', [True for k in range(dim)])

    num_pores = len(pores)
    boundary_pores = np.empty((0, dim + 1))
    offsets = [np.array([0, -1, 1]) for i in range(dim)]
    offset_combs = combvec(offsets)
    for ff in offset_combs:
        for p in pores:
            boundary_pores = np.append(boundary_pores, np.reshape(p + np.append(np.multiply(ff, dx), [0]), (1, dim + 1)), axis=0)

    return boundary_pores


# Adds pores images of those that lie on the boundary of the RVE
# x0 = corner with minimum x, y, and z coordinates
# dx = width, height, and length of domain
# pores = list of pores, format = x, y, z, r
def add_pore_images(pores, **kwargs):
    # Setup parameters
    dim = kwargs.get('dim', 3)
    x0 = np.array(kwargs.get('x0', np.zeros(dim)))
    dx = np.array(kwargs.get('dx', np.ones(dim)))
    periodic = kwargs.get('periodic', [True for k in range(dim)])

    num_pores = len(pores)

    indices = []
    boundary_pores = np.empty((0, dim + 1))
    for i in range(num_pores):
        offsets = [np.array([0]) for k in range(dim)]

        for k in range(dim):
            if periodic[k]:
                # Pore on - boundary
                if pores[i, k] - pores[i, dim] < x0[k]:
                    offsets[k] = np.append(offsets[k], 1.0)

                # Pore on + boundary
                if pores[i, k] + pores[i, dim] > x0[k] + dx[k]:
                    offsets[k] = np.append(offsets[k], -1.0)

        offset_combs = combvec(offsets)
        for ff in offset_combs:
            indices.append(i)
            boundary_pores = np.append(boundary_pores, np.reshape(pores[i, :] + np.append(np.multiply(ff, dx), [0]), (1, dim+1)), axis=0)

    return boundary_pores, indices


# Imports pores from a LAMMPS dump file
# pack_file = location of pack.dump file
# offset = translation of coordinates for simulation (default [0.0, 0.0, 0.0])
# add_boundary_images = boolean whether to include pores that lie on boundaries on other periodic surface 
# (e.g for generating the foam meshes in the particleFileWriter)
def import_pores(pack_file, **kwargs):
    # Setup parameters
    dim = kwargs.get('dim', 3)

    # x0 = np.array(kwargs.get('x0', np.zeros(dim)))
    # dx = np.array(kwargs.get('dx', np.ones(dim)))
    offset = np.array(kwargs.get('offset', np.zeros(dim)))
    scale = kwargs.get('scale', 1.0)
    add_boundary_images = kwargs.get('add_boundary_images', True)
    radius_scale = kwargs.get('radius_scale', 1.0)

    print("Reading pores from " + pack_file)

    file = open(pack_file, "r")
    lines = file.readlines()
    file.close()

    # Stores current data block being processed and list of voids
    curr_block = ""
    pores = []
    box = []
    curr_iter = 0
    num_atoms = 0
    count = 0
    for line in lines:
        count += 1
        # print("Line{}: {}".format(count, line.strip()))
        if line == "" or line == "\n":
            continue

        # Split string by spaces
        line_terms = line.split()

        # Determined what to do based on first entry in line (ITEM: = some header info, number
        if line_terms[0] == "ITEM:":
            curr_block = line_terms[1]

            if curr_block == "BOX":
                # For box get periodicity from same line
                periodic = [line_terms[3] == 'pp', line_terms[4] == 'pp', line_terms[5] == 'pp']
        elif curr_block == "ITERATION":
            curr_iter = int(line)
        elif curr_block == "NUMBER":
            num_atoms = int(line)
        elif curr_block == "BOX":
            box.append([float(i) for i in line_terms])
        elif curr_block == "ATOMS":
            atom = [float(x) for x in line_terms]
            # Output format radius, x, y, z
            pores.append([scale*(atom[2] + offset[0]), scale*(atom[3] + offset[1]), scale*(atom[4] + offset[2]), scale * radius_scale * atom[5]])

    pores = np.array(pores)  # convert python array to numpy array

    x0 = [scale*(box[0][0] + offset[0]), scale*(box[1][0] + offset[1]), scale*(box[2][0] + offset[2])]
    dx = [scale*box[0][1], scale*box[1][1], scale*box[2][1]]

    print("Box after offset and scaling:")
    print("\t" + str(x0[0]) + " " + str(dx[0]))
    print("\t" + str(x0[1]) + " " + str(dx[1]))
    print("\t" + str(x0[2]) + " " + str(dx[2]))
    print("Periodicity: " + str(periodic[0]) + " " + str(periodic[1]) + " " + str(periodic[2]))
    print("Pores from file: " + str(len(pores)))
    if add_boundary_images:
        output_pores = add_pore_images(pores, x0=x0, dx=dx, dim=dim, periodic=periodic)
        print("Pore with images: " + str(len(output_pores)))
    else:
        output_pores = pores

    # Reformt pores so radius is first, then x, y, z
    radii = output_pores[:,dim].copy()
    output_pores[:,1:dim+1] = output_pores[:,0:dim] # shift x, y, z right one column
    output_pores[:,0] = radii

    return np.array(output_pores)


# Generates densely packed random points in box domain with minimum separation
# Generalized for 2 and 3 dimensions
def poisson(spacing, **kwargs):
    # Setup parameters
    seed = kwargs.get('seed', np.random.randint(0, 1e8))
    dim = kwargs.get('dim', 3)
    x0 = np.array(kwargs.get('x0', np.zeros(dim)))
    dx = np.array(kwargs.get('dx', np.ones(dim)))
    trials = kwargs.get('trials', 30)
    periodic = kwargs.get('periodic', [True for k in range(dim)])

    # Grid cell size if determined given longest cell size diagonal is radius of spacing
    grid = Grid(spacing / math.sqrt(dim), dim=dim, x0=x0, dx=dx, periodic=periodic)

    # First seed point is at center of domain
    seed_pts = np.array([dx/2 + x0])  # Check that this makes the right size array
    pores = np.empty((0, dim+1))

    # Set seed for random number generation
    np.random.seed(seed)

    # Terminates when there are no more seed points left
    while seed_pts.shape[0] != 0:
        # Get the first point in the seed array
        curr_seed_pt = seed_pts[0, :]
        seed_pts = np.delete(seed_pts, 0, axis=0)

        num_rejected = 0
        while num_rejected < trials:
            # Pick random distance and random angle from seed point
            d = spacing*(1 + np.random.random() * 0.1)  # random spacing between seed and candidate point
            direction = random_direction(dim)
            candidate_pt = curr_seed_pt + d * direction

            accept = True  # Accept candidate point unless found otherwise

            # If periodic boundaries are not enabled test point must be inside domain
            is_outside, candidate_pt = grid.outside(candidate_pt)
            if is_outside:
                num_rejected = num_rejected + 1
                continue

            candidate_cell = grid.pos2cell(candidate_pt)
            if grid.occupied(candidate_cell):
                accept = False
            else:
                neighbor_pores = grid.cells_in_range(candidate_cell, 2*np.ones(dim).astype(int))

                # print( "Neighbors", neighbor_pores)
                # for n in neighbor_pores:
                #     print(n, pores[n,:-1], grid.pos2cell(pores[n,:-1]))
                # print("Done Neighbors")

                # Check every cell for occupant, if any are within spacing reject candidate point
                for n in neighbor_pores:

                    # Check distance between candidate point and point already there
                    ddx = candidate_pt - pores[n, :-1]

                    # Use shortest image distance if that point lies across periodic boundary
                    ddxx = np.logical_not(periodic) * ddx + periodic * np.minimum(np.absolute(ddx + dx),
                                                                                  np.minimum(np.absolute(ddx),
                                                                                             np.absolute(ddx - dx)))
                    # print(candidate_pt, pores[n, :-1], ddxx, np.sqrt(np.sum(ddxx**2)), np.sqrt(np.sum(ddxx**2)) < spacing)

                    if np.sum(ddxx**2) < spacing**2:
                        accept = False
                        break  # Once a conflicting point is found no need to search the remaining cells

            if accept:
                # Add pore to list of pores and seed points plus mark grid cell as occupied
                grid = grid.place(tuple(candidate_cell), pores.shape[0])
                new_pore = np.append(candidate_pt, spacing/2)
                pores = np.append(pores, np.reshape(new_pore, (1, dim+1)), axis=0)
                seed_pts = np.append(seed_pts, np.reshape(candidate_pt, (1, dim)), axis=0)
            else:
                num_rejected = num_rejected + 1

    return pores


# Creates pores which some minimum spacing via random sequential adsorption (RSA)
# N = number of desired particles
# spacing = minimum separation between spheres
# r = function alias for sampling radii (To be implemented, right now only monodisperse)
# seed = random seed for RSA
# x0 = corner with minimum x, y, and z coordinates
# dx = width, height, and length of domain
# trials = number of rejected trials before terminating RSA
def rsa(N, spacing, **kwargs):
    # Setup parameters
    seed = kwargs.get('seed', np.random.randint(0, 1e8))
    dim = kwargs.get('dim', 3)
    x0 = np.array(kwargs.get('x0', np.zeros(dim)))
    dx = np.array(kwargs.get('dx', np.ones(dim)))
    trials = kwargs.get('trials', 30)
    periodic = kwargs.get('periodic', [True for k in range(dim)])

    # Grid cell size if determined given longest cell size diagonal is radius of spacing
    grid = Grid(spacing / math.sqrt(dim), dim=dim)

    pores = np.empty((0, dim + 1))

    num_rejected = 0
    np.random.seed(seed)
    while pores.shape[0] < N and num_rejected <= trials:
        x = np.empty(dim + 1)
        x[dim] = spacing/2
        for k in range(dim):
            if periodic[k]:
                x[k] = np.random.uniform(0.0, dx[k])
            else:
                x[k] = np.random.uniform(0.0+spacing/2, dx[k]-spacing/2)

        candidate_cell = grid.pos2cell(x[0:-1])

        accept = True
        if grid.occupied(candidate_cell):
            accept = False
        else:
            neighbor_pores = grid.cells_in_range(candidate_cell, 2*np.ones(dim).astype(int))
            for p in neighbor_pores:
                d = x[0:dim] - pores[p, 0:-1]
                dd = d * np.logical_not(periodic) + periodic * np.minimum(np.abs(d), np.abs(d + dx))
                dd = d * np.logical_not(periodic) + periodic * np.minimum(np.abs(dd), np.abs(d - dx))

                if np.sum(np.square(dd)) < spacing ** 2:
                    accept = False
                    break

        if accept:
            grid = grid.place(candidate_cell, pores.shape[0])
            pores = np.append(pores, np.reshape(x, (1, dim + 1)), axis=0)
            num_rejected = 0
        else:
            num_rejected = num_rejected + 1

    for i in range(len(pores)):
        pores[i, 0:dim] = pores[i, 0:dim] + x0

    return pores



# ===========================================
# GEOMETRY OBJECTS
# ===========================================

#############################################
class union:
  'Geometry object for creating a union of two objects'
  def __init__(self,A,B):
    self.A = A
    self.B = B
    self.type = 'union'

    # copy properties from object A.
    self.v = A.v
    self.mat = A.mat
    self.group = A.group
    self.damage = A.damage

  def isInterior(self,x):
    if( self.A.isInterior(x)or self.B.isInterior(x) ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # add this.
    print('isSurface not supported for union!')
    return 0

#############################################
class intersection:
  'Geometry object for creating an intersection of two objects'
  # material velocity and damage are inherited from A.
  def __init__(self,A,B):
    self.A = A
    self.B = B
    self.type = 'intersection'
    
    # copy properties from object A.
    self.v = A.v
    self.mat = A.mat
    self.group = A.group
    self.damage = A.damage


  def getMat(self,x):
    return self.A.getMat(x)

  def isInterior(self,x):
    if( self.A.isInterior(x) and self.B.isInterior(x) ):
      return True
    else:
      return False

  def isSurface(self,x,l):
    # add this.
    if( self.A.isSurface(x,l) or self.B.isSurface(x,l) ):
      return defaultSurfaceFlag
    else:
      return 0

  # bounding box here is just for sorting:
  def xMin(self):
    # This might be used for binning objects
    return self.A.xMin()

  def xMax(self):
    # This might be used for binning objects
    return self.A.xMax()

#############################################
class difference:
  'Geometry object for creating a difference of two objects'
  def __init__(self,A,B):
    self.A = A
    self.B = B
    self.type = 'difference'

    # copy properties from object A.
    self.v = A.v
    self.mat = A.mat
    self.group = A.group
    self.damage = A.damage

  def isInterior(self,x):
    if( self.A.isInterior(x) and not self.B.isInterior(x) ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # add this.
    print('isSurface not supported for difference!')
    return 0

    # bounding box here is just for sorting:
  def xMin(self):
    # This might be used for binning objects
    return self.A.xMin()

  def xMax(self):
    # This might be used for binning objects
    return self.A.xMax()

#############################################
class reflect:
  'Geometry object for creating a rotated object around point x0 and axis a0 by angle alpha'
  def __init__(self,name,X,x0,a0):
    self.name = name
    self.X = X         # object to be rotated
    self.x0 = np.array(x0)           # center of reflection
    self.a0 = a0/np.linalg.norm(a0)  # unit vector defining reflection direction.
    self.type = 'reflect'

    # copy properties from child.
    self.v = X.v
    self.mat = X.mat

    self.group = X.group
    self.damage = X.damage

  def getMat(self,pt):
    aa = np.array([
      [self.a0[0]*self.a0[0], self.a0[0]*self.a0[1], self.a0[0]*self.a0[2]], 
      [self.a0[1]*self.a0[0], self.a0[1]*self.a0[1], self.a0[1]*self.a0[2]], 
      [self.a0[2]*self.a0[0], self.a0[2]*self.a0[1], self.a0[2]*self.a0[2]]
      ])
    x = pt - 2.0*( pt - self.x0 ).dot(aa)
    return self.X.getMat(x)

  def isInterior(self,pt):
    aa = np.array([
      [self.a0[0]*self.a0[0], self.a0[0]*self.a0[1], self.a0[0]*self.a0[2]], 
      [self.a0[1]*self.a0[0], self.a0[1]*self.a0[1], self.a0[1]*self.a0[2]], 
      [self.a0[2]*self.a0[0], self.a0[2]*self.a0[1], self.a0[2]*self.a0[2]]
      ])
    x = pt - 2.0*( pt - self.x0 ).dot(aa)
    if( self.X.isInterior(x) ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    aa = np.array([
      [self.a0[0]*self.a0[0], self.a0[0]*self.a0[1], self.a0[0]*self.a0[2]], 
      [self.a0[1]*self.a0[0], self.a0[1]*self.a0[1], self.a0[1]*self.a0[2]], 
      [self.a0[2]*self.a0[0], self.a0[2]*self.a0[1], self.a0[2]*self.a0[2]]
      ])
    x = pt - 2.0*( pt - self.x0 ).dot(aa)
    if( self.X.isSurface(x,l) ):
      return defaultSurfaceFlag
    else:
      return 0
    # bounding box here is just for sorting:

  def xMin(self):
    # This might be used for binning objects
    return np.array([-1.0e9,-1.0e9,-1.0e9])

  def xMax(self):
    # This might be used for binning objects
    return np.array([1.0e9,1.0e9,1.0e9])

#############################################
class rotate:
  'Geometry object for creating a rotated object around point x0 and axis a0 by angle alpha'
  def __init__(self,name,X,x0,a0,alpha):
    self.name = name
    self.X = X         # object to be rotated
    self.x0 = np.array(x0)       # center of rotation
    self.a0 = a0/np.linalg.norm(a0)       # vector defining rotation axis
    self.alpha = alpha # right handed rotation (radians) about axis
    self.type = 'rotate'

    # copy properties from child.
    self.v = X.v

    self.mat = X.mat

    self.group = X.group
    self.damage = X.damage

  def getMat(self,pt):
    aa = np.array([
      [self.a0[0]*self.a0[0], self.a0[0]*self.a0[1], self.a0[0]*self.a0[2]], 
      [self.a0[1]*self.a0[0], self.a0[1]*self.a0[1], self.a0[1]*self.a0[2]], 
      [self.a0[2]*self.a0[0], self.a0[2]*self.a0[1], self.a0[2]*self.a0[2]]
      ])
    A = np.array([
      [ 0         , -self.a0[2], self.a0[1]], 
      [ self.a0[2],           0,-self.a0[0]], 
      [-self.a0[1],  self.a0[0],          0]
      ])
    R = np.cos(self.alpha)*np.identity(3) + (1-np.cos(self.alpha))*aa + np.sin(self.alpha)*A
    x = self.x0 + (R.transpose()).dot(np.array(pt)-self.x0)
    return self.X.getMat(x)

  def isInterior(self,pt):
    aa = np.array([
      [self.a0[0]*self.a0[0], self.a0[0]*self.a0[1], self.a0[0]*self.a0[2]], 
      [self.a0[1]*self.a0[0], self.a0[1]*self.a0[1], self.a0[1]*self.a0[2]], 
      [self.a0[2]*self.a0[0], self.a0[2]*self.a0[1], self.a0[2]*self.a0[2]]
      ])
    A = np.array([
      [ 0         , -self.a0[2], self.a0[1]], 
      [ self.a0[2],           0,-self.a0[0]], 
      [-self.a0[1],  self.a0[0],          0]
      ])
    R = np.cos(self.alpha)*np.identity(3) + (1-np.cos(self.alpha))*aa + np.sin(self.alpha)*A
    x = self.x0 + (R.transpose()).dot(np.array(pt)-self.x0)
    if( self.X.isInterior(x) ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    aa = np.array([
      [self.a0[0]*self.a0[0], self.a0[0]*self.a0[1], self.a0[0]*self.a0[2]], 
      [self.a0[1]*self.a0[0], self.a0[1]*self.a0[1], self.a0[1]*self.a0[2]], 
      [self.a0[2]*self.a0[0], self.a0[2]*self.a0[1], self.a0[2]*self.a0[2]]
      ])
    A = np.array([
      [ 0         , -self.a0[2], self.a0[1]], 
      [ self.a0[2],           0,-self.a0[0]], 
      [-self.a0[1],  self.a0[0],          0]
      ])
    R = np.cos(self.alpha)*np.identity(3) + (1-np.cos(self.alpha))*aa + np.sin(self.alpha)*A
    x = self.x0 + (R.transpose()).dot(np.array(pt)-self.x0)
    if( self.X.isSurface(x,l) ):
      return defaultSurfaceFlag
    else:
      return 0

  def xMin(self):
    # This might be used for binning objects
    return np.array([-1.0e9,-1.0e9,-1.0e9])

  def xMax(self):
    # This might be used for binning objects
    return np.array([1.0e9,1.0e9,1.0e9])

#############################################
class shell:
  'Geometry object for creating a spherical shell defined by center and inner and outer radii'
  def __init__(self,name,x0,ri,ro,v,mat,group,damage):
    self.name = name
    self.type = 'sphere'
    self.x0 = np.array(x0)
    self.ri = ri
    self.ro = ro
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage

  def isInterior(self,pt):
    # is the point within the object
    x = np.array(pt)
    if ( self.ri <= np.linalg.norm(self.x0-x) < self.ro ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    x = np.array(pt)
    if ( ( np.linalg.norm(self.x0-x) > (self.ro - l) ) or ( np.linalg.norm(self.x0-x) < (self.ri + l) ) ):
      return defaultSurfaceFlag
    else:
      return 0

#############################################
class sphere:
  'Geometry object for creating a sphere defined by center and radius'
  def __init__(self,name,x0,r,v,mat,group,damage):
    self.name = name
    self.type = 'sphere'
    self.x0 = np.array(x0)
    self.r = r
    self.r2 = r*r
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage

  def isInterior(self,pt):
    # is the point within the object
    x = np.array(pt)
    # if ( np.linalg.norm(self.x0-x) < self.r ):
    if ( (self.x0[0]-x[0])*(self.x0[0]-x[0]) + (self.x0[1]-x[1])*(self.x0[1]-x[1]) + (self.x0[2]-x[2])*(self.x0[2]-x[2]) < self.r2 ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    x = np.array(pt)
    #if ( np.linalg.norm(self.x0-x) > (self.r - l) ):
    if ( (self.x0[0]-x[0])*(self.x0[0]-x[0]) + (self.x0[1]-x[1])*(self.x0[1]-x[1]) + (self.x0[2]-x[2])*(self.x0[2]-x[2]) > (self.r - l)*(self.r - l) ):
      return defaultSurfaceFlag
    else:
      return 0

  def xMin(self):
    # This might be used for binning objects
    return self.x0[0] - self.r

  def xMax(self):
    # This might be used for binning objects
    return self.x0[0] + self.r

#############################################
class sphere0:
  'Geometry object for creating a sphere defined by center and radius, but where surface is always false (useful for inclusions)'
  def __init__(self,name,x0,r,v,mat,group,damage):
    self.name = name
    self.type = 'sphere'
    self.x0 = np.array(x0)
    self.r = r
    self.r2 = r*r
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage

  def isInterior(self,pt):
    # is the point within the object
    x = np.array(pt)
    # if ( np.linalg.norm(self.x0-x) < self.r ):
    if ( (self.x0[0]-x[0])*(self.x0[0]-x[0]) + (self.x0[1]-x[1])*(self.x0[1]-x[1]) + (self.x0[2]-x[2])*(self.x0[2]-x[2]) < self.r2 ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
      return 0

  def xMin(self):
    # This might be used for binning objects
    return self.x0[0] - self.r

  def xMax(self):
    # This might be used for binning objects
    return self.x0[0] + self.r

#############################################

class directionWrapper:
  def __init__(self,name,subObject,matDir):
    self.subObject = subObject
    self.v = subObject.v
    self.mat = subObject.mat
    self.group = subObject.group
    self.damage = subObject.damage
    self.matDir = matDir

  def isInterior(self,pt):
     return self.subObject.isInterior(pt)

  def isSurface(self,pt,l):
   return self.subObject.isSurface(pt, l)

  def xMin(self):
    return self.subObject.xMin()

  def xMax(self):
    return self.subObject.xMax()

############################################

class surfaceFlagWrapper:
  def __init__(self,name,subObject,flagType):
    self.subObject = subObject
    self.v = subObject.v
    self.mat = subObject.mat
    self.damage = subObject.damage
    self.flagType = flagType

  def isInterior(self,pt):
    return self.subObject.isInterior(pt)

  def isSurface(self,pt,l):
    return self.subObject.isSurface(pt, l)

  def surfaceNormal(self, pt):
    return self.subObject.surfaceNormal(pt)

  def surfacePosition(self,pt):
    return self.subObject.surfacePosition(pt)

  def getGroup(self, pt):
    return self.subObject.getGroup(pt)

  def getMatDir(self, pt):
    return self.subObject.getMatDir(pt)

  def xMin(self):
    return self.subObject.xMin()

  def xMax(self):
    return self.subObject.xMax()

############################################

class strengthScaleWrapper:
  def __init__(self,name,subObject,strengthScale):
    self.subObject = subObject
    self.v = subObject.v
    self.mat = subObject.mat
    self.group = subObject.group
    self.damage = subObject.damage
    self.strengthScale = strengthScale

  def isInterior(self,pt):
     return self.subObject.isInterior(pt)

  def isSurface(self,pt,l):
   return self.subObject.isSurface(pt, l)

  def xMin(self):
    return self.subObject.xMin()

  def xMax(self):
    return self.subObject.xMax()

#############################################

class ellipsoid:
  'Geometry object for creating a grid aligned ellipsoid defined by center and three lengths'
  def __init__(self,name,x0,a,b,c,v,mat,group,damage):
    self.name = name
    self.type = 'ellipsoid'
    self.x0 = np.array(x0)
    self.a = a
    self.b = b
    self.c = c
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage
    #self.translate = np.array([0,0,0])
    #self.transform = np.array([[1,0,0],[0,1,0],[0,0,1]])

  def isInterior(self,pt):
    # is the point within the object
    x = pt[0] - self.x0[0]
    y = pt[1] - self.x0[1]
    z = pt[2] - self.x0[2]
    if ( (x/self.a)*(x/self.a) + (y/self.b)*(y/self.b) + (z/self.c)*(z/self.c)  < 1 ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?    
    x = pt[0] - self.x0[0]
    y = pt[1] - self.x0[1]
    z = pt[2] - self.x0[2]

    a = max( l/100, self.a - l)
    b = max( l/100, self.b - l)
    c = max( l/100, self.c - l)

    if ( (x/a)*(x/a) + (y/b)*(y/b) + (z/c)*(z/c) > 1 ):
      return defaultSurfaceFlag
    else:
      return 0


#############################################

class crystal:
  def __init__(self, name, center, axis, height, rmin, rmax, v, mat, group, damage):
    self.name = name
    self.center = center
    self.axis = axis
    self.height = height
    self.rmin = rmin
    self.rmax = rmax
    self.v = v
    self.mat = mat
    self.group = group
    self.damage = damage

    n0 = np.array(self.axis)
    n0 = n0 / np.linalg.norm(n0)
    self.matDir = n0

    n1 = np.array([3/(5*np.sqrt(2)), 2*np.sqrt(2)/5, 1/np.sqrt(2)]) 
    if np.all(n0 == n1):
      n1 = np.array([3/(5*np.sqrt(2)), 1/np.sqrt(2), 2*np.sqrt(2)/5])

    n1 = np.dot(np.identity(3)-np.tensordot(n0, n0, axes=0), n1)
    n1 = n1 / np.linalg.norm(n1)
    n2 = np.cross(n0, n1)

    # print(n0, n1, n2)
    # print(np.linalg.norm(n0), np.linalg.norm(n1), np.linalg.norm(n2))

    faceNormals = []

    # Face set 1
    for theta in np.linspace(0.0, 5*math.pi/3, 6): #math.pi/3):
      nn = np.cos(theta)*n1 + np.sin(theta)*n2
      nn = nn / np.linalg.norm(nn)
      faceNormals.append(nn)

    # Face set 2
    for theta in np.linspace(math.pi/6, 3*math.pi/2, 3): # 2*math.pi/3):
      nn = n0 + np.cos(theta)*n1 + np.sin(theta)*n2
      nn = nn / np.linalg.norm(nn)
      faceNormals.append(nn)

    # Face set 3
    for theta in np.linspace(math.pi/6, 3*math.pi/2, 3): #2*math.pi/3):
      nn = -n0 + np.cos(theta)*n1 + np.sin(theta)*n2
      nn = nn / np.linalg.norm(nn)
      faceNormals.append(nn)

    faceNormals.append(n0)
    faceNormals.append(-n0)

    self.faceNormals = faceNormals

    self.ci = np.vstack((np.zeros((6,3)), 
                         np.tile(self.height/2*n0,(3,1)),
                         np.tile(-self.height/2*n0,(3,1)),
                         np.zeros((2,3))))
    self.ri = (self.rmax-self.rmin)*np.random.rand(14,1)+self.rmin
    self.ri[-2:] = self.ri[-2:] + self.height

  def isInterior(self, pt):
    pc = pt - self.center
    for i in range(len(self.faceNormals)):
      if np.dot(self.faceNormals[i], pc-self.ci[i,:]) > self.ri[i]:
        return False
    return True

  def isSurface(self, pt, l):
    pc = pt - self.center
    for i in range(len(self.faceNormals)):
      if np.dot(self.faceNormals[i], pc-self.ci[i,:]) > self.ri[i] - l:
        return defaultSurfaceFlag
    return 0


#############################################
class indentor:
  'Geometry object for creating a grid-aligned indentor defined by angle and number of facets'
  def __init__(self,name,x0,nFaces,alpha,v,mat,group,damage):
    self.name = name
    self.type = 'indentor'
    self.x0 = np.array(x0)         # tip coordinate, assumes identor is aligned with z-axis with tip in -z direction
    self.v = v if callable(v) else np.array(v) # initial velocity
    self.nFaces = nFaces # number of facets
    self.alpha = np.radians(alpha)   # facet angle, input in degrees, stored in radians, 
    self.mat = mat
    self.group = group
    self.damage = damage

    # note: alpha=65.3 deg for berkovich with nFaces=3
    # note: alpha=68 deg for vickers with nFaces=4
    #print("indentor, alpha = ",self.alpha,", n=",self.nFaces)
    
    # create set of unit vectors orthogonal to each facet.
    theta = 2.*np.pi/nFaces # angle between normals about indentor axis
    self.n = []
    k=0.
    for i in range(0,nFaces):
      k=k+1.
      nz = -np.sin(self.alpha)
      r = np.cos(self.alpha)
      nx = r*np.cos(theta*(k-0.5))
      ny = r*np.sin(theta*(k-0.5))
      #print([nx,ny,nz])
      (self.n).append(np.array([nx,ny,nz]))

  def isInterior(self,pt):
    # is the point within the object    
    x = np.array(pt-self.x0)
    interior = True
    for i in range(0,self.nFaces):
      if ( x.dot(self.n[i]) > 0.0 ):
        interior = False
    return interior
    
  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    # check only applied to interior points.
    x = np.array(pt-self.x0)
    surface = 0
    for i in range(0,self.nFaces):
      if ( x.dot(self.n[i]) > -1.0*l ):
        surface = defaultSurfaceFlag
    return surface

  def xMin(self):
    # This might be used for binning objects
    return np.array([-1.0e9,-1.0e9,-1.0e9])

  def xMax(self):
    # This might be used for binning objects
    return np.array([1.0e9,1.0e9,1.0e9])

#############################################
class box:
  'Geometry object for creating a grid-aligned box defined by two corners.'
  def __init__(self,name,x0,x1,v,mat,group,damage):
    self.name = name
    self.type = 'box'
    self.x0 = np.array(x0)
    self.x1 = np.array(x1)
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage

  def isInterior(self,pt):
    # is the point within the object    
    x = np.array(pt)
    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])
    if ( (xmin<=x[0] and x[0]<xmax) and (ymin<=x[1] and x[1]<ymax) and (zmin<=x[2] and x[2]<zmax) ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    x = np.array(pt)
    xmin = min(self.x0[0],self.x1[0]) + l
    xmax = max(self.x0[0],self.x1[0]) - l
    ymin = min(self.x0[1],self.x1[1]) + l
    ymax = max(self.x0[1],self.x1[1]) - l
    zmin = min(self.x0[2],self.x1[2]) + l
    zmax = max(self.x0[2],self.x1[2]) - l
    if ( (xmin>=x[0] or x[0]>xmax) or (ymin>=x[1] or x[1]>ymax) or (zmin>=x[2] or x[2]>zmax) ):
      if abs(pt[1] - self.x0[1]) < abs(self.x0[1]-self.x1[1])/2:
        surfaceNormal = [ 0.0, -1.0, 0.0]
      else:
        surfaceNormal = [ 0.0, 1.0, 0.0]
      return defaultSurfaceFlag, surfaceNormal
    else:
      return 0, None 

  def surfacePosition(self, pt):
    x = np.array(pt[0:2])
    x0 = self.x0[0:2]
    x1 = self.x1[0:2]

    dx = np.vstack((x0-x,x1-x))

    dxI = np.argmin(np.absolute(dx),axis=0)
    # dxMin = np.array([dx[dxI[0]][0], dx[dxI[1]][1], dx[dxI[2]][2]])
    dxMin = np.array([dx[dxI[0]][0], dx[dxI[1]][1]])

    dxMin[0] = np.Inf # Temp fix for defining surface position at corners correctly in exmplae problem

    minI = np.argmin(np.absolute(dxMin))
    
    surfacePos = np.array([0.0, 0.0, 0.0])
    surfacePos[minI] = dx[dxI[minI]][minI]

    # print(x, x0, x1,dx, dxMin, minI, surfacePos)

    return surfacePos

  def surfaceArea(self, pt, px, py, pz):
    if abs(pt[1] - self.x1[1]) < px or abs(pt[1] - self.x0[1]) < px:
        return px*pz
    else:
        return 0.0

  def xMin(self):
    # This might be used for binning objects
    return min(self.x0[0], self.x1[0])

  def xMax(self):
    # This might be used for binning objects
    return max(self.x0[0], self.x1[0])


#############################################
class rotatedBox:
    'Geometry object for creating a grid-aligned box defined by two corners.'

    def __init__(self, name, x0, x1, v, mat, group, damage, rotationCenter, rotationMatrix):
        self.box = box('box', x0, x1, v, mat, group, damage)
        self.rotationCenter = np.array(rotationCenter)
        self.rotationMatrix = np.array(rotationMatrix)
        self.v = v
        self.mat = mat
        self.group = group
        self.damage = damage

    def isInterior(self, pt):
        pt = np.array(pt)
        ptr = np.matmul(np.transpose(self.rotationMatrix), pt - self.rotationCenter) + self.rotationCenter

        return self.box.isInterior(ptr)

    def isSurface(self, pt, l):
        pt = np.array(pt)
        ptr = np.matmul(np.transpose(self.rotationMatrix), pt - self.rotationCenter) + self.rotationCenter

        surfaceFlag, surfaceNormal = self.box.isSurface(ptr, l)

        if surfaceFlag != 0:
            cf = self.rotationMatrix  # compute cofactor to rotate surfaceNormal
            surfaceNormal = np.matmul(cf, surfaceNormal)

        return surfaceFlag, surfaceNormal

    def surfacePosition(self, pt):
        pt = np.array(pt)
        ptr = np.matmul(np.transpose(self.rotationMatrix), pt - self.rotationCenter) + self.rotationCenter

        surfacePos = self.box.surfacePosition(ptr)

        return np.matmul(self.rotationMatrix, surfacePos)

    def surfaceArea(self, pt, px, py, pz):
        return 0.0

    def xMin(self):
        # This might be used for binning objects
        return min(self.x0[0], self.x1[0])

    def xMax(self):
        # This might be used for binning objects
        return max(self.x0[0], self.x1[0])

#############################################
class polygon:
  'Geometry object for creating a polygon described by ordered vertices.'
  def __init__(self,name,plist,v,mat,group,damage):
    self.name = name
    self.type = 'polygon'
    self.plist = np.array(plist)
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage

  def ccw(self,x1,x2,x3):
    return (x3[1]-x1[1])*(x2[0]-x1[0]) - (x2[1]-x1[1])*(x3[0]-x1[0]) > -10**-16

  def intersect(self,A,B,C,D):
    return self.ccw(A,C,D) != self.ccw(B,C,D) and self.ccw(A,B,C) != self.ccw(A,B,D)

  def isInside(self,vertices,point):
    v_arr = np.array(vertices)
    xmin = min(v_arr[:,0])
    ymin = min(v_arr[:,1])
    xmax = max(v_arr[:,0])
    ymax = max(v_arr[:,1])
    dx = xmax - xmin
    dy = ymax - ymin
    outside = [xmin,ymin]
    outside[0] -= dx
    outside[1] -= dy
    nv = v_arr.shape[0]
    nIntersections = 0
    for i in range(nv):
      p1 = vertices[i]
      p2 = vertices[np.mod(i+1,nv)]
      if(self.intersect(p1,p2,outside,point)):
        nIntersections += 1
    if np.mod(nIntersections,2)==0:
      return False
    else:
      return True

  def isInterior(self,point):
    if(self.isInside(self.plist,point)):
      return True
    else:
      return False

  def isSurface(self,point,t):
    vertices = self.plist
    inner = vertices.copy()
    nv = inner.shape[0]
    for i in range(nv):
      im1 = np.mod(i-1,nv)
      ip1 = np.mod(i+1,nv)
      n1 = np.array(vertices[ip1]-vertices[i])
      n1 = n1/np.sqrt(n1[0]*n1[0]+n1[1]*n1[1])
      n2 = np.array(vertices[im1]-vertices[i])
      n2 = n2/np.sqrt(n2[0]*n2[0]+n2[1]*n2[1])
      theta = np.arccos(n1[0]*n2[0]+n1[1]*n2[1])
      n3 = np.array(n1+n2)
      n3 = n3/np.sqrt(n3[0]*n3[0]+n3[1]*n3[1])
      h = t*np.sqrt(1.0 + (1.0/np.tan(theta/2.0))**2.0 )
      testpoint = inner[i] + h*n3
      sign = 1.0 if self.isInterior(testpoint) else -1.0
      inner[i] = inner[i] + sign*h*n3
    if(not(self.isInside(inner,point)) and self.isInterior(point)):
      return defaultSurfaceFlag
    else:
      return 0

  def xMin(self):
    # This might be used for binning objects
    arr = np.array(self.plist)
    xmin = min(arr[:,0])
    return xmin

  def xMax(self):
    # This might be used for binning objects
    arr = np.array(self.plist)
    xmax = max(arr[:,0])
    return xmax

#############################################
class foam:
  'Geometry object for creating a grid-aligned box defined by two corners with spherical pores.'
  def __init__(self,name,x0,x1,pores,v,mat,group,damage):
    # This fills a box with spherical pores defined by the array pores
    # assumes this is a numpy array [[r,x,y,z],[r,x,y,z],...]
    # The surface flags should be correctly set and the searching of pores
    # for isInterior and isSurface use a KD-tree. which should be faster than
    # using particle file writer with each pore as its own object.
    self.name = name
    self.type = 'box'
    self.x0 = np.array(x0)
    self.x1 = np.array(x1)
    self.pores = pores # assumes this is a numpy array [[r,x,y,z],[r,x,y,z],...]
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage

    pts = []
    for p in pores:
      pts.append([ p[1], p[2], p[3] ])
    # neighbor list for points:
    self.kdt = KDTree(pts, leaf_size=np.ceil(len(pts)/2), metric='euclidean')

  def isInterior(self,pt):
    # is the point within the object    
    x = np.array(pt)
    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])
    if ( (xmin<=x[0] and x[0]<xmax) and (ymin<=x[1] and x[1]<ymax) and (zmin<=x[2] and x[2]<zmax) ):
      # in the box
      test = True

      dist, index = self.kdt.query(x.reshape(1,-1), k=5)
      for i in index[0]:
        p=self.pores[i]
        if ( ( x[0] - p[1] )**2 + ( x[1] - p[2] )**2 + ( x[2] - p[3] )**2 < p[0]**2):
          # in any pore
          test = False
      return test
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    x = np.array(pt)
    xmin = min(self.x0[0],self.x1[0]) + l
    xmax = max(self.x0[0],self.x1[0]) - l
    ymin = min(self.x0[1],self.x1[1]) + l
    ymax = max(self.x0[1],self.x1[1]) - l
    zmin = min(self.x0[2],self.x1[2]) + l
    zmax = max(self.x0[2],self.x1[2]) - l
    
    if ( (xmin>=x[0] or x[0]>xmax) or (ymin>=x[1] or x[1]>ymax) or (zmin>=x[2] or x[2]>zmax) ):
      # surface of the box
      return defaultSurfaceFlag
    else:
      # interior of the box
      test = 0

      dist, index = self.kdt.query(x.reshape(1,-1), k=5)
      for i in index[0]:
        p=self.pores[i]
        if ( ( x[0] - p[1] )**2 + ( x[1] - p[2] )**2 + ( x[2] - p[3] )**2 < ( p[0] + l )**2):
          # on the surface of a pore
          test = defaultSurfaceFlag
      return test

  def xMin(self):
    # This might be used for binning objects
    return min(self.x0[0], self.x1[0])

  def xMax(self):
    # This might be used for binning objects
    return max(self.x0[0], self.x1[0])

#############################################
# def voronoi_volumes(points, v0):
#   # computes volumes for a voronoi tesselation of given points
#   # or returns v0 for edge points without enough data to compute.
#     v = Voronoi(points)
#     vol = np.zeros(v.npoints)
#     for i, reg_num in enumerate(v.point_region):
#         indices = v.regions[reg_num]
#         if -1 in indices: # some regions can be opened
#             vol[i] = v0
#         else:
#             vol[i] = ConvexHull(v.vertices[indices]).volume
#     return vol

class voronoiWeibullBox:
  'Box wrapper for another object that will be used to assign voronoi-cell weibull distribution of strength scale'
  def __init__(self,name,subObject,x0,x1,flawSize, weibullVolume, weibullModulus, weibullSeed, vMin):
    self.object = subObject # all properties will be inhereted from this subObject, except getDamage() 
    self.name = name
    self.type = 'voronoiWeibullBox'
    self.x0 = np.array(x0) # should define box bigger than subObject
    self.x1 = np.array(x1) # should define box bigger than subObject
    self.v = subObject.v
    self.mat = subObject.mat
    self.group = subObject.group
    self.damage = subObject.damage

    np.random.seed(weibullSeed)

    # create random points used for voronoi tesselation:
    boxVolume = (x1[0]-x0[0])*(x1[1]-x0[1])*(x1[2]-x0[2])
    flawVolume = flawSize**3
    npts = int(np.ceil(boxVolume/flawVolume))
    xmin = x0[0]
    xmax = x1[0]
    ymin = x0[1]
    ymax = x1[1]
    zmin = x0[2]
    zmax = x1[2]
    pts=[]
    for i in range(0,npts):
      pt = [ np.random.uniform(xmin,xmax),np.random.uniform(ymin,ymax),np.random.uniform(zmin,zmax) ]
      pts.append(pt)

    # neighbor list for points:
    self.kdt = KDTree(pts, leaf_size=np.ceil(len(pts)/2), metric='euclidean')
    #voronoi tesselation of points
    vor=Voronoi(pts)
    # average volume, will be assigned to edge cells.
    v0=(xmax-xmin)*(ymax-ymin)*(zmax-zmin)/npts

    # compute volume of each voronoi cell
    vol = np.zeros(vor.npoints)
    for i, reg_num in enumerate(vor.point_region):
      indices = vor.regions[reg_num]
      if ( (-1 in indices) or ( vor.vertices[vor.regions[i]].shape[0] < 1 ) ): # some regions can be opened
        vol[i] = v0
      else:
        vol[i] = ConvexHull(vor.vertices[indices]).volume

        numInteriorVertices = 0
        numVertices = vor.vertices[vor.regions[i]].shape[0]
        #
        for v in vor.vertices[indices]:
          if subObject.isInterior(v):
            numInteriorVertices += 1
        #
        #print("v = ",vol[i],", nI = ",numInteriorVertices,", nV = ",numVertices)
        vol[i] = vol[i]*numInteriorVertices/numVertices

      vol[i] = max( vol[i], vMin )

    # define the value of strength scale that will be assigned to each cell's particles.
    cellStrengthScale=[]
    for i in range(0,npts):
      s = ( ( weibullVolume/vol[i] )*( np.log( np.random.uniform(1e-20,1.0) )/np.log(0.5) ) )**(1.0/weibullModulus)
      cellStrengthScale.append(s)

    self.cellStrengthScale = np.array(cellStrengthScale)

  def getMat(self,pt):
    return self.object.getMat(pt)

  def isInterior(self,pt):
    return self.object.isInterior(pt)

  def isSurface(self,pt,l):
    return self.object.isSurface(pt,l)

  def getStrengthScale(self,pt):
    x = np.array(pt)
    dist, index = self.kdt.query(x.reshape(1,-1), k=1)
    strengthScale = self.cellStrengthScale[index[0][0]]
    return strengthScale

  def xMin(self):
    # This might be used for binning objects
    return self.object.xMin()

  def xMax(self):
    # This might be used for binning objects
    return self.object.xMax()

class voronoiWeibullBox2D:
  'PlaneStrain Box wrapper for another object that will be used to assign voronoi-cell weibull distribution of strength scale'
  def __init__(self,name,subObject,x0,x1,flawSize, weibullVolume, weibullModulus, weibullSeed, thickness, vMin):
    
    # this is for use in 2D plane strain, where volume is determined based on input
    # thickness.  
    self.object = subObject # all properties will be inhereted from this subObject, except getDamage() 
    self.name = name
    self.type = 'voronoiWeibullBox'
    self.x0 = np.array(x0) # [x0, y0] should define box bigger than subObject
    self.x1 = np.array(x1) # [x1, y1] should define box bigger than subObject
    self.v = subObject.v
    self.mat = subObject.mat
    self.group = subObject.group
    self.damage = subObject.damage

    np.random.seed(weibullSeed)

    # create random points used for voronoi tesselation:
    boxArea = (x1[0]-x0[0])*(x1[1]-x0[1])
    flawArea = flawSize**2
    npts = int(np.ceil(boxArea/flawArea))
    xmin = x0[0]
    xmax = x1[0]
    ymin = x0[1]
    ymax = x1[1]

    pts=[]
    for i in range(0,npts):
      pt = [ np.random.uniform(xmin,xmax),np.random.uniform(ymin,ymax),0.0 ]
      pts.append(pt)
    pts=np.array(pts)

    # neighbor list for points:
    self.kdt = KDTree(pts, leaf_size=np.ceil(len(pts)/2), metric='euclidean')

    # 2-D voronoi tesselation of points
    vor=Voronoi(pts[:,0:2])
    # average volume, will be assigned to edge cells.
    A0=(xmax-xmin)*(ymax-ymin)/npts

    # compute volume of each voronoi cell
    vol = np.zeros(vor.npoints)
    for i, reg_num in enumerate(vor.point_region):
      indices = vor.regions[reg_num]
      if ( (-1 in indices) or ( vor.vertices[vor.regions[i]].shape[0] < 1 ) ): # some regions can be opened
        area = A0
      else:
        area = ConvexHull(vor.vertices[indices]).volume
        numInteriorVertices = 0
        numVertices = vor.vertices[vor.regions[i]].shape[0]
        #
        for v in vor.vertices[indices]:
          v3d = np.array([v[0],v[1],0.0])
          if subObject.isInterior(v3d):
            numInteriorVertices += 1
        #
        #print("v = ",vol[i],", nI = ",numInteriorVertices,", nV = ",numVertices)
        area = area*numInteriorVertices/numVertices

      vol[i] = max( area*thickness, vMin )

    # define the value of strength scale that will be assigned to each cell's particles.
    cellStrengthScale=[]
    for i in range(0,npts):
      s = ( ( weibullVolume/vol[i] )*( np.log( np.random.uniform(1e-20,1.0) )/np.log(0.5) ) )**(1.0/weibullModulus)
      cellStrengthScale.append(s)

    self.cellStrengthScale = np.array(cellStrengthScale)

  def getMat(self,pt):
    return self.object.getMat(pt)

  def isInterior(self,pt):
    return self.object.isInterior(pt)

  def isSurface(self,pt,l):
    return self.object.isSurface(pt,l)

  def getStrengthScale(self,pt):
    x = np.array(pt)
    dist, index = self.kdt.query(x.reshape(1,-1), k=1)
    strengthScale = self.cellStrengthScale[index[0][0]]
    return strengthScale

  def xMin(self):
    # This might be used for binning objects
    return self.object.xMin()

  def xMax(self):
    # This might be used for binning objects
    return self.object.xMax()

class gridWeibullBox2D:
  'same as voronoiWeibullBox2D, but uses a uniform distribution of points.'
  def __init__(self,name,subObject,x0,x1,flawSize, weibullVolume, weibullModulus, weibullSeed, thickness, vMin):
    
    # this is for use in 2D plane strain, where volume is determined based on input
    # thickness.  
    self.object = subObject # all properties will be inhereted from this subObject, except getDamage() 
    self.name = name
    self.type = 'gridWeibullBox'
    self.x0 = np.array(x0) # [x0, y0] should define box bigger than subObject
    self.x1 = np.array(x1) # [x1, y1] should define box bigger than subObject
    self.v = subObject.v
    self.mat = subObject.mat
    self.group = subObject.group
    self.damage = subObject.damage

    np.random.seed(weibullSeed)

    # create random points used for voronoi tesselation:
    boxArea = (x1[0]-x0[0])*(x1[1]-x0[1])
    flawArea = flawSize**2

    xmin = x0[0]
    xmax = x1[0]
    ymin = x0[1]
    ymax = x1[1]

    pts=[]

    for x in np.linspace(xmin,xmax,np.round( (xmax-xmin)/flawSize ) ):
      for y in np.linspace(ymin,ymax,np.round( (ymax-ymin)/flawSize ) ):
        pt = [ x, y, 0.0 ]
        pts.append(pt)
     
    pts=np.array(pts)
    npts=pts.shape[0]

    # neighbor list for points:
    self.kdt = KDTree(pts, leaf_size=np.ceil(len(pts)/2), metric='euclidean')

    # 2-D voronoi tesselation of points
    vor=Voronoi(pts[:,0:2])
    # average volume, will be assigned to edge cells.
    A0=(xmax-xmin)*(ymax-ymin)/npts

    # compute volume of each voronoi cell
    vol = np.zeros(vor.npoints)
    for i, reg_num in enumerate(vor.point_region):
      indices = vor.regions[reg_num]
      if ( (-1 in indices) or ( vor.vertices[vor.regions[i]].shape[0] < 1 ) ): # some regions can be opened
        area = A0
      else:
        area = ConvexHull(vor.vertices[indices]).volume
        numInteriorVertices = 0
        numVertices = vor.vertices[vor.regions[i]].shape[0]
        #
        for v in vor.vertices[indices]:
          v3d = np.array([v[0],v[1],0.0])
          if subObject.isInterior(v3d):
            numInteriorVertices += 1
        #
        #print("v = ",vol[i],", nI = ",numInteriorVertices,", nV = ",numVertices)
        area = area*numInteriorVertices/numVertices

      vol[i] = max( area*thickness, vMin )

    # define the value of strength scale that will be assigned to each cell's particles.
    cellStrengthScale=[]
    for i in range(0,npts):
      s = ( ( weibullVolume/vol[i] )*( np.log( np.random.uniform(1e-20,1.0) )/np.log(0.5) ) )**(1.0/weibullModulus)
      cellStrengthScale.append(s)

    self.cellStrengthScale = np.array(cellStrengthScale)

  def getMat(self,pt):
    return self.object.getMat(pt)

  def isInterior(self,pt):
    return self.object.isInterior(pt)

  def isSurface(self,pt,l):
    return self.object.isSurface(pt,l)

  def getStrengthScale(self,pt):
    x = np.array(pt)
    dist, index = self.kdt.query(x.reshape(1,-1), k=1)
    strengthScale = self.cellStrengthScale[index[0][0]]
    return strengthScale

  def xMin(self):
    # This might be used for binning objects
    return self.object.xMin()

  def xMax(self):
    # This might be used for binning objects
    return self.object.xMax()

#############################################
class twoFieldBox:
  'Geometry object for creating a grid-aligned box defined by two corners with random assignment of group 1 or group 2.'
  def __init__(self,name,x0,x1,v,mat,group1,group2,damage):
    self.name = name
    self.type = 'box'
    self.x0 = np.array(x0)
    self.x1 = np.array(x1)
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = [group1,group2]
    self.damage = damage

  def isInterior(self,pt):
    # is the point within the object    
    x = np.array(pt)
    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])
    if ( (xmin<=x[0] and x[0]<xmax) and (ymin<=x[1] and x[1]<ymax) and (zmin<=x[2] and x[2]<zmax) ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    x = np.array(pt)
    xmin = min(self.x0[0],self.x1[0]) + l
    xmax = max(self.x0[0],self.x1[0]) - l
    ymin = min(self.x0[1],self.x1[1]) + l
    ymax = max(self.x0[1],self.x1[1]) - l
    zmin = min(self.x0[2],self.x1[2]) + l
    zmax = max(self.x0[2],self.x1[2]) - l
    if ( (xmin>=x[0] or x[0]>xmax) or (ymin>=x[1] or x[1]>ymax) or (zmin>=x[2] or x[2]>zmax) ):
      return defaultSurfaceFlag
    else:
      return 0  

  def xMin(self):
    # This might be used for binning objects
    return min(self.x0[0], self.X1[0])

  def xMax(self):
    # This might be used for binning objects
    return max(self.x0[0], self.X1[0])

#############################################
class cylinder:
  'Geometry object for creating a cylinder defined by two points and radius'
  def __init__(self,name,x1,x2,r,v,mat,group,damage):
    self.name = name
    self.type = 'cylinder'    
    self.x1 = np.array(x1)
    self.x2 = np.array(x2)
    self.r = r
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage
    self.center = self.x2-self.x1 # Temporary to get surface normals of 2D disks

    # for now we will default to cylinder having axis-aligned fiber direction,
    # may make more general later, but this allows for backwards compatibility
    h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    self.matDir = (self.x2-self.x1)/h   # unit vector along cylinder axis

  def isInterior(self,pt):    
    x = np.array(pt)
    h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    a = (self.x2-self.x1)/h             # unit vector along cylinder axis
    z = np.dot(x-self.x1,a)                                   # z-coordinate of test point
    r = np.linalg.norm( (x-self.x1) - np.dot(x-self.x1,a)*a )  # r coordinate of test point

    if ( (z>=0 and z<h) and (r < self.r) ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface? assumes particle is already interior
    x = np.array(pt)
    h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    a = (self.x2-self.x1)/h             # unit vector along cylinder axis
    z = np.dot(x-self.x1,a)                                   # z-coordinate of test point
    r = np.linalg.norm( (x-self.x1) - np.dot(x-self.x1,a)*a )  # r coordinate of test point

    if ( (z<l or z>(h-l)) or (r > self.r-l) ):
      surfaceNormal = x - self.x1
      surfaceNormal[2] = 0
      surfaceNormal = surfaceNormal / np.linalg.norm(surfaceNormal)
      # print(pt, self.center, self.x2, self.x1, self.x2-self.x1)
      # print(surfaceNormal, np.linalg.norm(surfaceNormal))
      return defaultSurfaceFlag, surfaceNormal
    else:
      return 0, None

  def xMin(self):
    # This might be used for binning objects
    return min(self.x1[0], self.x2[0])-self.r

  def xMax(self):
    # This might be used for binning objects
    return max(self.x1[0], self.x2[0])+self.r

#############################################
class whiskers:
  'Geometry object for creating a cylinder defined by two points and radius'
  def __init__(self,name,x1,x2,r,numWhiskers,volFracWhiskers,v,mat,group,damage):
    self.name = name
    self.type = 'cylinder'    
    self.x1 = np.array(x1)
    self.x2 = np.array(x2)
    self.r = r
    self.nw = numWhiskers
    self.VF = volFracWhiskers
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage

    self.h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    self.e1 = (self.x2-self.x1)/self.h             # unit vector along cylinder axis
    # random vectors orthogonal to e1 and each other
    self.e2 = np.random.rand(3)
    self.e2 = self.e2 - np.dot(self.e2,self.e1)*self.e1
    self.e2 = self.e2/np.linalg.norm(self.e2)
    self.e3 = np.cross(self.e1,self.e2)



  def isInterior(self,pt):    
    x = np.array(pt)
    z = np.dot(x-self.x1,self.e1)                                   # z-coordinate of test point
    pMinusX1InPlane = (x-self.x1) - np.dot(x-self.x1,self.e1)*self.e1
    r = np.linalg.norm( pMinusX1InPlane )  # r coordinate of test point
    

    ex = np.dot( self.e2,pMinusX1InPlane )
    ey = np.dot( self.e3,pMinusX1InPlane )

    theta = np.arctan2(ex,ey)

    dtheta = 2.0*np.pi/self.nw

    if ( ( (z>=0 and z<self.h) and (r < self.r) ) and ( (theta%dtheta)/dtheta < self.VF) ):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface? assumes particle is already interior
    x = np.array(pt)
    z = np.dot(x-self.x1,self.e1)                                   # z-coordinate of test point
    pMinusX1InPlane = (x-self.x1) - np.dot(x-self.x1,self.e1)*self.e1
    r = np.linalg.norm( pMinusX1InPlane )  # r coordinate of test point

    if ( (z<l or z>(self.h-l)) or (r > self.r-l) ):
      return defaultSurfaceFlag
    else:
      return 0

  def xMin(self):
    # This might be used for binning objects
    return min(self.x1[0], self.x2[0])-self.r

  def xMax(self):
    # This might be used for binning objects
    return max(self.x1[0], self.x2[0])+self.r

#############################################
class toroid:
  'Geometry object for creating a circular toroid defined by a point, direction ring radius and revolved circle radius'
  def __init__(self,name,x0,n,r1,r2,v,mat,group,damage):
    self.name = name
    self.type = 'toroid'    
    self.x0 = np.array(x0) # center point
    self.n = np.array(n) # axes unit vector
    self.r1 = r1         # ring radius
    self.r2 = r2         # revolved circle radius
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage

  def isInterior(self,pt):    
    x = np.array(pt)
    a = x - self.x0
    z = a.dot(self.n)
    r = np.linalg.norm(a - z*self.n)

    if (np.sqrt( (r-self.r1)*(r-self.r1) + z*z ) < self.r2):
      return True
    else:
      return False

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface? assumes particle is already interior
    x = np.array(pt)
    a = x - self.x0
    z = a.dot(self.n)
    r = np.linalg.norm(a - z*self.n)

    if (np.sqrt( (r-self.r1)*(r-self.r1) + z*z ) > self.r2 - l):
      return defaultSurfaceFlag
    else:
      return 0

  def xMin(self):
    # This might be used for binning objects
    return self.x0[0]-self.r1-self.r2

  def xMax(self):
    # This might be used for binning objects
    return self.x0[0]+self.r1+self.r2

#############################################
class fill:
  'Geometry object for creating particles to fill the background.'
  def __init__(self,name,mat,group,damage):
    print('Warning: Fill should only be used with untransformed objects.')
    self.v = np.array([0.,0.,0.])
    self.name = name
    self.mat = mat
    self.group = group
    self.damage = damage

  def isInterior(self,x):
    return True

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    return 0

  def xMin(self):
    # This might be used for binning objects
    return -1.0e9

  def xMax(self):
    # This might be used for binning objects
    return 1.0e9

#############################################

# TODO: add a v0, v1 for specifying an initial velocity gradient.
# 

class VCCTL:
  'Geometry object for creating a 3D object from a VCCTL voxelized dataset'
  def __init__(self,name,data,ni,nj,nk,x0,x1,v,map,group,damage):

    # data        text file with header stripped, containing one line per voxel and integer
    #             values indicating phase
    # ni, nj, nk  number of voxels in the x,y, and z directions.
    # x0, x1      coordinates of the -,+ corners of the object in the domain
    # map         dictionary of mappings from index to mat#: dict([(1, 2), (3, 4)])
    # mat         default material
    # group       contact group
    
    self.name = name
    self.type = 'VCCTL'    
    self.data = data
    self.ni = ni
    self.nj = nj
    self.nk = nk
    self.x0 = np.array(x0)
    self.x1 = np.array(x1)
    self.v = v if callable(v) else np.array(v)
    self.map = map
    self.group = group
    self.damage = damage

  def isInterior(self,pt):
    # is the point within the object
    x = np.array(pt)

    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])
    if ( (xmin<=x[0] and x[0]<xmax) and (ymin<=x[1] and x[1]<ymax) and (zmin<=x[2] and x[2]<zmax) ):
      i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      k = int( np.floor( self.nk*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
      i=max(min(self.ni-1,i),0)
      j=max(min(self.nj-1,j),0)
      k=max(min(self.nk-1,k),0)

      n = i*self.nj*self.nk + j*self.nk + k

      #if(n>=len(self.data)):
      #  print('n = ',n,', x0=',self.x0,', x1=',self.x1,',x=',x,'[i.j.k]=',[i,j,k],',[ni,nj,nk]=',[self.ni,self.nj,self.nk])
      index = int(self.data[n])
      mat = self.map.get(index)

      if (mat >-1):
        return True   # interior to the domain and not porosity
      else:
        return False  # internal porosity
    else:
      return False    # outside of domain

  def getMat(self,pt):
    # assumes the point within the object
    x = np.array(pt)

    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])

    i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
    j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
    k = int( np.floor( self.nk*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
    i=max(min(self.ni-1,i),0)
    j=max(min(self.nj-1,j),0)
    k=max(min(self.nk-1,k),0)
    n = i*self.nj*self.nk + j*self.nk + k

    index = self.data[n]
    mat = self.map[index]
    return mat # will error if map doesn't contain index as key

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    x = np.array(pt)
    xmin = min(self.x0[0],self.x1[0]) + l
    xmax = max(self.x0[0],self.x1[0]) - l
    ymin = min(self.x0[1],self.x1[1]) + l
    ymax = max(self.x0[1],self.x1[1]) - l
    zmin = min(self.x0[2],self.x1[2]) + l
    zmax = max(self.x0[2],self.x1[2]) - l
    if ( (xmin>=x[0] or x[0]>xmax) or (ymin>=x[1] or x[1]>ymax) or (zmin>=x[2] or x[2]>zmax) ):
      return defaultSurfaceFlag
    else:
      # check if point is next to porosity
      i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      k = int( np.floor( self.nk*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
      i=max(min(self.ni-1,i),0)
      j=max(min(self.nj-1,j),0)
      k=max(min(self.nk-1,k),0)
      
      imin = int( max(0,i-1) )
      imax = int( min(self.ni-1,i+1) )
      jmin = int( max(0,j-1) )
      jmax = int( min(self.nj-1,j+1) )
      kmin = int( max(0,k-1) )
      kmax = int( min(self.nk-1,k+1) )
      poreSurface = False
      for i in range(imin,imax):
        for j in range(jmin,jmax):
          for k in range(kmin,kmax):
            n = i*self.nj*self.nk + j*self.nk + k
            index = self.data[n]
            mat = self.map[index] # will error if map doesn't contain index as key

            if (mat <0):
              return defaultSurfaceFlag
      return 0  

#############################################


class CT:
  'Geometry object for creating a 3D object from a CT voxelized dataset'
  def __init__(self,name,data,ni,nj,nk,x0,x1,v,map,group,damage):

    # data are generated by a mathematica script, flattening a 3D array with
    # ordering [z,y,x]

    # This is similar (possibly identical) to the VCCTL object, but kept separate to allow
    # for specialization without compromising backwards compatibility.
    # data        text file with header stripped, containing one line per voxel and integer
    #             values indicating phase
    # ni, nj, nk  number of voxels in the x,y, and z directions.
    # x0, x1      coordinates of the -,+ corners of the object in the domain
    # map         dictionary of mappings from index to mat#: dict([(1, 2), (3, 4)])
    # mat         default material
    # group       contact group
    
    self.name = name
    self.type = 'CT'    
    self.data = data
    self.ni = ni
    self.nj = nj
    self.nk = nk
    self.x0 = np.array(x0)
    self.x1 = np.array(x1)
    self.v = v if callable(v) else np.array(v)
    self.map = map
    self.group = group
    self.damage = damage

  def isInterior(self,pt):
    # is the point within the object
    x = np.array(pt)

    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])
    if ( (xmin<=x[0] and x[0]<xmax) and (ymin<=x[1] and x[1]<ymax) and (zmin<=x[2] and x[2]<zmax) ):
      i = int( np.floor( self.ni*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      k = int( np.floor( self.nk*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )

      i=max(min(self.ni-1,i),0)
      j=max(min(self.nj-1,j),0)
      k=max(min(self.nk-1,k),0)

      n = i*self.nj*self.nk + j*self.nk + k

      #if(n>=len(self.data)):
      #  print('n = ',n,', x0=',self.x0,', x1=',self.x1,',x=',x,'[i.j.k]=',[i,j,k],',[ni,nj,nk]=',[self.ni,self.nj,self.nk])
      index = int(self.data[n])
      mat = self.map.get(index)

      if (mat >= 0):
        return True   # interior to the domain and not porosity
      else:
        return False  # internal porosity
    else:
      return False    # outside of domain

  def getMat(self,pt):
    # assumes the point within the object
    x = np.array(pt)

    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])

    i = int( np.floor( self.ni*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
    j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
    k = int( np.floor( self.nk*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
    i=max(min(self.ni-1,i),0)
    j=max(min(self.nj-1,j),0)
    k=max(min(self.nk-1,k),0)
    n = i*self.nj*self.nk + j*self.nk + k

    index = int(self.data[n])
    mat = self.map.get(index)
    return mat # will error if map doesn't contain index as key

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    x = np.array(pt)
    xmin = min(self.x0[0],self.x1[0]) + l
    xmax = max(self.x0[0],self.x1[0]) - l
    ymin = min(self.x0[1],self.x1[1]) + l
    ymax = max(self.x0[1],self.x1[1]) - l
    zmin = min(self.x0[2],self.x1[2]) + l
    zmax = max(self.x0[2],self.x1[2]) - l
    if ( (xmin>=x[0] or x[0]>xmax) or (ymin>=x[1] or x[1]>ymax) or (zmin>=x[2] or x[2]>zmax) ):
      return defaultSurfaceFlag
    else:
      # check if point is next to porosity
      i = int( np.floor( self.ni*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      k = int( np.floor( self.nk*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
      i=max(min(self.ni-1,i),0)
      j=max(min(self.nj-1,j),0)
      k=max(min(self.nk-1,k),0)
      
      imin = int( max(0,i-1) )
      imax = int( min(self.ni-1,i+1) )
      jmin = int( max(0,j-1) )
      jmax = int( min(self.nj-1,j+1) )
      kmin = int( max(0,k-1) )
      kmax = int( min(self.nk-1,k+1) )
      poreSurface = False
      for i in range(imin,imax):
        for j in range(jmin,jmax):
          for k in range(kmin,kmax):
            n = i*self.nj*self.nk + j*self.nk + k
            index = int(self.data[n])
            mat = self.map.get(index) # will error if map doesn't contain index as key
            if (mat <0):
              return defaultSurfaceFlag
      return 0  

#############################################

class bitmap:
  'Geometry object for creating a 2D object from an image file'
  def __init__(self,name,data,ni,nj,x0,x1,v,map,group,damage):

    # data        binary array
    # ni, nj, nk  number of voxels in the x,y, and z directions.
    # x0, x1      coordinates of the -,+ corners of the object in the domain
    # map         dictionary of mappings from index to mat#: dict([(1, 2), (3, 4)])
    # mat         default material
    # group       contact group
    
    self.name = name
    self.type = 'bitmap'    
    self.data = data
    self.ni = ni
    self.nj = nj
    self.x0 = np.array(x0)
    self.x1 = np.array(x1)
    self.v = v if callable(v) else np.array(v)
    self.map = map
    self.group = group
    self.damage = damage

  def isInterior(self,pt):
    # is the point within the object
    x = np.array(pt)
    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])
    if ( (xmin<=x[0] and x[0]<xmax) and (ymin<=x[1] and x[1]<ymax) and (zmin<=x[2] and x[2]<zmax) ):

      i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      #n = i*self.nj + j
      index = int(self.data[j,i])

      if (index >0):
        return True   # interior to the domain and not porosity
      else:
        return False  # internal porosity
    else:
      return False    # outside of domain

  def getMat(self,pt):
    # assumes the point within the object
    x = np.array(pt)

    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])

    i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
    j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
    #n = i*self.nj + j
    index = int(self.data[j,i])

    #index = self.data[n]
    return self.map[index] # will error if map doesn't contain index as key

  def isSurface(self,pt,l):
    # is the particle within distance l of the suface?
    x = np.array(pt)
    xmin = min(self.x0[0],self.x1[0]) + l
    xmax = max(self.x0[0],self.x1[0]) - l
    ymin = min(self.x0[1],self.x1[1]) + l
    ymax = max(self.x0[1],self.x1[1]) - l
    zmin = min(self.x0[2],self.x1[2]) + l
    zmax = max(self.x0[2],self.x1[2]) - l
    if ( (xmin>=x[0] or x[0]>xmax) or (ymin>=x[1] or x[1]>ymax) or (zmin>=x[2] or x[2]>zmax) ):
      return defaultSurfaceFlag
    else:
      # check if point is next to porosity
      i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )

      imin = int( max(0,i-1) )
      imax = int( min(self.ni-1,i+1) )
      jmin = int( max(0,j-1) )
      jmax = int( min(self.nj-1,j+1) )

      poreSurface = False
      for i in range(imin,imax):
        for j in range(jmin,jmax):
          #n = i*self.nj + j
          index = int(self.data[j,i])
          if (self.data[n]==0):
            return defaultSurfaceFlag
      return 0

################################################################################

class ctScene:
  'Geometry object for creating a 3D object from a CT voxelized dataset'
  def __init__(self,name,data,x0,x1,v,mat,group,damage):

    # data is a numpy array with dimensions equal to num voxels in y, x, and z
    # values indicate is void or material
    # x0, x1      coordinates of the -,+ corners of the object in the domain [x-, y-, z-]
    # mat         default material type
    # group       contact group
    # v           velocity vector
    # damage      scalar damage value for material
    # ni, nj, nk  number of voxels in the x,y, and z directions.

    self.name = name
    self.type = 'ctScene'    
    self.data = data
    self.ni = np.shape(data)[1]
    self.nj = np.shape(data)[0]
    self.nk = np.shape(data)[2]
    self.x0 = np.array(x0)
    self.x1 = np.array(x1)
    self.v = v if callable(v) else np.array(v)
    self.mat = mat
    self.group = group
    self.damage = damage
    self.dx = (self.x1[0]-self.x0[0])/self.ni
    self.dy = (self.x1[1]-self.x0[1])/self.nj
    self.dz = (self.x1[2]-self.x0[2])/self.nk

  def castToIJK(self, pt):

    stepsFromX0x = int(np.floor((pt[0] - self.x0[0])/self.dx))
    stepsFromX0y = int(np.floor((pt[1] - self.x0[1])/self.dy))
    stepsFromX0z = int(np.floor((pt[2] - self.x0[2])/self.dz))

    return [stepsFromX0x, stepsFromX0y, stepsFromX0z]

  def inBoundingBox(self, pt):
    if (self.x0[0] <= pt[0] <= self.x1[0]) and (self.x0[1] <= pt[1] <= self.x1[1]) and (self.x0[2] <= pt[2] <= self.x1[2]):
        return True
    else:
        return False

  def isInterior(self,pt):
    if self.inBoundingBox( pt):
        ijk = self.castToIJK(pt)
        return int(self.data[(self.nj-1) - int(ijk[1]), int(ijk[0]), int(ijk[2])] != 0)
    else: 
        return int(False)

  def isSurface(self,pt,l):
    for i in [-self.dx, 0, self.dx]:
        for j in [-self.dy, 0, self.dy]:
            for k in [-self.dz, 0, self.dz]:
                if not self.isInterior(pt+[i,j,k]):
                    return int(True)    

    return int(False)