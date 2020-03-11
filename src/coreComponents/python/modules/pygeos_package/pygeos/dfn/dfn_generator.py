"""
DFN generator for GEOS.
Original author: Pengcheng Fu, Version 2.0 (Oct 2017)
Ported version: Christopher Sherman (Mar 2020)
"""

import sys
import numpy as np
from lxml import etree
from scipy.stats import truncnorm
from pygeos.dfn import vtk_writer


def print_grid_dimensions():
  """
  Print the available grid dimensions for a DFN
    0 means the xz plane
    2 means the yz plane
    1 and 3 are the planes in between.
    Note: for gridDimension 0 and 2, only the box mode is valid.
  """
  print('')
  print('  y')
  print('  ^')
  print('  |\\3  /1')
  print('  | \\ /')
  print('2 |  x')
  print('  | / \\')
  print('  |/   \\')
  print('  +------>x')
  print('     0')
  print('')


def BoundedPowerLaw(xmin, xmax, power):
  """
  Truncated power law distribution (http://mathworld.wolfram.com/RandomNumber.html)

  @arg xmin Minimum value
  @arg xmax Maximum value
  @arg power Distribution exponent
  """
  x = (np.random.random() * (xmax**(power + 1) - xmin**(power + 1)) + xmin**(power + 1))**(1.0 / (power + 1))
  return x


def TruncatedNorm(mean, stdev, xmin, xmax):
  """
  Truncated normal distribution

  @arg mean Mean
  @arg stdev Standard deviation
  @arg xmin Minimum value
  @arg xmax Maximum value
  """
  if (xmin == xmax):
    return mean
  else:
    a = (xmin - mean) / stdev
    b = (xmax - mean) / stdev
    return truncnorm.rvs(a, b) * stdev + mean


def Intersect(ci, li, ai, cj, lj, aj, margin):
  """
  Check for fracture intersections

  @arg ci Center of first fracture
  @arg li Length of first fracture
  @arg ai Angle of first fracture (radians, ccws from x-axis)
  @arg cj Center of second fracture
  @arg lj Length of second fracture
  @arg aj Angle of second fracture (radians, ccws from x-axis)
  @arg margin If fractures are within this distance, count them as intersecting
  """

  if ((ci[0] - cj[0])**2 + (ci[1] - cj[1])**2 > (li + lj)**2 / 4):
    return False
  xi0 = ci[0] - np.cos(ai) * (li + margin) * 0.5
  xi1 = ci[0] + np.cos(ai) * (li + margin) * 0.5
  yi0 = ci[1] - np.sin(ai) * (li + margin) * 0.5
  yi1 = ci[1] + np.sin(ai) * (li + margin) * 0.5
  xj0 = cj[0] - np.cos(aj) * (lj + margin) * 0.5
  xj1 = cj[0] + np.cos(aj) * (lj + margin) * 0.5
  yj0 = cj[1] - np.sin(aj) * (lj + margin) * 0.5
  yj1 = cj[1] + np.sin(aj) * (lj + margin) * 0.5

  Poly = [[xi0, yi0], [xi1, yi1], [xj0, yj0]]
  A0 = PolyArea(Poly)
  Poly = [[xi0, yi0], [xi1, yi1], [xj1, yj1]]
  A1 = PolyArea(Poly)
  Poly = [[xj0, yj0], [xj1, yj1], [xi0, yi0]]
  A2 = PolyArea(Poly)
  Poly = [[xj0, yj0], [xj1, yj1], [xi1, yi1]]
  A3 = PolyArea(Poly)
  if (A0 * A1 <= 0.0 and A2 * A3 <= 0.0):
    return True
  else:
    return False


def PolyArea(Poly):
  """
  Calculate polygon area

  @arg Poly Nx2 array of coordinate pairs
  """
  area = 0.0
  for i in range(-1, len(Poly) - 1):
    area += Poly[i][0] * Poly[i + 1][1] - Poly[i + 1][0] * Poly[i][1]
  return 0.5 * area


def EvenSpacing(x):
  """
  Check whether grid spacing within a mesh region are uniform

  @arg x Array of axis coordinates in region
  """
  intv = [abs(x[i + 1] - x[i]) for i in range(len(x) - 1)]
  meanIntv = np.mean(intv)

  if len(x) == 2:
    return meanIntv
  else:
    sdIntv = np.std(intv)
    print(sdIntv)

    if sdIntv <= 1e-6:
      return meanIntv
    else:
      return -meanIntv


def NearestIndex(a, x):
  """
  Find the index of the closest value in an array

  @arg a 1D array
  @arg x Target value
  """
  return np.argmin(abs(np.array(a) - x))


def BuildNodeset(parent, name, fractureCenter, fractureLength, fractureHeight, fractureAngle, margin):
  """
  Build a nodeset for a single fracture

  @arg parent The parent xml element for the new nodeset definition
  @arg name Name of the nodeset
  @arg fractureCenter Center of the nodeset
  @arg fractureLength Length of the nodeset
  @arg fractureHeight Height of the nodeset
  @arg fractureAngle Angle of the nodeset (radians, ccws from x-axis)
  @arg margin Margin of the nodeset definition
  """

  newSet = etree.Element("BoundedPlane", name=name)

  newSet.set('lengthVector', '%1.6e, %1.6e, 0' % (np.cos(fractureAngle), np.sin(fractureAngle)))
  newSet.set('widthVector', '0, 0, 1')
  newSet.set('origin', '%1.6e, %1.6e, %1.6e' % (fractureCenter[0], fractureCenter[1], fractureCenter[2]))
  newSet.set('normal', '%1.6e, %1.6e, 0' % (np.cos(fractureAngle+0.5*np.pi), np.sin(fractureAngle+0.5*np.pi)))
  newSet.set('dimensions', '{%1.6e, %1.6e}' % (fractureLength+margin, fractureHeight+margin))

  parent.insert(-1, newSet)


def WriteDFNxml(inputFileName, outputFileName, fractureCenter, fractureLength, fractureHeight, fractureAngle, fractureSet, margin):
  """
  Write a new XML file with the generated DFN

  @arg inputFileName Name of the origianl XML file
  @arg outputFileName Name of the new XML file
  @arg fractureCenter Nx3 List of nodeset centers
  @arg fractureLength List of nodeset lengths
  @arg fractureHeight List of nodeset heights
  @arg fractureAngle List of nodeset angles (radians, ccws from x-axis)
  @arg fractureSet List of nodeset ID's
  @arg margin Margin for nodeset definitions
  """

  parser = etree.XMLParser(remove_comments=True, remove_blank_text=True)
  tree = etree.parse(inputFileName, parser=parser)
  root = tree.getroot()

  # Write fracture definitions
  GeometryBlock = root.findall('Geometry')[0]
  fracNames = ['Frac_%06d' % (ii) for ii in range(0, len(fractureCenter))]
  for ii in range(0, len(fractureCenter)):
    BuildNodeset(GeometryBlock, fracNames[ii], fractureCenter[ii], fractureLength[ii], fractureHeight[ii], fractureAngle[ii], margin)

  # Write an initial condition to mark the DFN faces
  InitDFN = etree.Element("FieldSpecification", name='init_dfn')
  InitDFN.set('fieldName', 'ruptureState')
  InitDFN.set('initialCondition', '1')
  InitDFN.set('objectPath', 'faceManager')
  InitDFN.set('scale', '1')
  InitDFN.set('setNames', '{' + ', '.join(fracNames) + '}')
  FieldSpecificationsBlock = root.findall('FieldSpecifications')[0]
  FieldSpecificationsBlock.insert(-1, InitDFN)

  tree.write(outputFileName, pretty_print=True)


def parse_attributes(source, targets, types, defaults):
  """
  Parse values from an XML node

  @arg source XML node
  @arg targets list of attribute names
  @arg types list of attribute types (used for string conversion)
  @arg defaults dict defining attribute default values
  """

  attribute_dict = {}
  for k, kt in zip(targets, types):
    if k in source.attrib:
      tmp = source.get(k)
      if '{' in tmp:
        attribute_dict[k] = [kt(x) for x in tmp[1:-1].split(',')]
      else:
        attribute_dict[k] = kt(tmp)

    elif k in defaults:
      attribute_dict[k] = defaults[k]

    else:
      raise Exception('Required attribute not supplied: %s' % (k))

  return attribute_dict


def parse_geosx_xml(inputFileName):
  """
  Parse a geosx format input file

  @arg inputFileName Name of the source xml document

  This function will attempt to read the following from file.
  Mandatory:
     name: string, name.
     mode: string, box or strip.
     gridDimension: integer, 0, 1, 2, or 3. 0 means the xz plane, 2 means the yz plane, 1 and 3 are the planes in between.
                    for gridDimension 0 and 2, only the box mode is valid.
     xmin, xmax: three real values separated by white spaces. They are coordinates in the unskewed coordinate system. For the box mode, they define the two corners of the box.
                 For the the strip mode, xmin[0] and xmax[2] define the left and right ends of the set's intersecton with y=skewCenter[1]
     coverageRatio: real, the target coverage ratio on each grid plane. A fracture plane typically would saturate at a coverage of 0.4, depending on the minGap value.
     minGap: integer, the minimum gap size, in terms of grids, between two neighbor fractures on the same grid plane.
     lmin, lmax, power: define the length distribution, it follows a power-law distribution.
     meanAR, stdevAR, minAR, maxAR: defines the aspect ratio distribution.
  Optional:
     seed: integer, >=0, random seed
     planeExtensionFactor: integer, default is 10. We extend the grid plane by n * grid spacing to avoid artefacts caused by truncating the domain
     leftExtension, 1 or 0, defualt 0. When we are in the strip mode, we apply this to the leftmost set to extend it to cover the triangle area that would be otherwise missed
     rightExtension, 1 or 0, same idea.
  """

  # Parse file
  results = {'dfn': {}}
  parser = etree.XMLParser(remove_comments=True, remove_blank_text=True)
  tree = etree.parse(inputFileName, parser=parser)
  root = tree.getroot()

  # Read mesh parameters
  meshNode = root.findall('Mesh')[0]
  intentalMeshNode = meshNode.findall('InternalMesh')[0]
  mesh_targets = ['xCoords', 'yCoords', 'zCoords', 'nx', 'ny', 'nz', 'skewCenter', 'skewAngle']
  mesh_types = [float, float, float, int, int, int, float, float]
  mesh_defaults = {'skewCenter': [0.0, 0.0, 0.0], 'skewAngle': 0.0}
  results['mesh'] = parse_attributes(intentalMeshNode, mesh_targets, mesh_types, mesh_defaults)

  # Calculate grid ticks (note: these are before skewing)
  for k in ['x', 'y', 'z']:
    c = results['mesh']['%sCoords' % (k)]
    n = results['mesh']['n%s' % (k)]

    tmp = [np.linspace(c[ii], c[ii+1], n[ii]+1) for ii in range(0, len(n))]
    results['mesh']['%sticks' % (k)] = np.unique(np.concatenate(tmp))

  # Find center
  results['mesh']['yCenterID'] = NearestIndex(results['mesh']['xticks'], results['mesh']['skewCenter'][1])

  # Parse DFN definitions
  dfn_targets = ['mode', 'gridDimension', 'xmin', 'xmax', 'coverageRatio', 'minGap', 'lmin',
                 'lmax', 'power', 'meanAR', 'stdevAR', 'maxAR', 'minAR', 'seed',
                 'planeExtensionFactor', 'leftExtension', 'rightExtension']
  dfn_types = [str, int, float, float, float, int, float,
               float, float, float, float, float, float, int,
               int, int, int]
  dfn_defaults = {'seed': -1, 'planeExtensionFactor': 10, 'leftExtension': 0, 'rightExtension': 0}

  for DFNSetsNode in root.findall('DFNSets'):
    for SetNode in DFNSetsNode.findall('DFNSet'):
      name = SetNode.get('name')
      results['dfn'][name] = parse_attributes(SetNode, dfn_targets, dfn_types, dfn_defaults)

  return results


def generate_from_xml(inputFileName, outputFileName, nIter=10000, margin=0.01, calculatePercolation=True):
  """
  This function will parse a geosx input file and parse any embedded DFN instructions.

  @arg inputFileName Input file name
  @arg outputFileName Output filename
  @arg nIter Max number of iterations for DFN construction
  @arg margin Margin for DFN sets
  @arg calculatePercolation Calculate percolation for generated DFN
  """

  # Setup
  dfn_config = parse_geosx_xml(inputFileName)
  fractures = {k: [] for k in ['length', 'height', 'angle', 'center', 'set', 'planeID', 'intersections']}
  planeID = -1
  setID = -1

  # Construct DFN
  for kd in dfn_config['dfn'].keys():
    gridDimension = dfn_config['dfn'][kd]['gridDimension']
    coverageRatio = dfn_config['dfn'][kd]['coverageRatio']
    minGap = dfn_config['dfn'][kd]['minGap']
    lmin = dfn_config['dfn'][kd]['lmin']
    lmax = dfn_config['dfn'][kd]['lmax']
    power = dfn_config['dfn'][kd]['power']
    meanAR = dfn_config['dfn'][kd]['meanAR']
    stdevAR = dfn_config['dfn'][kd]['stdevAR']
    maxAR = dfn_config['dfn'][kd]['maxAR']
    minAR = dfn_config['dfn'][kd]['minAR']
    seed = dfn_config['dfn'][kd]['seed']
    planeExtensionFactor = dfn_config['dfn'][kd]['planeExtensionFactor']
    leftExtension = dfn_config['dfn'][kd]['leftExtension']
    rightExtension = dfn_config['dfn'][kd]['rightExtension']
    mode = dfn_config['dfn'][kd]['mode']
    xmin = dfn_config['dfn'][kd]['xmin']
    xmax = dfn_config['dfn'][kd]['xmax']

    # Check inputs
    if mode == 'strip' and (gridDimension == 0 or gridDimension == 2):
      raise Exception(kd + ': For gridDimension = 0 or 2, there is not reason to use the strip mode.  Use box mode instead.')
    if mode == 'box' and (leftExtension > 0 or rightExtension > 0):
      raise Exception(kd + ': leftExtension and right Extension only apply to the strip mode, not the box mode.')

    # Set seed
    if (seed >= 0):
      np.random.seed(seed)

    # Grab a segment of the mesh
    subMesh = []
    for ii, kx in zip([0, 1, 2], ['x', 'y', 'z']):
      t = dfn_config['mesh']['%sticks' % (kx)]
      Ia = np.where((t >= xmin[ii]) & (t <= xmax[ii]))
      subMesh.append(t[Ia])

    if EvenSpacing(subMesh[0]) < 0.0 or EvenSpacing(subMesh[1]) < 0.0 or EvenSpacing(subMesh[2]) < 0.0:
      sys.exit('ERROR: mesh grids in the DFN box must be evenly. This is not satisfied in set ' + kd)

    # Select geometric parameters
    a0, a1, da, b0, b1, db, c0, c1, dc = subMesh[1][0], subMesh[1][-1], EvenSpacing(subMesh[1]), subMesh[2][0], subMesh[2][-1], EvenSpacing(subMesh[2]), subMesh[0][0], subMesh[0][-1], EvenSpacing(subMesh[0])
    if (gridDimension == 0):
      a0, a1, da, b0, b1, db, c0, c1, dc = subMesh[0][0], subMesh[0][-1], EvenSpacing(subMesh[0]), subMesh[2][0], subMesh[2][-1], EvenSpacing(subMesh[2]), subMesh[1][0], subMesh[1][-1], EvenSpacing(subMesh[1])

    if (gridDimension == 1 or gridDimension == 3):
      xID = NearestIndex(dfn_config['mesh']['xticks'][0], c0)
      if (xID + dfn_config['mesh']['yCenterID']) % 2 == 0:
        c0 -= dc
      dc *= 2

      nIntYPlus = max(0, round((a1 - dfn_config['mesh']['skewCenter'][1]) / da / 2) + 1)
      nIntYMinus = max(0, round((dfn_config['mesh']['skewCenter'][1] - a0) / da / 2) + 1)

      if mode == 'box':
        if (gridDimension == 1):
          c0 -= dc * nIntYPlus
          c1 += dc * nIntYMinus
        elif (gridDimension == 3):
          c0 -= dc * nIntYMinus
          c1 += dc * nIntYPlus

      if mode == 'strip' and leftExtension > 0:
        if (gridDimension == 1):
          c0 -= dc * nIntYPlus
        elif (gridDimension == 3):
          c0 -= dc * nIntYMinus
      if mode == 'strip' and rightExtension > 0:
        if (gridDimension == 1):
          c1 += dc * nIntYMinus
        elif (gridDimension == 3):
          c1 += dc * nIntYPlus

      a0 -= da * planeExtensionFactor  # Extension of the domain to avoid boundary effects
      a1 += da * planeExtensionFactor

    na = int((a1 - a0) / da)
    nb = int((b1 - b0) / db)
    nc = int(round(c1 - c0) / dc) + 1

    setID += 1
    totalCount = 0
    rejectionCount = 0
    occupiedCellCout = 0

    for ic in range(0, nc):
      planeID += 1
      c = c0 + ic * dc

      # strikeAngle is from x=0, counterclockwise; not the geological definition from north
      lengthProjFactor = 1.0
      if (gridDimension == 0):
        strikeAngle = 0.0
        lengthProjFactor = 1.0
      elif (gridDimension == 1):
        strikeAngle = np.arctan(1 / (dc / 2 / da - np.tan(np.radians(dfn_config['mesh']['skewAngle']))))
        lengthProjFactor = np.sin(strikeAngle)
      elif (gridDimension == 2):
        strikeAngle = np.pi / 2 + np.radians(dfn_config['mesh']['skewAngle'])
        lengthProjFactor = np.sin(strikeAngle)
      elif (gridDimension == 3):
        strikeAngle = np.pi - np.arctan(1 / (dc / 2 / da + np.tan(np.radians(dfn_config['mesh']['skewAngle']))))
        lengthProjFactor = abs(np.sin(strikeAngle))
      else:
        sys.exit('ERROR: gridDimension must be between 0 and 3.')

      # Generate the length array and sort it.  We place the big fractures first and small ones later.
      totalArea = 0.0

      lw = []
      area = []
      while (totalArea < (a1 - a0) * (b1 - b0) * coverageRatio):
        la = int(round(BoundedPowerLaw(lmin * lengthProjFactor, lmax * lengthProjFactor, power) / da)) * da  # This is the projected length
        aspectRatio = TruncatedNorm(meanAR, stdevAR, minAR, maxAR)
        # A fracture needs to be at least 2 elements wide to topologically split the mesh.
        w = max(2, int(round(la / lengthProjFactor * aspectRatio / db))) * db
        lw.append([la, w])
        area.append(la * w)
        totalArea += la * w  # Projected area

      # Sort length-width pair by area
      zipped = zip(area, lw)
      zipped = sorted(zipped)
      lw = [a[1] for a in zipped]
      lw[:] = lw[::-1]  # reverse the order

      backCells = np.ndarray(shape=(na, nb), dtype=int)
      backCells *= 0

      totalCount += len(lw)
      for i in range(0, len(lw)):
        placed = 0
        la = int(lw[i][0] / da)
        lb = int(round(lw[i][1] / db))
        count = 0
        while (placed == 0):
          count += 1

          if (count > nIter):
            print("rejecting ", i)
            rejectionCount += 1
            break
          ia = np.random.randint(0, na - int(lw[i][0] / da))
          ib = np.random.randint(0, nb - int(lw[i][1] / db))
          placed = 1

          if (np.sum(backCells[max(ia - minGap, 0):min(ia + la + minGap, na), max(ib - minGap, 0):min(ib + lb + minGap, nb)]) > 0):
            placed = 0

          if (placed == 1):
            backCells[ia:ia + la, ib:ib + lb] += np.random.randint(2, 10)

            fractures['length'].append(la * da / lengthProjFactor)
            fractures['height'].append(lb * db)
            fractures['angle'].append(strikeAngle)
            fractures['set'].append(kd)
            if (gridDimension >= 1):
              x = c
              y = a0 + ia * da + la * da * 0.5
              z = b0 + ib * db + lb * db * 0.5
              x += y / np.tan(strikeAngle)
            else:
              x = a0 + ia * da + la * da * 0.5
              y = c
              z = b0 + ib * db + lb * db * 0.5
              x -= (y - dfn_config['mesh']['skewCenter'][1]) * np.tan(np.radians(dfn_config['mesh']['skewAngle']))

            fractures['center'].append([x, y, z])
            fractures['planeID'].append(setID)
            fractures['intersections'].append(0)

      occupiedCellCout += np.count_nonzero(backCells)
    print('In set ' + kd + ' rejected ' + str(rejectionCount) + ' out of ' + str(totalCount) + ' fractures.' + 'We have ' + str(nc) + ' planes.')
    print('Set ' + kd + ' average coverage ratio is: ' + str(occupiedCellCout * 1.0 / na / nb / nc))

  print('Finished, we have ', len(fractures['center']), ' fractures')

  # Calculate percolation number
  if (calculatePercolation >= 1):
    nPairs = len(fractures['center']) * (len(fractures['center']) - 1) / 2
    print('Calculating percolation number, there are ' + str(nPairs) + ' to check')

    pairCount = 0
    print('Finished    \r')
    for i in range(0, len(fractures['center'])):
      for j in range(i + 1, len(fractures['center'])):
        pairCount += 1

        if (pairCount % (int(nPairs / 10)) == 0 or (pairCount % (int(nPairs / 100)) == 0 and pairCount < 0.1 * nPairs)):
          print(str(pairCount * 100 / nPairs).zfill(2) + '%\r')

        if (fractures['angle'][i] != fractures['angle'][j]):
          zA0 = fractures['center'][i][2] - fractures['height'][i] * 0.5
          zA1 = fractures['center'][i][2] + fractures['height'][i] * 0.5
          zB0 = fractures['center'][j][2] - fractures['height'][j] * 0.5
          zB1 = fractures['center'][j][2] + fractures['height'][j] * 0.5

          if (not(zA0 >= zB1 or zB0 >= zA1)):
            if Intersect(fractures['center'][i], fractures['length'][i], fractures['angle'][i], fractures['center'][j], fractures['length'][j], fractures['angle'][j], margin):
              fractures['intersections'][i] += 1
              fractures['intersections'][j] += 1
  else:
    print('Not calculating percolation number as instructed.')
  vtk_writer.WriteVTK('dfn_preview.vtk', fractures['center'], fractures['length'], fractures['height'], fractures['angle'], fractures['planeID'], fractures['intersections'])
  WriteDFNxml(inputFileName, outputFileName, fractures['center'], fractures['length'], fractures['height'], fractures['angle'], fractures['set'], margin)
  totalIntersection = np.sum(fractures['intersections'])
  print('The percolation number is ', totalIntersection * 1.0 / len(fractures['center']))


