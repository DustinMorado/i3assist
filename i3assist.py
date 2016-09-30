""" I3 helper module.

Written by Dustin R. Morado
Version: 0.1.0
Date: 27.09.2016

This is a module so that IPython and Python can be used to facilitate the use
of I3 for subtomogram averaging.
"""

import math
import numpy

class Euler:
  """ ZXZ extrinsic euler angle rotation object.

  Describes a rotation of coordinate axes by phi around Z, then theta around
  the new X axis, and finally phi around the new Z.
  """
  def __init__(self, phi=0.0, theta=0.0, psi=0.0, unit="deg"):
    """ Initializes euler object.

    Keyword arguments:
    phi -- First rotation of the Z coordinate axis.
    theta -- Second rotation of the new X coordinate axis.
    psi -- Third rotation of the new Z coordinate axis.
    unit -- Unit of rotations, either "deg" or "rad" for degrees or radians.
    """
    self.unit = unit
    self.phi = phi
    self.theta = theta
    self.psi = psi

  @property
  def unit(self):
    """Describes whether Euler angles are in degrees or radians."""
    return self.__unit

  @unit.setter
  def unit(self, value):
    if value == "deg" or value == "rad":
      self.__unit = value
    else:
      raise ValueError("Unit must be either \"deg\" or \"rad\".")

  @property
  def phi(self):
    """First rotation of the Z coordinate axis."""
    return self.__phi

  @phi.setter
  def phi(self, phi):
    """Set phi.

    Keyword arguments:
    phi -- Rotation
    """
    self.__phi = float(phi)

  @property
  def theta(self):
    """First rotation of the Z coordinate axis."""
    return self.__theta

  @theta.setter
  def theta(self, theta):
    """Set theta.

    Keyword arguments:
    theta -- Rotation
    """
    self.__theta = float(theta)

  @property
  def psi(self):
    """First rotation of the Z coordinate axis."""
    return self.__psi

  @psi.setter
  def psi(self, psi, unit="deg"):
    """Set psi.

    Keyword arguments:
    psi -- Rotation
    """
    self.__psi = float(psi)

  def degrees(self):
    """Convert Euler angles from radians to degrees."""
    if self.unit == "rad":
      self.unit = "deg"
      self.phi = math.degrees(self.phi)
      self.theta = math.degrees(self.theta)
      self.psi = math.degrees(self.psi)
  
  def radians(self):
    """Convert Euler angles from degrees to radians."""
    if self.unit == "deg":
      self.unit = "rad"
      self.phi = math.radians(self.phi)
      self.theta = math.radians(self.theta)
      self.psi = math.radians(self.psi)

  def copy(self):
    """Return copy of Euler angle object."""
    return Euler(phi=self.phi, theta=self.theta, psi=self.psi, unit=self.unit)

  def getVector(self):
    """Return numpy array of euler angles."""
    return numpy.array([self.phi, self.theta, self.psi], numpy.float_)

  def __str__(self):
    return "{: f}\t{: f}\t{: f}".format(self.phi, self.theta, self.psi)

  def getRotationMatrix(self):
    """Return Rotation Matrix object using the euler object's rotation."""
    return RotationMatrix(phi=self.phi, theta=self.theta, psi=self.psi,
                          unit=self.unit)

  def normalize(self):
    """Normalize Euler angles in the range needed by old I3."""
    rotationMatrix = self.getRotationMatrix()
    normalizedEuler = rotationMatrix.getEuler(self.unit)
    self.phi = normalizedEuler.phi
    self.theta = normalizedEuler.theta
    self.psi = normalizedEuler.psi

  def invert(self):
    """Invert Euler angles."""
    rotationMatrix = self.getRotationMatrix()
    rotationMatrix.invert()
    inverseEuler = rotationMatrix.getEuler(self.unit)
    self.phi = inverseEuler.phi
    self.theta = inverseEuler.theta
    self.psi = inverseEuler.psi

class RotationMatrix:
  """ ZXZ extrinsic rotation matrix object.

  Describes a rotation of coordinate axes by phi around Z, then theta around
  the new X axis, and finally phi around the new Z.
  """
  def __init__(self, phi=0., theta=0., psi=0., unit="deg"):
    """ Initializes rotation matrix object.

    Keyword arguments:
    phi -- First rotation of the Z coordinate axis.
    theta -- Second rotation of the new X coordinate axis.
    psi -- Third rotation of the new Z coordinate axis.
    unit -- Unit of rotations, either "deg" or "rad" for degrees or radians.
    """
    angles = numpy.array([phi, theta, psi], numpy.float_)
    if unit == "deg":
      angles = numpy.radians(angles)
    elif unit == "rad":
      pass
    else:
      raise ValueError("Unit must be either \"deg\" or \"rad\".")
    self.matrix = angles
  
  @property
  def matrix(self):
    """ZXZ Rotation matrix."""
    return self.__matrix

  @matrix.setter
  def matrix(self, angles):
    cos = numpy.cos(angles)
    sin = numpy.sin(angles)
    self.__matrix = numpy.array(
        [[ cos[2] * cos[0] - sin[2] * cos[1] * sin[0],
           cos[2] * sin[0] + sin[2] * cos[1] * cos[0],
           sin[2] * sin[1]],
         [-sin[2] * cos[0] - cos[2] * cos[1] * sin[0],
          -sin[2] * sin[0] + cos[2] * cos[1] * cos[0],
           cos[2] * sin[1]],
         [ sin[1] * sin[0], -sin[1] * cos[0], cos[1]]], numpy.float_)
  
  def getVector(self):
    """Return single linear array of rotation matrix."""
    return self.matrix.reshape(1,9)[0]

  def __str__(self):
    """Print rotation matrix in i3 format (one line a00 a01 a02 a10 ...)."""
    return ("{: f} " * 8 + "{: f}").format(*self.getVector())

  def invert(self):
    """Invert rotation matrix."""
    self.__matrix = self.__matrix.T

  def transpose(self):
    """Transpose (invert) rotation matrix."""
    self.__matrix = self.__matrix.T

  def getEuler(self, unit="deg"):
    """Returns Euler object in radians equivalent of Rotation Matrix."""
    absSineTheta = math.sqrt(self.__matrix[0,2] ** 2 + self.__matrix[1,2] ** 2)
    if not numpy.isclose(absSineTheta, 0.0):
      phi = math.atan2(self.__matrix[2,0], -self.__matrix[2,1])
      psi = math.atan2(self.__matrix[0,2], self.__matrix[1,2])
      
      if numpy.isclose(math.sin(phi), 0.0):
        thetaSign = math.copysign(1.0, (self.__matrix[1,2] / math.cos(psi)))
      else:
        thetaSign = math.copysign(1.0, (self.__matrix[0,2] / math.sin(psi)))
      theta = math.atan2(thetaSign * absSineTheta, self.__matrix[2,2])
    else:
      psi = 0.0
      phi = math.atan2(-self.__matrix[1,0], self.__matrix[0,0])
      if math.copysign(1.0, self.__matrix[2,2]) > 0:
        theta = 0.0
      else:
        theta = math.pi

    euler = Euler(phi=phi, theta=theta, psi=psi, unit="rad")
    if unit == "deg":
      euler.degrees()
      return euler
    elif unit == "rad":
      return euler
    else:
      raise ValueError("Unit must be either \"deg\" or \"rad\".")

  def copy(self):
    euler = self.getEuler("rad")
    return RotationMatrix(euler.phi, euler.theta, euler.psi, "rad")

  def fromTransform(self, transformMatrix):
    """Set matrix from string describig rotation matrix in a TRF file."""
    self.__matrix = numpy.array(transformMatrix, numpy.float_).reshape(3,3)

class GridSearch:
  """I3 Local Rotational Grid Search."""
  def __init__(
      self, thetaMax=0., thetaStep=0., psiMax=0., psiStep=0., do180=False):
    self.thetaMax = thetaMax
    self.thetaStep = thetaStep
    self.psiMax = psiMax
    self.psiStep = psiStep
    self.do180 = do180
    self.rotations = (self.thetaMax, self.thetaStep,
                      self.psiMax, self.psiStep,
                      self.do180)
  
  @property
  def thetaMax(self):
    """Maximum half-angle of nutation."""
    return self.__thetaMax

  @thetaMax.setter
  def thetaMax(self, value):
    if value < 0:
      raise ValueError("Theta maximum cannot be less than zero.")
    elif value > 180:
      raise ValueError("Theta maximum cannot be greater than 180 degrees.")
    else:
      self.__thetaMax = float(value)

  @property
  def thetaStep(self):
    """Angular step of nutation."""
    return self.__thetaStep

  @thetaStep.setter
  def thetaStep(self, value):
    if value < 0:
      raise ValueError("Theta step cannot be less than zero.")
    else:
      self.__thetaStep = float(value)
    
  @property
  def psiMax(self):
    """Maximum absolute spin angle (searched from plus/minus)."""
    return self.__psiMax

  @psiMax.setter
  def psiMax(self, value):
    if value < 0:
      raise ValueError("Psi maximum cannot be less than zero.")
    elif value > 180:
      raise ValueError("Psi maximum cannot be greater than 180 degrees.")
    else:
      self.__psiMax = float(value)

  @property
  def psiStep(self):
    """Angular step of spin angle."""
    return self.__psiStep

  @psiStep.setter
  def psiStep(self, value):
    if value < 0:
      raise ValueError("Psi step cannot be less than zero.")
    else:
      self.__psiStep = float(value)
    
  @property
  def do180(self):
    """Determine whether or not to search the flip side of the spin angle."""
    return self.__do180

  @do180.setter
  def do180(self, value):
    if value == True or value == False:
      self.__do180 = value
    else:
      raise ValueError("Do180 must be True or False.")

  @property
  def rotations(self):
    return self.__rotations

  @rotations.setter
  def rotations(self, params):
    thetaMax, thetaStep, psiMax, psiStep, do180 = params
    self.__rotations = [Euler(0., 0., 0.)]
    if do180:
      self.__rotations.append(Euler(0., 0., 180.))

    for psiDeg in numpy.arange(psiStep, psiMax, psiStep):
      self.__rotations.append(Euler(0., 0., psiDeg))
      self.__rotations.append(Euler(0., 0., -psiDeg))
      if do180:
        self.__rotations.append(Euler(0., 0., 180 + psiDeg))
        self.__rotations.append(Euler(0., 0., 180 - psiDeg))

    if psiMax < 180:
      self.__rotations.append(Euler(0., 0., psiMax))
      self.__rotations.append(Euler(0., 0., -psiMax))
      if do180:
        self.__rotations.append(Euler(0., 0., 180 + psiMax))
        self.__rotations.append(Euler(0., 0., 180 - psiMax))

    if thetaStep > 0 and thetaMax > 0:
      thetaDeg = thetaStep
      while thetaDeg <= thetaMax:
        thetaRad = math.radians(thetaDeg)
        phiNum = math.sin(thetaRad) * (180.0 / thetaStep)
        phiNum = 4 if phiNum < 4 else math.ceil(phiNum)
        for i in range(phiNum):
          phiDeg = i * 360.0 / phiNum
          phiRad = math.radians(phiDeg)
          psiRad0 = math.atan2(-math.sin(phiRad),
                               math.cos(thetaRad) * math.cos(phiRad))
          psiDeg0 = math.degrees(psiRad0)

          self.__rotations.append(Euler(phiDeg, thetaDeg, psiDeg0))
          if do180:
            self.__rotations.append(Euler(phiDeg, thetaDeg, psiDeg0+180))

          for psiDeg in numpy.arange(psiStep, psiMax, psiStep):
            self.__rotations.append(Euler(phiDeg, thetaDeg, psiDeg0 + psiDeg))
            self.__rotations.append(Euler(phiDeg, thetaDeg, psiDeg0 - psiDeg))
            if do180:
              self.__rotations.append(
                  Euler(phiDeg, thetaDeg, 180 + psiDeg0 + psiDeg))
              self.__rotations.append(
                  Euler(phiDeg, thetaDeg, 180 + psiDeg0 - psiDeg))
          if psiMax < 180:
            self.__rotations.append(Euler(phiDeg, thetaDeg, psiDeg0 + psiMax))
            self.__rotations.append(Euler(phiDeg, thetaDeg, psiDeg0 - psiMax))
            if do180:
              self.__rotations.append(
                  Euler(phiDeg, thetaDeg, 180 + psiDeg0 + psiMax))
              self.__rotations.append(
                  Euler(phiDeg, thetaDeg, 180 + psiDeg0 - psiMax))
        thetaDeg += thetaStep

  def __len__(self):
    return len(self.rotations)

  def __iter__(self):
    return iter(self.rotations)

  def __getitem__(self, key):
    return self.rotations[key]

  def updateRotations(self):
    self.rotations = (self.thetaMax, self.thetaStep,
                      self.psiMax, self.psiStep,
                      self.do180)

class Transform:
  """A single line in a trf file describing a single raw motifs parameters."""
  def __init__(self, transformLine):
    transformContents = transformLine.split()

    if len(transformContents) < 16:
      raise ValueError("Transform line must have at least 16 fields.")
    elif len(transformContents) == 16:
      self.score = None
      self.classNumber = None
    elif len(transformContents) == 17:
      self.score = transformContents[16]
      self.classNumber = None
    elif len(transformContents) >= 18:
      self.score = transformContents[16]
      self.classNumber = transformContents[17]

    self.subset = transformContents[0]
    self.coordinates = transformContents[1:4]
    self.shifts = transformContents[4:7]
    self.rotation = transformContents[7:16]

  @property
  def subset(self):
    return self.__subset

  @subset.setter
  def subset(self, value):
    self.__subset = value

  @property
  def coordinates(self):
    """An integer column vector of subvolume center relative to tomogram."""
    return self.__coordinates

  @coordinates.setter
  def coordinates(self, coords):
    self.__coordinates = numpy.array(coords, numpy.int_).reshape(3,1)

  @property
  def shifts(self):
    """A float column vector of subvolume displacement from center coords."""
    return self.__shifts

  @shifts.setter
  def shifts(self, displacements):
    self.__shifts = numpy.array(displacements, numpy.float_).reshape(3,1)

  @property
  def rotation(self):
    """A RotationMatrix object describing subvolume coordinate system."""
    return self.__rotation

  @rotation.setter
  def rotation(self, matrix):
    self.__rotation = RotationMatrix()
    self.__rotation.fromTransform(matrix)

  @property
  def score(self):
    """Cross-correlation coefficient."""
    return self.__score

  @score.setter
  def score(self, value):
    self.__score = float(value) if value is not None else None

  @property
  def classNumber(self):
    return self.__classNumber

  @classNumber.setter
  def classNumber(self, value):
    self.__classNumber = int(value) if value is not None else None

  def __str__(self):
    result = self.subset
    result += ' {:d} {:d} {:d} '.format(*self.coordinates[:,0])
    result += '{: 14.10f} {: 14.10f} {: 14.10f} '.format(*self.shifts[:,0])
    result += ('{: 18.16f} ' * 8 + '{: 18.16f}').format(
        *self.rotation.getVector())

    if self.score is not None:
      result += ' {: 14.10f}'.format(self.score)
    elif self.classNumber is not None:
      result += ' {:d}'.format(self.classNumber)
    return result

  def scale(self, scaleFactor):
    particleCenter = self.coordinates + self.shifts
    particleCenter *= scaleFactor
    newShifts, newCoordinates = numpy.modf(particleCenter)
    self.shifts = newShifts[:,0]
    self.coordinates = newCoordinates.astype(numpy.int_)[:,0]

  def addShift(self, x=0.0, y=0.0, z=0.0):
    referenceShift = numpy.array([[x], [y], [z]], numpy.float_)
    particleCenter = self.coordinates + self.shifts
    inverseRotationMatrix = self.rotation.copy()
    inverseRotationMatrix.invert()
    particleShift = inverseRotationMatrix.matrix.dot(referenceShift)
    particleCenter += particleShift
    newShifts, newCoordinates = numpy.modf(particleCenter)
    self.shifts = newShifts[:,0]
    self.coordinates = newCoordinates.astype(numpy.int_)[:,0]

  def copy(self):
    return Transform(str(self))

class TransformList:
  def __init__(self, transforms = []):
    self.filename = ''
    self.transforms = transforms

  @property
  def filename(self):
    return self.__filename

  @filename.setter
  def filename(self, value):
    self.__filename = value

  @property
  def transforms(self):
    return self.__transforms

  @transforms.setter
  def transforms(self, trfList):
    self.__transforms = trfList

  def fromFile(self, filename):
    with open(filename) as transformFile:
      self.transforms = [ Transform(trf) for trf in transformFile ]
    self.filename = filename
  
  def __len__(self):
    return len(self.transforms)

  def __iter__(self):
    return iter(self.transforms)

  def __getitem__(self, key):
    return self.transforms[key]

  def sortByScore(self):
    if self[0].score == None:
      raise ValueError("Transform list does not contain scores.")
    else:
      self.transforms.sort(key=lambda trf: trf.score)

  def sortByClassNumber(self):
    if self[0].classNumber == None:
      raise ValueError("Transform list does not contatin class numbers.")
    else:
      self.transforms.sort(key=lambda trf: trf.classNumber)

  def getBySubset(self, subset):
    trfList = []
    for trf in self:
      if trf.subset == subset:
        trfList.append(trf.copy())
    return TransformList(trfList)

  def getByClassNumber(self, classNumber):
    if self[0].classNumber == None:
      raise ValueError("Transform list does not contatin class numbers.")

    trfList = []
    for trf in self:
      if trf.classNumber == classNumber:
        trfList.append(trf.copy())
    return TransformList(trfList)
    
  def toFile(self, filename):
    with open(filename, 'w') as transformFile:
      for trf in self.transforms:
        transformFile.write(str(trf) + '\n')

  def scale(self, scaleFactor):
    for trf in self:
      trf.scale(scaleFactor)

  def addShift(self, x=0.0, y=0.0, z=0.0):
    for trf in self:
      trf.addShift(x, y, z)

  def copy(self):
    trfList = []
    if len(self) > 0:
      trfList = [ trf.copy() for trf in self ]
    return TransformList(trfList)
