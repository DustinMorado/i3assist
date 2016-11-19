#. -*- coding: utf-8 -*-
"""A module to help using I3 for sub-tomogram averaging.

Written By Dustin Reed Morado
Last updated 04.11.2016
"""
from __future__ import print_function, division
import numpy
import numpy.linalg

class Euler(object):
    """Describes a particle's orientation using ZXZ extrinsic euler angles.

    The first rotation (phi) is about the z-axis (z).
    The second rotation (theta) is about the new x-axis (x').
    The third rotation (psi) is about the new z-axis (z'') after the second
    rotation.

    All rotations are of the coordinate axes and not of point coordinates.

    Args:
        phi (:obj:`float`, optional): First rotation about the z-axis (z).
        theta (:obj:`float`, optional): Second rotation about the new x-axis
            (x').
        psi (:obj:`float`, optional): Third rotation about the new z-axis
            (z'').
        unit (:obj:`str`, optional): Whether rotations are in degrees "deg" or
            radians "rad".

    """
    def __init__(self, phi=0., theta=0., psi=0., unit="deg"):
        self.__unit = None
        self.unit = unit
        self.__angles = None
        self.angles = numpy.array([phi, theta, psi], numpy.float_)

    @property
    def unit(self):
        """:obj:`str`: Describes the unit of rotations as degrees or radians.

        Args:
            value (:obj:`str`): "deg" to set unit as degrees or "rad" to set
                unit as radians.

        """
        return self.__unit

    @unit.setter
    def unit(self, value):
        if value == "deg" or value == "rad":
            self.__unit = value
        else:
            raise ValueError("Unit must be either \"deg\" or \"rad\".")

    @property
    def angles(self):
        """:obj:`list` of :obj:`float`: Array of rotations in order.

        Array is in the order phi, theta, psi

        Args:
            value (:obj:`list` of :obj:`float`): List-type array of the three
                rotations.

        """
        return self.__angles

    @angles.setter
    def angles(self, value):
        temp = numpy.array(value, numpy.float_)
        self.__angles = temp.reshape(3,)

    @property
    def phi(self):
        """:obj:`float`: First rotation about the z-axis (z).

        Args:
            value (:obj:`float`): New phi rotation.

        """
        return self.__angles[0]

    @phi.setter
    def phi(self, value):
        self.__angles[0] = value

    @property
    def theta(self):
        """:obj:`float`: Second rotation about the new x-axis (x').

        Args:
            value (:obj:`float`): New theta rotation.

        """
        return self.__angles[1]

    @theta.setter
    def theta(self, value):
        self.__angles[1] = value

    @property
    def psi(self):
        """:obj:`float`: Third rotation about the new z-axis (z").

        Args:
            value (:obj:`float`): New psi rotation.

        """
        return self.__angles[2]

    @psi.setter
    def psi(self, value):
        self.__angles[2] = value

    def copy(self):
        """Returns a copy of the euler object.

        Returns:
            :obj:`i3assist.Euler`: Copy of Euler object.

        """
        return Euler(phi=self.phi, theta=self.theta, psi=self.psi,
                     unit=self.unit)

    def degrees(self, inplace=False):
        """Converts Euler object to describe rotations in units of degrees.

        Args:
            inplace (:obj:`boolean`, optional): If True conversion will be
                done in-place.

        Returns:
            :obj:`i3assist.Euler`: Copy or self object with rotations in degree
                units.

        """
        if inplace is True:
            if self.unit == "rad":
                self.unit = "deg"
                self.angles = numpy.degrees(self.angles)
            return self
        else:
            result = self.copy()
            if self.unit == "rad":
                result.unit = "deg"
                result.angles = numpy.degrees(result.angles)
            return result

    def radians(self, inplace=False):
        """Converts Euler object to describe rotations in units of radians.

        Args:
            inplace (:obj:`boolean`, optional): If True conversion will be done
                in-place.

        Returns:
            :obj:`i3assist.Euler`: Copy or self object with rotations in radian
                units.

        """
        if inplace is True:
            if self.unit == "deg":
                self.unit = "rad"
                self.angles = numpy.radians(self.angles)
            return self
        else:
            result = self.copy()
            if self.unit == "deg":
                result.unit = "rad"
                result.angles = numpy.radians(result.angles)
            return result

    def __str__(self):
        return "{: f} {: f} {: f}".format(*self.angles)

    def pos_string(self):
        """Returns string of angles in old I3 pos format.

            Returns:
                :obj:`str`: Old I3 pos format string of euler angles.

        """
        return "{: 10.4f}{: 10.4f}{: 10.4f}".format(*self.angles)

    def to_matrix(self):
        """Returns the RotationMatrix object equivalent of Euler object.

        Returns:
            :obj:`i3assist.RotationMatrix`: Rotation matrix describing the
                equivalent euler angles.

        """
        return RotationMatrix(angles=self.angles, unit=self.unit)

    def trf_string(self):
        """Returns string of angles in new I3 trf format.

            Returns:
                :obj:`str`: New I3 trf format string of rotation matrix

        """
        return self.to_matrix().trf_string()

    def normalize(self, inplace=False):
        """ Returns a copy of euler object with rotations in old I3 bounds.

        Old I3 requires that euler angles be given within the bounds:
        phi: -180 to 180
        theta: 0 to 180
        psi: -180 to 180

        There are many equivalent euler angle formats and we can normalize the
        angles using these identities, but these ranges correspond to the ranges
        given by the default C math arc trig functions, which are used in going
        from rotation matrix to Euler angles so we just convert to rotation
        matrix and then back to Euler angles.

        Args:
            inplace: If set to True the Euler angles will be changed in-place
                instead of returning a copy.

        Returns:
            :obj:`i3assist.Euler`: Copy or self object with rotations
                normalized to ranges required by old I3.

        """
        result = self.to_matrix().to_euler(self.unit)
        if inplace is True:
            self.angles = result.angles
            del result
            return self
        else:
            return result

    def invert(self, inplace=False):
        """ Returns a copy of euler object with inverted rotations.

        Again we can simply take the negative of psi, theta, phi to invert the
        rotatations but it is more theoretically clear to convert the angles to
        a rotation matrix, transpose it, and convert it back to angles. We also
        get the added benefit of normalizing the angles.

        Args:
            inplace: If set to True the Euler angles will be changed in-place
                instead of returning a copy.

        Returns:
            :obj:`i3assist.Euler`: Copy or self object describing the opposite
                rotation.

        """
        result = self.to_matrix().invert().to_euler(self.unit)
        if inplace is True:
            self.angles = result.angles
            del result
            return self
        else:
            return result

    def __add__(self, augend):
        if not isinstance(augend, Euler):
            raise TypeError("Augend must also be an Euler object instance.")

        addend_matrix = self.to_matrix()
        augend_matrix = augend.to_matrix()
        sum_matrix = augend_matrix + addend_matrix
        return sum_matrix.to_euler(self.unit)

class RotationMatrix(object):
    """Describes a particle's orientation using ZXZ passive rotation matrix.

    The rotation matrix is the composition of three rotations matrices with
    the first and third being about the Z-axis and the second about the X-axis.

    All rotations are passive (alias) transformations of the coordinate axes
    and not of point coordinates (active / alibi).

    Args:
        angles (:obj:`list` of :obj:`float`, optional): Euler angles describing
            the rotation.
        unit (:obj:`str`, optional): Whether rotations are in degrees "deg" or
            radians "rad"

    """
    def __init__(self, angles=None, unit="deg"):
        self.__matrix = None
        if angles is None:
            angles = numpy.array([0., 0., 0.], numpy.float_)

        angles_array = numpy.array(angles, numpy.float_)
        if unit == "deg" and angles_array.size == 3:
            angles_array = numpy.radians(angles)
        elif unit == "rad" or unit == "deg":
            pass
        else:
            raise ValueError("Unit must be either \"deg\" or \"rad\".")

        self.matrix = angles_array

    @property
    def matrix(self):
        """:obj:`numpy.ndarray`: (3,3) matrix describing the rotation.

        Matrix is ordered in the standard format with the first index
        describing the row and the second index describing the column.

        Args:
            value (:obj:`numpy.ndarray`): Euler angles in radians or Rotation
                matrix.

        """
        return self.__matrix

    @matrix.setter
    def matrix(self, value):
        if value.size == 3:
            cos = numpy.cos(value)
            sin = numpy.sin(value)
            matrix = numpy.identity(3, numpy.float_)
            matrix[0, 0] = cos[2] * cos[0] - sin[2] * cos[1] * sin[0]
            matrix[0, 1] = cos[2] * sin[0] + sin[2] * cos[1] * cos[0]
            matrix[0, 2] = sin[2] * sin[1]
            matrix[1, 0] = -sin[2] * cos[0] - cos[2] * cos[1] * sin[0]
            matrix[1, 1] = -sin[2] * sin[0] + cos[2] * cos[1] * cos[0]
            matrix[1, 2] = cos[2] * sin[1]
            matrix[2, 0] = sin[1] * sin[0]
            matrix[2, 1] = -sin[1] * cos[0]
            matrix[2, 2] = cos[1]
            self.__matrix = matrix
        elif value.size == 9:
            matrix = numpy.array(value, numpy.float_).reshape(3, 3)
            if not numpy.allclose(
                    matrix.dot(matrix.T), numpy.identity(3, numpy.float_)):
                raise ValueError("Given matrix is not orthonormal.")
            elif not numpy.isclose(numpy.linalg.det(matrix), 1.0):
                raise ValueError("Given matrix is not a proper rotation.")
            else:
                self.__matrix = matrix
        else:
            raise TypeError("Input must be array of angles or rotation matrix")

    def __str__(self):
        return ("{: f} " * 8 + "{: f}").format(*self.matrix.reshape(-1))

    def invert(self, inplace=False):
        """Returns rotation matrix describing opposite rotation.

        Rotation matrices are orthogonal and hence their transpose is their
        inverse.

        Args:
            inplace (:obj:`boolean`, optional): If True the operation will be
                done on the object in-place instead of returning a copy.

        Returns:
            :obj:`i3assist.RotationMatrix`: Copy or self object of the inverted
                rotation matrix.

        """
        if inplace is True:
            self.matrix = self.matrix.T
            return self
        else:
            result = RotationMatrix()
            result.matrix = self.matrix.T
            return result

    def transpose(self, inplace=False):
        """Returns rotation matrix describing opposite rotation.

        Rotation matrices are orthogonal and hence their transpose is their
        inverse.

        Args:
            inplace (:obj:`boolean`, optional): If True the operation will be
                done on the object in-place instead of returning a copy.

        Returns:
            :obj:`i3assist.RotationMatrix`: Copy or self object of the inverted
                rotation matrix.

        """
        if inplace is True:
            self.matrix = self.matrix.T
            return self
        else:
            result = RotationMatrix()
            result.matrix = self.matrix.T
            return result

    def to_euler(self, unit="deg"):
        """Converts rotation matrix to equivalent euler angles.

        Algorithm from Chapter 1 of "Computational Methods for
        Three-Dimensional Microscopy Reconstruction" Ed. Joachim Frank, Gabor
        Herman.

        Args:
            unit (:obj:`str`, optional): Whether to return the Euler angles in
                degrees or radians.

        Returns:
            :obj:`i3assist.Euler`: Euler angle equivalent of rotation matrix

        """
        abs_sine_theta = self.matrix[0, 2] ** 2 + self.matrix[1, 2] ** 2
        abs_sine_theta = numpy.sqrt(abs_sine_theta)

        if not numpy.isclose(abs_sine_theta, 0.0):
            phi = numpy.arctan2(self.matrix[2, 0], -self.matrix[2, 1])
            psi = numpy.arctan2(self.matrix[0, 2], self.matrix[1, 2])

            if numpy.isclose(numpy.sin(phi), 0.0):
                theta_sign = self.matrix[1, 2] / numpy.cos(psi)
                theta_sign = numpy.copysign(1.0, theta_sign)
            else:
                theta_sign = self.matrix[0, 2] / numpy.sin(psi)
                theta_sign = numpy.copysign(1.0, theta_sign)

            theta = numpy.arctan2(
                theta_sign * abs_sine_theta, self.matrix[2, 2])

        else:
            psi = 0.0
            phi = numpy.arctan2(-self.matrix[1, 0], self.matrix[0, 0])

            if numpy.copysign(1.0, self.matrix[2, 2]) > 0:
                theta = 0.0
            else:
                theta = numpy.pi

        euler = Euler(phi=phi, theta=theta, psi=psi, unit="rad")

        if unit == "deg":
            return euler.degrees()
        elif unit == "rad":
            return euler
        else:
            raise ValueError("Unit must be either \"deg\" or \"rad\".")

    def copy(self):
        """Returns a copy of the rotation matrix."""
        return RotationMatrix(angles=self.matrix, unit="rad")

    def pos_string(self):
        """Returns string of angles in old I3 pos format.

            Returns:
                :obj:`str`: Old I3 pos format string of euler angles.

        """
        euler = self.to_euler()
        return euler.pos_string()

    def trf_string(self):
        """Returns string of angles in new I3 trf format.

            Returns:
                :obj:`str`: New I3 trf format string of rotation matrix

        """
        matrix_array = self.matrix.reshape(-1)
        return ("{: 19.15f}" * 9).format(*matrix_array)

    def __add__(self, augend):
        sum_matrix = augend.matrix.dot(self.matrix)
        result = RotationMatrix()
        result.matrix = sum_matrix
        return result

class GridSearch(object):
    """Describes the local grid search implemented in new I3.

    See the explanation in MRASRCH and the I3 subvolume tutorial for more
    information on how the grid search is implemented. But overall the grid is
    defined by four parameters: Nutation (theta) maximum and step and Spin
    (psi) maximum and step. Finally there is a do_180 parameter to support the
    old I3 eulerFG scripts, but this is not available in new I3.

    Args:
        theta_max (:obj:`float`, optional): Maximum half-angle of a cone of
            nutation about the north pole of the unit sphere.
        theta_step (:obj:`float`, optional): Angular increment of nutation.
        psi_max (:obj:`float`, optional): Maximum absolute angle of spin about
            the orientation axis of the particle. Searched in both directions.
        psi_step (:obj:`float`, optional): Angular increment of spin.
        do_180 (:obj:`boolean`, optional): If True the spins opposite of the
            current orientation's facing will also be searched.

    """
    def __init__(
            self, theta_max=0.0, theta_step=0.0, psi_max=0.0, psi_step=0.0,
            do_180=False):
        self.__theta_max = None
        self.__theta_step = None
        self.__psi_max = None
        self.__psi_step = None
        self.__do_180 = None
        self.__rotations = None

        self.theta_max = theta_max
        self.theta_step = theta_step
        self.psi_max = psi_max
        self.psi_step = psi_step
        self.do_180 = do_180
        self.rotations = (
            self.theta_max, self.theta_step, self.psi_max, self.psi_step,
            self.do_180)

    @property
    def rotations(self):
        """:obj:`list` of :obj:`i3assist.Euler` Rotations searched in i3.

        Args:
            params (:obj:`tuple` of :obj:`float` and :obj:`boolean`): A tuple
                with theta_max, theta_step, psi_max, psi_step, and do_180.

        """
        return self.__rotations

    @rotations.setter
    def rotations(self, params):
        theta_max, theta_step, psi_max, psi_step, do_180 = params
        self.__rotations = [Euler(0., 0., 0.)]
        if do_180:
            self.__rotations.append(Euler(0., 0., 180.))

        if psi_step > 0. and psi_max > 0.:
            for psi_deg in numpy.arange(psi_step, psi_max, psi_step):
                self.__rotations.append(Euler(0., 0., psi_deg))
                self.__rotations.append(Euler(0., 0., -psi_deg))
                if do_180:
                    self.__rotations.append(Euler(0., 0., 180 + psi_deg))
                    self.__rotations.append(Euler(0., 0., 180 - psi_deg))

            if psi_max < 180:
                self.__rotations.append(Euler(0., 0., psi_max))
                self.__rotations.append(Euler(0., 0., -psi_max))
                if do_180:
                    self.__rotations.append(Euler(0., 0., 180 + psi_max))
                    self.__rotations.append(Euler(0., 0., 180 - psi_max))

        if theta_step > 0 and theta_max > 0:
            theta_deg = theta_step
            while theta_deg <= theta_max:
                theta_rad = numpy.radians(theta_deg)
                phi_num = numpy.sin(theta_rad) * (180.0 / theta_step)
                phi_num = 4 if phi_num < 4 else int(numpy.ceil(phi_num))

                for i in range(phi_num):
                    phi_deg = i * 360.0 / phi_num
                    phi_rad = numpy.radians(phi_deg)
                    psi_rad0 = numpy.arctan2(
                        -numpy.sin(phi_rad),
                        numpy.cos(theta_rad) * numpy.cos(phi_rad))
                    psi_deg0 = numpy.degrees(psi_rad0)

                    self.__rotations.append(
                        Euler(phi_deg, theta_deg, psi_deg0))
                    if do_180:
                        self.__rotations.append(
                            Euler(phi_deg, theta_deg, psi_deg0+180))

                    if psi_step > 0. and psi_max > 0.:
                        for psi_deg in numpy.arange(
                                psi_step, psi_max, psi_step):
                            self.__rotations.append(
                                Euler(phi_deg, theta_deg, psi_deg0 + psi_deg))
                            self.__rotations.append(
                                Euler(phi_deg, theta_deg, psi_deg0 - psi_deg))
                            if do_180:
                                self.__rotations.append(
                                    Euler(phi_deg, theta_deg,
                                          180 + psi_deg0 + psi_deg))
                                self.__rotations.append(
                                    Euler(phi_deg, theta_deg,
                                          180 + psi_deg0 - psi_deg))
                        if psi_max < 180:
                            self.__rotations.append(
                                Euler(phi_deg, theta_deg, psi_deg0 + psi_max))
                            self.__rotations.append(
                                Euler(phi_deg, theta_deg, psi_deg0 - psi_max))
                            if do_180:
                                self.__rotations.append(
                                    Euler(phi_deg, theta_deg,
                                          180 + psi_deg0 + psi_max))
                                self.__rotations.append(
                                    Euler(phi_deg, theta_deg,
                                          180 + psi_deg0 - psi_max))
                theta_deg += theta_step

    @property
    def theta_max(self):
        """:obj:`float` Max half-angle of nutation about the north pole.

        Args:
            value (:obj:`float`): half-angle in the range 0 to 180.

        """
        return self.__theta_max

    @theta_max.setter
    def theta_max(self, value):
        if value < 0:
            raise ValueError("Theta maximum cannot be less than zero.")
        elif value > 180:
            raise ValueError("Theta maximum cannot be greater than 180.")
        else:
            self.__theta_max = float(value)
            if self.rotations is not None:
                self.rotations = (
                    self.theta_max, self.theta_step, self.psi_max,
                    self.psi_step, self.do_180)

    @property
    def theta_step(self):
        """:obj:`float` Angular increment of nutation.

        Args:
            value (:obj:`float`): Angular increment in the range 0 to 180.

        """
        return self.__theta_step

    @theta_step.setter
    def theta_step(self, value):
        if value < 0:
            raise ValueError("Theta step cannot be less than zero.")
        else:
            self.__theta_step = float(value)
            if self.rotations is not None:
                self.rotations = (
                    self.theta_max, self.theta_step, self.psi_max,
                    self.psi_step, self.do_180)

    @property
    def psi_max(self):
        """:obj:`float` Max absolute angle of spin about the particle z-axis.

        Args:
            value (:obj:`float`): Max absolute spin angle in range 0 to 180.

        """
        return self.__psi_max

    @psi_max.setter
    def psi_max(self, value):
        if value < 0:
            raise ValueError("Psi maximum cannot be less than zero.")
        elif value > 180:
            raise ValueError("Psi maximum cannot be greater than 180 degrees.")
        else:
            self.__psi_max = float(value)
            if self.rotations is not None:
                self.rotations = (
                    self.theta_max, self.theta_step, self.psi_max,
                    self.psi_step, self.do_180)

    @property
    def psi_step(self):
        """:obj:`float` Angular increment of spin.

        Args:
            value (:obj:`float`): Spin angular increment in the range 0 to 180.

        """
        return self.__psi_step

    @psi_step.setter
    def psi_step(self, value):
        if value < 0:
            raise ValueError("Psi step cannot be less than zero.")
        else:
            self.__psi_step = float(value)
            if self.rotations is not None:
                self.rotations = (
                    self.theta_max, self.theta_step, self.psi_max,
                    self.psi_step, self.do_180)

    @property
    def do_180(self):
        """:obj:`boolean` Whether to search spins opposite particle facing.

        Args:
            value (:obj:`bool`): True to search opposite facing spin angles.

        """
        return self.__do_180

    @do_180.setter
    def do_180(self, value):
        self.__do_180 = bool(value)
        if self.rotations is not None:
            self.rotations = (
                self.theta_max, self.theta_step, self.psi_max, self.psi_step,
                self.do_180)

    def __len__(self):
        return len(self.rotations)

    def __iter__(self):
        return iter(self.rotations)

    def __getitem__(self, key):
        return self.rotations[key]

class Transform(object):
    """A single particle transform in new I3.

    For more information refer to the subvolume tutorial document for I3.

    Args:
        transformLine (:obj:`str`): A string with the transform data.

    """
    def __init__(self, transform_line):
        transform_contents = transform_line.split()
        self.__score = None
        self.__class_number = None
        self.__subset = None
        self.__coordinates = None
        self.__shifts = None
        self.__rotation = None

        if len(transform_contents) < 16:
            raise ValueError("Transform line must have at least 16 fields.")
        elif len(transform_contents) == 16:
            self.score = None
            self.class_number = None
        elif len(transform_contents) == 17:
            self.score = transform_contents[16]
            self.class_number = None
        elif len(transform_contents) >= 18:
            self.score = transform_contents[16]
            self.class_number = transform_contents[17]

        self.subset = transform_contents[0]
        self.coordinates = [int(x) for x in transform_contents[1:4]]
        self.shifts = [float(x) for x in transform_contents[4:7]]
        self.rotation = [float(x) for x in transform_contents[7:16]]

    @property
    def subset(self):
        """:obj:`str` Subset identifier for particle.

        Args:
            value (:obj:`str`): Subset name around 10 characters.

        """
        return self.__subset

    @subset.setter
    def subset(self, value):
        self.__subset = str(value)

    @property
    def coordinates(self):
        """:obj:`list` of :obj:`int` Particles integer coordinates in tomogram.

        Coordinates are stored here as column arrays to help with using them
        with the rotation matrices.

        Args:
            value (:obj:`list` of :obj:`int`): List with 3 elements for
                particles coordinates relative to the tomogram map it's
                extracted from.

        """
        return self.__coordinates

    @coordinates.setter
    def coordinates(self, value):
        self.__coordinates = numpy.array(value, numpy.int_).reshape(3, 1)

    @property
    def shifts(self):
        """:obj:`list` of :obj:`float` Particle displacements from coordinates.

        Displacements are stored here as column arrays to help with using them
        with the rotation matrices. Again as for rotations translations are
        alibi translations of the coordinate system and not the actual points.

        Args:
            value (:obj:`list' of :obj`float`): List with 3 elements for
                particles shifts relative to the center of the reference.

        """
        return self.__shifts

    @shifts.setter
    def shifts(self, value):
        self.__shifts = numpy.array(value, numpy.float_).reshape(3, 1)

    @property
    def rotation(self):
        """:obj:`i3assist.RotationMatrix` Particles rotation matrix to orient.

        Args:
            value (:obj:`list` of :obj:`float`): Rotation matrix.

        """
        return self.__rotation

    @rotation.setter
    def rotation(self, value):
        self.__rotation = RotationMatrix(angles=value, unit="rad")

    @property
    def score(self):
        """:obj:`float` Correlation score of particle alignment to reference.

        Args:
            value (:obj:`float`): Correlation score.

        """
        return self.__score

    @score.setter
    def score(self, value):
        self.__score = float(value) if value is not None else None

    @property
    def class_number(self):
        """:obj:`int` Class number that the particle belongs to.

        Args:
            value (:obj:`int`): Particle Class.

        """
        return self.__class_number

    @class_number.setter
    def class_number(self, value):
        self.__class_number = int(value) if value is not None else None

    def __str__(self):
        result = self.subset
        result += ' {:d} {:d} {:d} '.format(*self.coordinates[:, 0])
        result += '{: 19.14f} {: 19.14f} {: 19.14f} '.format(
            *self.shifts[:, 0])
        result += self.rotation.trf_string()

        if self.score is not None:
            result += ' {: 14.10f}'.format(self.score)
        elif self.class_number is not None:
            result += ' {:d}'.format(self.class_number)
        return result

    def to_pos(self):
        """Return the equivalent old I3 position.

        Returns:
            :obj:`i3assist.Position`: Old I3 equivalent orientation.

        """
        result = '{:d} {:d} {:d} 0 '.format(*self.coordinates.reshape(-1))
        inv_shifts = self.shifts * -1
        result += '{:10.4f} {:10.4f} {:10.4f} '.format(*inv_shifts.reshape(-1))
        inv_euler = self.rotation.invert().to_euler()
        result += inv_euler.pos_string()
        if self.class_number is not None:
            result += ' {:d} '.format(self.class_number)
        else:
            result += ' 0 '

        if self.score is not None:
            result += ' {:10.5f} {:10.5f}'.format(self.score, 0.0)
        else:
            result += ' {:10.5f} {:10.5f}'.format(0.0, 0.0)

        return Position(result)

    def copy(self):
        """Returns a copy of the transform."""
        return Transform(str(self))

    def scale(self, scale_factor, inplace=False):
        """Scale the transform to handle binning.

        Args:
            scale_factor (:obj:`float`): Amount to scale transform by.
            inplace: If True the operation will be done in-place and modify
                the transform instead of returning a copy.

        Returns:
            :obj:`i3assist.Transform`: Copy or self scaled transform object.

        """
        particle_center = (self.coordinates + self.shifts) * scale_factor
        new_shifts, new_coordinates = numpy.modf(particle_center)
        if inplace is True:
            self.shifts = new_shifts
            self.coordinates = new_coordinates.astype(numpy.int_)
            return self
        else:
            result = self.copy()
            result.shifts = new_shifts
            result.coordinates = new_coordinates.astype(numpy.int_)
            return result

    def add_shift(self, shift_x=0.0, shift_y=0.0, shift_z=0.0, inplace=False):
        """Adjusts the particles defined center by an arbitrary vector.

        Args:
            shift_x (:obj:`float`, optional): Amount to shift in x.
            shift_y (:obj:`float`, optional): Amount to shift in y.
            shift_z (:obj:`float`, optional): Amount to shift in z.
            inplace (:obj:`float`, optional): If True the operation will be
                done in-place and modify thet transform instead of returning a
                copy.

        Returns:
            :obj:`i3assist.Transform`: Copy or self shifted transform object.

        """
        reference_shift = numpy.array(
            [shift_x, shift_y, shift_z], numpy.float_).reshape(3, 1)
        particle_center = self.coordinates + self.shifts
        transform_matrix = self.rotation.invert().matrix
        transformed_shift = transform_matrix.dot(reference_shift)
        particle_center += transformed_shift
        new_shifts, new_coordinates = numpy.modf(particle_center)

        if inplace is True:
            self.shifts = new_shifts
            self.coordinates = new_coordinates.astype(numpy.int_)
            return self
        else:
            result = self.copy()
            result.shifts = new_shifts
            result.coordinates = new_coordinates.astype(numpy.int_)
            return result

    def add_rotation(self, rotation, inplace=False):
        """Adds a rotation in addition to particles current orientation.

        Args:
            rotation (:obj:`i3assist.RotationMatrix`): Rotation to add.
            inplace: If True the operation will be done in-place and modify
                the transform instead of returning a copy.

        Returns:
            :obj:`i3assist.Transform`: Copy or self rotated transform object.

        """
        if not isinstance(rotation, RotationMatrix):
            raise ValueError("Rotation must be RotationMatrix type.")
        new_rotation = self.rotation + rotation
        if inplace is True:
            self.rotation = new_rotation.matrix
            return self
        else:
            result = self.copy()
            result.rotation = new_rotation.matrix
            return result

class TransformList(object):
    """Describes a full transform file with as a list of Transforms.

    For more information refer to the subvolume tutorial document for I3.

    Args:
        filename (:obj:`str`): Filename of trf file.
        transforms (:obj:`list` of :obj:`i3assist.Transform`): List of
            transforms.

    """
    def __init__(self, filename='', transforms=None):
        self.__filename = None
        self.filename = filename
        self.__transforms = None
        self.transforms = transforms

    @property
    def filename(self):
        """Filename associated with transform list.

        Args:
            value (:obj:`str`): Filename of trf file.

        """
        return self.__filename

    @filename.setter
    def filename(self, value):
        self.__filename = str(value)

    @property
    def transforms(self):
        """List of Transform objects.

        Args:
            value (:obj:`list` of :obj:`i3assist.Transform`): List of
                transforms.

        """
        return self.__transforms

    @transforms.setter
    def transforms(self, value):
        self.__transforms = value

    def from_file(self, filename):
        """Loads list of transforms from a trf file.

        Args:
            filename (:obj:`str`): Filename of trf file.

        """
        with open(filename) as transform_file:
            self.transforms = [Transform(trf) for trf in transform_file]
        self.filename = filename

    def __len__(self):
        return len(self.transforms)

    def __iter__(self):
        return iter(self.transforms)

    def __getitem__(self, key):
        return self.transforms[key]

    def sort_by_score(self, inplace=False):
        """Sorts a transform list by correlation coefficient.

        Args:
            inplace: If True the list is not returned and the sorting is done
                in-place. Otherwise a sorted copy is returned.

        Returns:
            :obj:`i3assist.TransformList`: Transform list sorted by score.

        """
        if self[0].score is None:
            raise ValueError("Transform list does not contain scores.")
        else:
            if inplace is True:
                self.transforms.sort(key=lambda trf: trf.score)
            else:
                sorted_transforms = sorted(self.transforms,
                                           key=lambda trf: trf.score)
                return TransformList(filename=self.filename,
                                     transforms=sorted_transforms)

    def sort_by_class(self, inplace=False):
        """Sorts a transform list by class numbers.

        Args:
            inplace: If True the list is not returned and the sorting is done
                in-place. Otherwise a sorted copy is returned.

        Returns:
            :obj:`i3assist.TransformList`: Transform list sorted by class.

        """
        if self[0].classNumber is None:
            raise ValueError("Transform list does not contatin class numbers.")
        else:
            if inplace is True:
                self.transforms.sort(key=lambda trf: trf.class_number)
            else:
                sorted_transforms = sorted(self.transforms,
                                           key=lambda trf: trf.class_number)
                return TransformList(filename=self.filename,
                                     transforms=sorted_transforms)


    def get_by_subset(self, subset):
        """Gets a subset of a transform list based on the subset field.

        Args:
            subset (:obj:`str`): The subset field to search for.

        Returns:
            :obj:`i3assist.TransformList`: Subset of self that matches the.
                subset requested.

        """
        trf_list = []
        for trf in self:
            if trf.subset == subset:
                trf_list.append(trf.copy())
        return TransformList(filename=self.filename,
                             transforms=trf_list)

    def get_by_class(self, class_number):
        """Gets a subset of a transform list based on the class number.

        Args:
            class_number (:obj:`int`): The class number to search for.

        Returns:
            :obj:`i3assist.TransformList`: Subset of self that matches the.
                class number requested.

        """
        if self[0].class_number is None:
            raise ValueError("Transform list does not contatin class numbers.")

        trf_list = []
        for trf in self:
            if trf.class_number == class_number:
                trf_list.append(trf.copy())
        return TransformList(filename=self.filename,
                             transforms=trf_list)

    def to_file(self, filename):
        """Writes out a Transform list to a trf file.

        Args:
            filename (:obj:`str`): Name of file to write to.

        """
        with open(filename, 'w') as transform_file:
            for trf in self.transforms:
                transform_file.write(str(trf) + '\n')

    def scale(self, scale_factor, inplace=False):
        """Scales all of the transforms in a list.

        Args:
            scale_factor (:obj:`float`): Amount which to scale the transforms.
            inplace: If True, scaling will be done in-place instead of
                returning a copy.

        Returns:
            :obj:`i3assist.TransformList`: Scaled transform list.

        """
        if inplace is True:
            for trf in self:
                trf.scale(scale_factor)
        else:
            trf_list = [trf.copy().scale(scale_factor) for trf in self]
            return TransformList(filename=self.filename,
                                 transforms=trf_list)

class Position(object):
    """A single particle transform in old I3.

    The fields in the position are broken down as follows:
        1. Particle's X coordinate in the tomogram *not used in alignment*
        2. Particle's Y coordinate in the tomogram *not used in alignment*
        3. Particle's Z coordinate in the tomogram *not used in alignment*
        4. Particle's 4D index in the stack of particles
        5. Particle's displacement (shift) along X from center of volume.
        6. Particle's displacement (shift) along Y from center of volume.
        7. Particle's displacement (shift) along Z from center of volume.
        8. First euler angle in degrees about Z axis in range [-180, 180].
        9. Second euler angle in degrees about X axis in range [0, 180].
        10. Third euler angle in degrees about Z axis in range [-180, 180].
        11. Class number the particle belongs to.
        12. Max correlation coefficient between the particle and reference
        13. Unknown reserved value.

    Args:
        positionLine (:obj:`str`): A string with the transform data.

    """
    def __init__(self, position_line):
        position_contents = position_line.split()
        self.__score = None
        self.__class_number = None
        self.__coordinates = None
        self.__stack_index = None
        self.__shifts = None
        self.__rotation = None

        if len(position_contents) < 13:
            raise ValueError("Position line must have at least 13 fields.")

        self.coordinates = [int(x) for x in position_contents[0:3]]
        self.stack_index = position_contents[3]
        self.shifts = [float(x) for x in position_contents[4:7]]
        self.rotation = [float(x) for x in position_contents[7:10]]
        self.class_number = int(position_contents[10])
        self.score = position_contents[11]

    @property
    def coordinates(self):
        """:obj:`list` of :obj:`int` Particles integer coordinates in tomogram.

        Coordinates are stored here as column arrays to help with using them
        with the rotation matrices.

        Args:
            value (:obj:`list` of :obj:`int`): List with 3 elements for
                particles coordinates relative to the raw tomogram it's
                extracted from.

        """
        return self.__coordinates

    @coordinates.setter
    def coordinates(self, value):
        self.__coordinates = numpy.array(value, numpy.int_).reshape(3, 1)

    @property
    def stack_index(self):
        """:obj:`int` Fourth dimension coordinate of particle in stack.

        Old I3 used 4D stacks of particles as input. To keep the ability to
        transform coordinates with rotation matrices and shifts, the fourth
        dimension coordinate is given its own attribute.

        Args:
            value (:obj:`int`): Particles' fourth dimension coordinate

        """
        return self.__stack_index

    @stack_index.setter
    def stack_index(self, value):
        self.__stack_index = int(value)

    @property
    def shifts(self):
        """:obj:`list` of :obj:`float` Particle displacements from box center.

        Displacements are stored here as column arrays to help with using them
        with the rotation matrices. Again as for rotations, translations are
        alibi translations of the coordinate system and not the actual points.

        Args:
            value (:obj:`list` of :obj:`float`): List with 3 elements for
                particle's shifts relative to their box center.

        """
        return self.__shifts

    @shifts.setter
    def shifts(self, value):
        self.__shifts = numpy.array(value, numpy.float_).reshape(3, 1)

    @property
    def rotation(self):
        """:obj:`i3assist.Euler` Euler angles rotating reference to particle.

        In Old I3 the rotations are opposite the ones given in New I3. This
        could mean that previously alias rotations were used, but I think it is
        more likely that the rotations just describe the rotation of the
        reference.

        Args:
            value (:obj:`list` of :obj:`float`): Euler angles phi, theta, psi

        """
        return self.__rotation

    @rotation.setter
    def rotation(self, value):
        self.__rotation = Euler(phi=value[0], theta=value[1], psi=value[2])

    @property
    def score(self):
        """:obj:`float` Correlation score of particle alignment to reference.

        Args:
            value (:obj:`float`): Correlation score.

        """
        return self.__score

    @score.setter
    def score(self, value):
        self.__score = float(value)

    @property
    def class_number(self):
        """:obj:`int` Class number that the particle belongs to.

        Args:
            value (:obj:`int`): Particle Class.

        """
        return self.__class_number

    @class_number.setter
    def class_number(self, value):
        self.__class_number = int(value)

    def __str__(self):
        result = '{:d} {:d} {:d} '.format(*self.coordinates.reshape(-1))
        result += '{:d} '.format(self.stack_index)
        result += '{:10.4f} {:10.4f} {:10.4f} '.format(
            *self.shifts.reshape(-1))
        result += self.rotation.pos_string()
        result += ' {:d} '.format(self.class_number)
        result += '{:10.5f} {:10.5f}'.format(self.score, 0.0)
        return result

    def copy(self):
        """Returns a copy of the position."""
        return Position(str(self))

    def to_trf(self):
        """Returns the equivalent new I3 transform.

        Returns:
           :obj:`i3assist.Transform`: New I3 equivalent orientation.

        """
        result = 'SUBSET '
        result += '{:d} {:d} {:d} '.format(*self.coordinates.reshape(-1))
        inv_shifts = self.shifts * -1
        result += '{: 19.14f} {: 19.14f} {: 19.14f} '.format(
            *inv_shifts.reshape(-1))
        inv_rotmat = self.rotation.to_matrix().invert()
        result += inv_rotmat.trf_string()
        result += ' {: 14.10f} '.format(self.score)
        result += ' {:d}'.format(self.class_number)
        return Transform(result)

class PositionList(object):
    """Describes a full position file as a list of Positions.

    Refer to the documentation for Position objects

    Args:
        filename (:obj:`str`): Filename of pos file.
        positions (:obj:`list` of :obj:`i3assist.Position`): List of
            positions.

    """
    def __init__(self, filename='', positions=None):
        self.__filename = None
        self.filename = filename
        self.__positions = None
        self.positions = positions

    @property
    def filename(self):
        """Filename associated with position list.

        Args:
            value (:obj:`str`): Filename of pos file.

        """
        return self.__filename

    @filename.setter
    def filename(self, value):
        self.__filename = str(value)

    @property
    def positions(self):
        """List of Position objects.

        Args:
            value (:obj:`list` of :obj:`i3assist.Position`): List of
                positions.

        """
        return self.__positions

    @positions.setter
    def positions(self, value):
        self.__positions = value

    def from_file(self, filename):
        """Loads list of positions from a trf file.

        Args:
            filename (:obj:`str`): Filename of pos file.

        """
        with open(filename) as position_file:
            self.positions = [Position(pos) for pos in position_file]
        self.filename = filename

    def __len__(self):
        return len(self.positions)

    def __iter__(self):
        return iter(self.positions)

    def __getitem__(self, key):
        return self.positions[key]

    def to_file(self, filename):
        """Writes out a Position list to a pos file.

        Args:
            filename (:obj:`str`): Name of file to write to.

        """
        with open(filename, 'w') as position_file:
            for pos in self.positions:
                position_file.write(str(pos) + '\n')

# vim: textwidth=80 tabstop=8 expandtab shiftwidth=4 softtabstop=4:
