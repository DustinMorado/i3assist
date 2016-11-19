import unittest
import i3assist
import numpy

class TestEulerMethods(unittest.TestCase):

    def test_default_deg_initialization_unit(self):
        default_euler = i3assist.Euler()
        self.assertEqual(default_euler.unit, "deg")

    def test_default_deg_initialization_phi(self):
        default_euler = i3assist.Euler()
        self.assertAlmostEqual(default_euler.phi, 0.)

    def test_default_deg_initialization_theta(self):
        default_euler = i3assist.Euler()
        self.assertAlmostEqual(default_euler.theta, 0.)

    def test_default_deg_initialization_psi(self):
        default_euler = i3assist.Euler()
        self.assertAlmostEqual(default_euler.psi, 0.)

    def test_default_rad_initialization_unit(self):
        default_euler = i3assist.Euler(unit="rad")
        self.assertEqual(default_euler.unit, "rad")

    def test_default_rad_initialization_phi(self):
        default_euler = i3assist.Euler(unit="rad")
        self.assertAlmostEqual(default_euler.phi, 0.)

    def test_default_rad_initialization_theta(self):
        default_euler = i3assist.Euler(unit="rad")
        self.assertAlmostEqual(default_euler.theta, 0.)

    def test_default_rad_initialization_psi(self):
        default_euler = i3assist.Euler(unit="rad")
        self.assertAlmostEqual(default_euler.psi, 0.)

    def test_custom_deg_initialization(self):
        custom_euler = i3assist.Euler(phi=56.0, theta=53.0, psi=-304.0,
                                      unit="deg")
        self.assertEqual(custom_euler.unit, "deg")
        self.assertAlmostEqual(custom_euler.phi, 56.0)
        self.assertAlmostEqual(custom_euler.theta, 53.0)
        self.assertAlmostEqual(custom_euler.psi, -304.0)

    def test_custom_rad_initialization(self):
        custom_euler = i3assist.Euler(phi=5.27, theta=1.23, psi=-4.64,
                                      unit="rad")
        self.assertEqual(custom_euler.unit, "rad")
        self.assertAlmostEqual(custom_euler.phi, 5.27)
        self.assertAlmostEqual(custom_euler.theta, 1.23)
        self.assertAlmostEqual(custom_euler.psi, -4.64)

    def test_default_deg_copy(self):
        initial_euler = i3assist.Euler()
        copy_euler = initial_euler.copy()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(initial_euler.unit, copy_euler.unit)
        self.assertAlmostEqual(initial_euler.phi, copy_euler.phi)
        self.assertAlmostEqual(initial_euler.theta, copy_euler.theta)
        self.assertAlmostEqual(initial_euler.psi, copy_euler.psi)

    def test_default_rad_copy(self):
        initial_euler = i3assist.Euler(unit="rad")
        copy_euler = initial_euler.copy()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(initial_euler.unit, copy_euler.unit)
        self.assertAlmostEqual(initial_euler.phi, copy_euler.phi)
        self.assertAlmostEqual(initial_euler.theta, copy_euler.theta)
        self.assertAlmostEqual(initial_euler.psi, copy_euler.psi)

    def test_custom_deg_copy(self):
        initial_euler = i3assist.Euler(phi=56.0, theta=53.0, psi=-304.0,
                                       unit="deg")
        copy_euler = initial_euler.copy()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(initial_euler.unit, copy_euler.unit)
        self.assertAlmostEqual(initial_euler.phi, copy_euler.phi)
        self.assertAlmostEqual(initial_euler.theta, copy_euler.theta)
        self.assertAlmostEqual(initial_euler.psi, copy_euler.psi)

    def test_custom_rad_copy(self):
        initial_euler = i3assist.Euler(phi=5.27, theta=1.23, psi=-4.64,
                                         unit="rad")
        copy_euler = initial_euler.copy()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(initial_euler.unit, copy_euler.unit)
        self.assertAlmostEqual(initial_euler.phi, copy_euler.phi)
        self.assertAlmostEqual(initial_euler.theta, copy_euler.theta)
        self.assertAlmostEqual(initial_euler.psi, copy_euler.psi)

    def test_default_deg_degrees_copy(self):
        initial_euler = i3assist.Euler()
        copy_euler = initial_euler.degrees()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(initial_euler.unit, copy_euler.unit)
        self.assertAlmostEqual(initial_euler.phi, copy_euler.phi)
        self.assertAlmostEqual(initial_euler.theta, copy_euler.theta)
        self.assertAlmostEqual(initial_euler.psi, copy_euler.psi)

    def test_default_deg_degrees_inplace(self):
        initial_euler = i3assist.Euler()
        initial_unit = initial_euler.unit
        initial_phi = initial_euler.phi
        initial_theta = initial_euler.theta
        initial_psi = initial_euler.psi
        initial_euler.degrees(inplace=True)
        self.assertEqual(initial_euler.unit, initial_unit)
        self.assertAlmostEqual(initial_euler.phi, initial_phi)
        self.assertAlmostEqual(initial_euler.theta, initial_theta)
        self.assertAlmostEqual(initial_euler.psi, initial_psi)

    def test_default_rad_degrees_copy(self):
        initial_euler = i3assist.Euler(unit="rad")
        copy_euler = initial_euler.degrees()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(copy_euler.unit, "deg")
        self.assertAlmostEqual(copy_euler.phi,
                               numpy.degrees(initial_euler.phi))
        self.assertAlmostEqual(copy_euler.theta,
                               numpy.degrees(initial_euler.theta))
        self.assertAlmostEqual(copy_euler.psi,
                               numpy.degrees(initial_euler.psi))

    def test_default_rad_degrees_inplace(self):
        initial_euler = i3assist.Euler(unit="rad")
        initial_unit = initial_euler.unit
        initial_phi = initial_euler.phi
        initial_theta = initial_euler.theta
        initial_psi = initial_euler.psi
        initial_euler.degrees(inplace=True)
        self.assertEqual(initial_euler.unit, "deg")
        self.assertAlmostEqual(initial_euler.phi, numpy.degrees(initial_phi))
        self.assertAlmostEqual(initial_euler.theta,
                               numpy.degrees(initial_theta))
        self.assertAlmostEqual(initial_euler.psi, numpy.degrees(initial_psi))
    
    def test_custom_deg_degrees_copy(self):
        initial_euler = i3assist.Euler(phi=56.0, theta=53.0, psi=-304.0,
                                       unit="deg")
        copy_euler = initial_euler.degrees()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(initial_euler.unit, copy_euler.unit)
        self.assertAlmostEqual(initial_euler.phi, copy_euler.phi)
        self.assertAlmostEqual(initial_euler.theta, copy_euler.theta)
        self.assertAlmostEqual(initial_euler.psi, copy_euler.psi)

    def test_custom_deg_degrees_inplace(self):
        initial_euler = i3assist.Euler(phi=56.0, theta=53.0, psi=-304.0,
                                       unit="deg")
        initial_unit = initial_euler.unit
        initial_phi = initial_euler.phi
        initial_theta = initial_euler.theta
        initial_psi = initial_euler.psi
        initial_euler.degrees(inplace=True)
        self.assertEqual(initial_euler.unit, initial_unit)
        self.assertAlmostEqual(initial_euler.phi, initial_phi)
        self.assertAlmostEqual(initial_euler.theta, initial_theta)
        self.assertAlmostEqual(initial_euler.psi, initial_psi)

    def test_custom_rad_degrees_copy(self):
        initial_euler = i3assist.Euler(phi=5.27, theta=1.23, psi=-4.64,
                                       unit="rad")
        copy_euler = initial_euler.degrees()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(copy_euler.unit, "deg")
        self.assertAlmostEqual(copy_euler.phi,
                               numpy.degrees(initial_euler.phi))
        self.assertAlmostEqual(copy_euler.theta,
                               numpy.degrees(initial_euler.theta))
        self.assertAlmostEqual(copy_euler.psi,
                               numpy.degrees(initial_euler.psi))

    def test_custom_rad_degrees_inplace(self):
        initial_euler = i3assist.Euler(phi=5.27, theta=1.23, psi=-4.64,
                                       unit="rad")
        initial_unit = initial_euler.unit
        initial_phi = initial_euler.phi
        initial_theta = initial_euler.theta
        initial_psi = initial_euler.psi
        initial_euler.degrees(inplace=True)
        self.assertEqual(initial_euler.unit, "deg")
        self.assertAlmostEqual(initial_euler.phi, numpy.degrees(initial_phi))
        self.assertAlmostEqual(initial_euler.theta,
                               numpy.degrees(initial_theta))
        self.assertAlmostEqual(initial_euler.psi, numpy.degrees(initial_psi))

    def test_default_deg_radians_copy(self):
        initial_euler = i3assist.Euler()
        copy_euler = initial_euler.radians()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(copy_euler.unit, "rad")
        self.assertAlmostEqual(copy_euler.phi,
                               numpy.radians(initial_euler.phi))
        self.assertAlmostEqual(copy_euler.theta,
                               numpy.radians(initial_euler.theta))
        self.assertAlmostEqual(copy_euler.psi,
                               numpy.radians(initial_euler.psi))

    def test_default_deg_radians_inplace(self):
        initial_euler = i3assist.Euler()
        initial_unit = initial_euler.unit
        initial_phi = initial_euler.phi
        initial_theta = initial_euler.theta
        initial_psi = initial_euler.psi
        initial_euler.radians(inplace=True)
        self.assertEqual(initial_euler.unit, "rad")
        self.assertAlmostEqual(initial_euler.phi, numpy.radians(initial_phi))
        self.assertAlmostEqual(initial_euler.theta,
                               numpy.radians(initial_theta))
        self.assertAlmostEqual(initial_euler.psi, numpy.radians(initial_psi))

    def test_default_rad_radians_copy(self):
        initial_euler = i3assist.Euler(unit="rad")
        copy_euler = initial_euler.radians()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(initial_euler.unit, copy_euler.unit)
        self.assertAlmostEqual(initial_euler.phi, copy_euler.phi)
        self.assertAlmostEqual(initial_euler.theta, copy_euler.theta)
        self.assertAlmostEqual(initial_euler.psi, copy_euler.psi)

    def test_default_rad_radians_inplace(self):
        initial_euler = i3assist.Euler(unit="rad")
        initial_unit = initial_euler.unit
        initial_phi = initial_euler.phi
        initial_theta = initial_euler.theta
        initial_psi = initial_euler.psi
        initial_euler.radians(inplace=True)
        self.assertEqual(initial_euler.unit, initial_unit)
        self.assertAlmostEqual(initial_euler.phi, initial_phi)
        self.assertAlmostEqual(initial_euler.theta, initial_theta)
        self.assertAlmostEqual(initial_euler.psi, initial_psi)
    
    def test_custom_deg_radians_copy(self):
        initial_euler = i3assist.Euler(phi=56.0, theta=53.0, psi=-304.0,
                                       unit="deg")
        copy_euler = initial_euler.radians()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(copy_euler.unit, "rad")
        self.assertAlmostEqual(copy_euler.phi,
                               numpy.radians(initial_euler.phi))
        self.assertAlmostEqual(copy_euler.theta,
                               numpy.radians(initial_euler.theta))
        self.assertAlmostEqual(copy_euler.psi,
                               numpy.radians(initial_euler.psi))

    def test_custom_deg_radians_inplace(self):
        initial_euler = i3assist.Euler(phi=56.0, theta=53.0, psi=-304.0,
                                       unit="deg")
        initial_unit = initial_euler.unit
        initial_phi = initial_euler.phi
        initial_theta = initial_euler.theta
        initial_psi = initial_euler.psi
        initial_euler.radians(inplace=True)
        self.assertEqual(initial_euler.unit, "rad")
        self.assertAlmostEqual(initial_euler.phi, numpy.radians(initial_phi))
        self.assertAlmostEqual(initial_euler.theta,
                               numpy.radians(initial_theta))
        self.assertAlmostEqual(initial_euler.psi, numpy.radians(initial_psi))

    def test_custom_rad_radians_copy(self):
        initial_euler = i3assist.Euler(phi=5.27, theta=1.23, psi=-4.64,
                                       unit="rad")
        copy_euler = initial_euler.radians()
        self.assertIsNot(initial_euler, copy_euler)
        self.assertEqual(initial_euler.unit, copy_euler.unit)
        self.assertAlmostEqual(initial_euler.phi, copy_euler.phi)
        self.assertAlmostEqual(initial_euler.theta, copy_euler.theta)
        self.assertAlmostEqual(initial_euler.psi, copy_euler.psi)

    def test_custom_rad_radians_inplace(self):
        initial_euler = i3assist.Euler(phi=5.27, theta=1.23, psi=-4.64,
                                       unit="rad")
        initial_unit = initial_euler.unit
        initial_phi = initial_euler.phi
        initial_theta = initial_euler.theta
        initial_psi = initial_euler.psi
        initial_euler.radians(inplace=True)
        self.assertEqual(initial_euler.unit, initial_unit)
        self.assertAlmostEqual(initial_euler.phi, initial_phi)
        self.assertAlmostEqual(initial_euler.theta, initial_theta)
        self.assertAlmostEqual(initial_euler.psi, initial_psi)

    def test_default_deg_to_string(self):
        initial_euler = i3assist.Euler()
        self.assertEqual(str(initial_euler),
                         "{: f} {: f} {: f}".format(0., 0., 0.))

    def test_default_rad_to_string(self):
        initial_euler = i3assist.Euler(unit="rad")
        self.assertEqual(str(initial_euler),
                         "{: f} {: f} {: f}".format(0., 0., 0.))

    def test_custom_deg_to_string(self):
        initial_euler = i3assist.Euler(phi=56.0, theta=53.0, psi=-304.0,
                                       unit="deg")
        self.assertEqual(str(initial_euler),
                         "{: f} {: f} {: f}".format(56., 53., -304.))

    def test_custom_rad_to_string(self):
        initial_euler = i3assist.Euler(phi=5.27, theta=1.23, psi=-4.64,
                                       unit="rad")
        self.assertEqual(str(initial_euler),
                         "{: f} {: f} {: f}".format(5.27, 1.23, -4.64))

    def test_default_deg_pos_string(self):
        initial_euler = i3assist.Euler()
        self.assertEqual(initial_euler.pos_string(),
                         "{: 10.4f}{: 10.4f}{: 10.4f}".format(0., 0., 0.))

    def test_default_rad_pos_string(self):
        initial_euler = i3assist.Euler(unit="rad")
        self.assertEqual(initial_euler.pos_string(),
                         "{: 10.4f}{: 10.4f}{: 10.4f}".format(0., 0., 0.))

    def test_custom_deg_pos_string(self):
        initial_euler = i3assist.Euler(phi=56.0, theta=53.0, psi=-304.0,
                                       unit="deg")
        self.assertEqual(initial_euler.pos_string(),
                         "{: 10.4f}{: 10.4f}{: 10.4f}".format(56., 53., -304.))

    def test_custom_rad_pos_string(self):
        initial_euler = i3assist.Euler(phi=5.27, theta=1.23, psi=-4.64,
                                       unit="rad")
        self.assertEqual(initial_euler.pos_string(),
                         "{: 10.4f}{: 10.4f}{: 10.4f}".format(5.27, 1.23,
                                                              -4.64))

if __name__ == "__main__":
    unittest.main()
