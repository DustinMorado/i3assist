import unittest
import i3assist
import numpy

class TestRotationMatrixMethods(unittest.TestCase):

    def test_default_deg_initialization(self):
        initial_matrix = i3assist.RotationMatrix()
        self.assertTrue(numpy.allclose(initial_matrix.matrix,
                                       numpy.identity(3)))

    def test_default_rad_initialization(self):
        initial_matrix = i3assist.RotationMatrix()
        self.assertTrue(numpy.allclose(initial_matrix.matrix,
                                       numpy.identity(3)))

    def test_custom_deg_vector_initialization(self):
        initial_matrix = i3assist.RotationMatrix([56.0, 53.0, -304.0])
        i3euler_output = numpy.array(
            [-0.1009327461287336, 0.7425885137345521, 0.6620988446058652,
             -0.742588513734552, -0.4991177229766853, 0.4465913096781866,
             0.6620988446058653, -0.4465913096781864, 0.6018150231520484],
            numpy.float_)
        self.assertTrue(numpy.allclose(initial_matrix.matrix.reshape(-1),
                                       i3euler_output))

    def test_custom_rad_vector_initialization(self):
        initial_matrix = i3assist.RotationMatrix([5.27, 1.23, -4.64],
                                                 unit="rad")
        i3euler_output = numpy.array(
            [0.2445932562414799, 0.2377722286041004, 0.9400204818544271,
             -0.5482869348823414, 0.8335075294132077, -0.0681662339345359,
             -0.7997221867864596, -0.4987279475988188, 0.3342377271424531],
            numpy.float_)
        self.assertTrue(numpy.allclose(initial_matrix.matrix.reshape(-1),
                                       i3euler_output))

    def test_custom_deg_matrix_initialization(self):
        initial_array = [0.2445932562414799, 0.2377722286041004,
                         0.9400204818544271, -0.5482869348823414,
                         0.8335075294132077, -0.0681662339345359,
                         -0.7997221867864596, -0.4987279475988188,
                         0.3342377271424531]
        initial_matrix = i3assist.RotationMatrix(initial_array)
        self.assertTrue(numpy.allclose(initial_matrix.matrix.reshape(-1),
                                       initial_array))

    def test_custom_rad_matrix_initialization(self):
        initial_array = [0.2445932562414799, 0.2377722286041004,
                         0.9400204818544271, -0.5482869348823414,
                         0.8335075294132077, -0.0681662339345359,
                         -0.7997221867864596, -0.4987279475988188,
                         0.3342377271424531]
        initial_matrix = i3assist.RotationMatrix(initial_array, unit="rad")
        self.assertTrue(numpy.allclose(initial_matrix.matrix.reshape(-1),
                                       initial_array))

    def test_default_to_string(self):
        pass

    def test_custom_vector_to_string(self):
        pass

    def test_custom_matrix_to_string(self):
        pass

if __name__ == "__main__":
    unittest.main()
