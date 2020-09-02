import Simulation
import numpy
import configparser
import pytest

# Tests concerning coefficient_hi_generator()

# RANDOM COEFFICIENTS


def test_shape_random_coefficient_hi_generator():
    """ This test ensures that the random hi coefficients are generated with
        a proper shape.

        Asserts:
            That the returning value of the function
            random_coefficient_hi_generator() is a numpy array of shape
            (number_of_particles,).
    """
    assert numpy.shape(
        Simulation.random_coefficient_hi_generator(
            3)) == (
        3,)


def test_real_random_coefficient_hi_generator():
    """ This test ensures that the random hi coefficients are real numbers.

        Asserts:
            That the returning value of the function
            random_coefficient_hi_generator() is an array of real numbers.
    """
    i = 0
    while i < 3:
        assert isinstance(
            Simulation.random_coefficient_hi_generator(
                3)[i], float)
        i += 1

# NON RANDOM COEFFICIENTS


def test_shape_coefficient_hi_reader():
    """ This test ensures that the hi coefficients read from the file
        Simulation_parameters.ini are generated with a proper shape.

        Asserts:
            That the returning value of the function coefficient_hi_reader() is a
            numpy array of shape (number_of_particles,).
    """
    parameters = configparser.ConfigParser()
    parameters.read("Simulation_parameters.ini")
    number_of_particles = int(
        parameters.get(
            "Parameters",
            "Number of Particles"))
    hi = parameters.get("Coefficients", "hi_coefficients")
    assert numpy.shape(
        Simulation.coefficient_hi_reader(
            number_of_particles, hi)) == (
        number_of_particles,)


def test_real_coefficient_hi_reader():
    """ This test ensures that the hi coefficients read from the file
        Simulation_parameters.ini are generated as real numbers.

        Asserts:
            That the returning value of the function coefficient_hi_reader() is an
            array of real numbers.
    """
    parameters = configparser.ConfigParser()
    parameters.read("Simulation_parameters.ini")
    number_of_particles = int(
        parameters.get(
            "Parameters",
            "Number of Particles"))
    hi = parameters.get("Coefficients", "hi_coefficients")
    i = 0
    while i < number_of_particles:
        assert isinstance(
            Simulation.coefficient_hi_reader(
                number_of_particles, hi)[i], float)
        i += 1
# -----------------------------------------------------------------------------

# Tests concerning coefficient_Jij_generator()

# RANDOM COEFFICIENTS


def test_shape_random_coefficient_Jij_generator():
    """ This test ensures that the random Jij coefficients are generated with
        a proper shape.

        Asserts:
            That the returning value of the function
            random_coefficient_Jij_generator() is a numpy array of shape
            (number_of_particles,number_of_particles).
    """
    assert numpy.shape(
        Simulation.random_coefficient_Jij_generator(
            3)) == (
        3, 3)


def test_matrix_random_diagonal_coefficient_Jij_generator():
    """ This test ensures that the Jii coefficients are 0.

        Asserts:
            That the returning value of the function
            random_coefficient_Jij_generator() is a matrix with 0s in the diagonal.
    """
    i = 0
    generated_matrix = Simulation.random_coefficient_Jij_generator(3)
    while i < 3:
        j = 0
        while j < 3:
            if i == j:
                assert generated_matrix[i, j] == 0.
            j += 1
        i += 1


def test_matrix_real_random_coefficient_Jij_generator():
    """ This test ensures that the random Jij coefficients are real numbers.

        Asserts:
            That the returning value of the function
            random_coefficient_Jij_generator() is a matrix of real numbers.
    """
    i = 0
    generated_matrix = Simulation.random_coefficient_Jij_generator(3)
    while i < 3:
        j = 0
        while j < 3:
            if i != j:
                assert isinstance(generated_matrix[i, j], float)
            j += 1
        i += 1


def test_symmetry_random_coefficient_Jij_generator():
    """ This test ensures that the coefficients Jji=Jij.

        Asserts:
            That the returning value of the function
            random_coefficient_Jij_generator() is a symmetric matrix.
    """
    i = 0
    generated_matrix = Simulation.random_coefficient_Jij_generator(3)
    while i < 3:
        j = 0
        while j < 3:
            if i != j:
                assert generated_matrix[i, j] == generated_matrix[j, i]
            j += 1
        i += 1

# NON RANDOM COEFFICIENTS


def test_shape_coefficient_Jij_reader():
    """ This test ensures that the Jij coefficients read from the file
        Simulation_parameters.ini are generated with a proper shape.

        Asserts:
            That the returning value of the function coefficient_Jij_reader()
            is a matrix of shape (number_of_particles,number_of_particles).
    """
    parameters = configparser.ConfigParser()
    parameters.read("Simulation_parameters.ini")
    number_of_particles = int(
        parameters.get(
            "Parameters",
            "Number of Particles"))
    Jij = parameters.get("Coefficients", "Jij_coefficients")
    assert numpy.shape(
        Simulation.coefficient_Jij_reader(
            number_of_particles,
            Jij)) == (
        number_of_particles,
        number_of_particles)


def test_matrix_diagonal_coefficient_Jij_reader():
    """ This test ensures that the Jii coefficients are 0.

        Asserts:
            That the returning value of the function coefficient_Jij_reader()
            is a matrix with 0s in the diagonal.
    """
    parameters = configparser.ConfigParser()
    parameters.read("Simulation_parameters.ini")
    number_of_particles = int(
        parameters.get(
            "Parameters",
            "Number of Particles"))
    Jij = parameters.get("Coefficients", "Jij_coefficients")
    i = 0
    generated_matrix = Simulation.coefficient_Jij_reader(
        number_of_particles, Jij)
    while i < number_of_particles:
        j = 0
        while j < number_of_particles:
            if i == j:
                assert generated_matrix[i, j] == 0.
            j += 1
        i += 1


def test_matrix_real_coefficient_Jij_reader():
    """ This test ensures that the Jij coefficients are real numbers.

        Asserts:
            That the returning value of the function coefficient_Jij_reader()
            is a matrix of real numbers.
    """
    parameters = configparser.ConfigParser()
    parameters.read("Simulation_parameters.ini")
    number_of_particles = int(
        parameters.get(
            "Parameters",
            "Number of Particles"))
    Jij = parameters.get("Coefficients", "Jij_coefficients")
    i = 0
    generated_matrix = Simulation.coefficient_Jij_reader(
        number_of_particles, Jij)
    while i < number_of_particles:
        j = 0
        while j < number_of_particles:
            if i != j:
                assert isinstance(generated_matrix[i, j], float)
            j += 1
        i += 1


def test_symmetry_coefficient_Jij_reader():
    """ This test ensures that the coefficients Jji=Jij.

        Asserts:
            That the returning value of the function coefficient_Jij_reader()
            is a symmetric matrix.
    """
    parameters = configparser.ConfigParser()
    parameters.read("Simulation_parameters.ini")
    number_of_particles = int(
        parameters.get(
            "Parameters",
            "Number of Particles"))
    Jij = parameters.get("Coefficients", "Jij_coefficients")
    i = 0
    generated_matrix = Simulation.coefficient_Jij_reader(
        number_of_particles, Jij)
    while i < number_of_particles:
        j = 0
        while j < number_of_particles:
            if i != j:
                assert generated_matrix[i, j] == generated_matrix[j, i]
            j += 1
        i += 1

# -----------------------------------------------------------------------------


def test_btest_1():
    """This test ensures that the function btest(i,n) works properly.

       Asserts:
            That the returning value of btest(100,2) is 1.0, which means that
            the third bit of the number 100 is 1.
    """
    assert Simulation.btest(100, 2) == 1.0


def test_btest_2():
    """This test ensures that the function btest(i,n) works properly.

       Asserts:
            That the returning value of btest(99,2) is -1.0, which means that
            the third bit of the number 99 is 0.
    """
    assert Simulation.btest(99, 2) == -1, 0

# -----------------------------------------------------------------------------


def test_Hamiltonian_1():
    """This test ensures that the target Hamiltonian is generated as it should.

       Asserts:
           That the returning value of Hamiltonian_1(2,hi,Jij), given certain
           hi and Jij coefficient values, is the one that is expected.
    """
    hi = numpy.array((1, 1))
    Jij = numpy.array(((0, 1), (1, 0)))
    expected_result = numpy.array(
        ((0., 0., 0., 0.), (0., 2., 0., 0.), (0., 0., 2., 0.), (0., 0., 0., 4.)))
    actual_result = Simulation.Hamiltonian_1(2, hi, Jij)
    assert actual_result.all() == expected_result.all()

# -----------------------------------------------------------------------------


def test_Hamiltonian_0():
    """This test ensures that the initial Hamiltonian is generated as it should.

       Asserts:
           That the returning value of Hamiltonian_0(2), is the one that is
           expected.
    """
    expected_result = numpy.array(
        ((0., 1., 1., 0.), (1., 0., 0., 1.), (1., 0., 0., 1.), (0., 1., 1., 0.)))
    actual_result = Simulation.Hamiltonian_0(2)
    assert actual_result.all() == expected_result.all()

# -----------------------------------------------------------------------------


def test_gap_two_smallest_in_array():
    """This test ensures that the function gap_two_smallest_in_array()
       calculates the difference between the two smallest elements in an array
       properly.
    """
    array = [0, 2, 5, 9, 14]
    expected_result = 2
    actual_result = Simulation.gap_two_smallest_in_array(array)
    assert actual_result == expected_result

# -----------------------------------------------------------------------------


def test_position_smallest_in_array():
    """This test ensures that the function position_smallest_in_array()
       locates the smallest element in an array properly.
    """
    array = [2, 5, 9, 0, 14]
    expected_result = 3
    actual_result = Simulation.position_smallest_in_array(array)
    assert actual_result == expected_result

# -----------------------------------------------------------------------------


def test_gap_quantum_simulation():
    """This test ensures that the evolution of the gap is properly calculated
       during the simulation.

       Asserts:
           That the function quantum_simulation(2,2,H0,H1) returns the expected
           value of the evolution of the gap, in a simulation with 3 steps
           (initial, middle and final), 2 particles and certain hi and Jij
           coefficients.
    """
    hi = numpy.array((1, 1))
    Jij = numpy.array(((0, 1), (1, 0)))
    H0 = Simulation.Hamiltonian_0(2)
    H1 = Simulation.Hamiltonian_1(2, hi, Jij)
    expected_result = numpy.array((2., 0.4811943, 2.))
    actual_result, _ = Simulation.quantum_simulation(2, 2, H0, H1)
    assert actual_result.all() == expected_result.all()

# -----------------------------------------------------------------------------


def test_eigenvectors_quantum_simulation():
    """This test ensures that the evolution of the probability of the states
       is properly calculated during the simulation.

       Asserts:
           That the function quantum_simulation(2,2,H0,H1) returns the expected
           value of the evolution of the probability of the states, in a
           simulation with 3 steps (initial, middle and final), 2 particles and
           certain hi and Jij coefficients.
    """
    hi = numpy.array((1, 1))
    Jij = numpy.array(((0, 1), (1, 0)))
    H0 = Simulation.Hamiltonian_0(2)
    H1 = Simulation.Hamiltonian_1(2, hi, Jij)
    expected_result = numpy.array(
        ((0.25,
          0.25,
          0.25,
          0.25),
         (0.17956835,
          0.39396156,
          0.39396156,
          0.03250853),
            (0.,
             1.,
             0.,
             0.)))
    _, actual_result = Simulation.quantum_simulation(2, 2, H0, H1)
    assert actual_result.all() == expected_result.all()

# -----------------------------------------------------------------------------


def test_computational_time():
    """This test ensures that the time required for the Adiabatic Quantum
       Computation is properly calculated.

       Asserts:
           That the function results_of_simulation(100,2,H0,H1) returns the
           expected value of the time, in a simulation with 100 steps, 2
           particles and certain hi and Jij coefficients.
    """
    hi = numpy.array((1, 2))
    Jij = numpy.array(((0, 1), (1, 0)))
    H0 = Simulation.Hamiltonian_0(2)
    H1 = Simulation.Hamiltonian_1(2, hi, Jij)
    expected_result = 0.6571
    actual_result = Simulation.results_of_simulation(100, 2, H0, H1)
    assert actual_result == expected_result
