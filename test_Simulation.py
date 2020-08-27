import Simulation
import numpy
import pytest

# Tests concerning coefficient_hi_generator()

# RANDOM COEFFICIENTS 
# Test the hi coefficient generator function to see if the return value is
# a properly shaped numpy array.


def test_shape_random_coefficient_hi_generator():
    assert numpy.shape(Simulation.coefficient_hi_generator(3, 1)) == (3,)

# Test the hi coefficient generator function to see if the return value is
# an array of real numbers.


def test_real_random_coefficient_hi_generator():
    i = 0
    while i < 3:
        assert isinstance(Simulation.coefficient_hi_generator(3, 1)[i], float)
        i += 1

# NON RANDOM COEFFICIENTS
# Test the hi coefficient generator function to see if the return value is
# a properly shaped numpy array.


def test_shape_coefficient_hi_generator():
    parameters = open("Simulation_parameters.txt", "r")
    lines = parameters.readlines()
    number_of_particles = int(lines[0])
    parameters.close()
    assert numpy.shape(
        Simulation.coefficient_hi_generator(
            number_of_particles, 0)) == (
        number_of_particles,)

# Test the hi coefficient generator function to see if the return value is
# an array of real numbers.


def test_real_coefficient_hi_generator():
    parameters = open("Simulation_parameters.txt", "r")
    lines = parameters.readlines()
    number_of_particles = int(lines[0])
    parameters.close()
    i = 0
    while i < number_of_particles:
        assert isinstance(
            Simulation.coefficient_hi_generator(
                number_of_particles, 0)[i], float)
        i += 1
# -----------------------------------------------------------------------------

# Tests concerning coefficient_Jij_generator()

# RANDOM COEFFICIENTS
# Test the Jij generator function to see if the matrix has the apropriate
# dimension.


def test_shape_random_coefficient_Jij_generator():
    assert numpy.shape(Simulation.coefficient_Jij_generator(3, 1)) == (3, 3)

# Test if the Jij coefficient matrix has 0s in the diagonal.


def test_matrix_random_diagonal_coefficient_Jij_generator():
    i = 0
    generated_matrix = Simulation.coefficient_Jij_generator(3, 1)
    while i < 3:
        j = 0
        while j < 3:
            if i == j:
                assert generated_matrix[i, j] == 0.
            j += 1
        i += 1

# Test if the Jij coefficient matrix is composed of real numbers elsewhere.


def test_matrix_real_random_coefficient_Jij_generator():
    i = 0
    generated_matrix = Simulation.coefficient_Jij_generator(3, 1)
    while i < 3:
        j = 0
        while j < 3:
            if i != j:
                assert isinstance(generated_matrix[i, j], float)
            j += 1
        i += 1

# Test if the Jij coefficient matrix is symmetric.


def test_symmetry_random_coefficient_Jij_generator():
    i = 0
    generated_matrix = Simulation.coefficient_Jij_generator(3, 1)
    while i < 3:
        j = 0
        while j < 3:
            if i != j:
                assert generated_matrix[i, j] == generated_matrix[j, i]
            j += 1
        i += 1

# NON RANDOM COEFFICIENTS
# Test the Jij generator function to see if the matrix has the apropriate
# dimension.


def test_shape_coefficient_Jij_generator():
    parameters = open("Simulation_parameters.txt", "r")
    lines = parameters.readlines()
    number_of_particles = int(lines[0])
    parameters.close()
    assert numpy.shape(
        Simulation.coefficient_Jij_generator(
            number_of_particles,
            0)) == (
        number_of_particles,
        number_of_particles)

# Test if the Jij coefficient matrix has 0s in the diagonal.


def test_matrix_diagonal_coefficient_Jij_generator():
    parameters = open("Simulation_parameters.txt", "r")
    lines = parameters.readlines()
    number_of_particles = int(lines[0])
    parameters.close()
    i = 0
    generated_matrix = Simulation.coefficient_Jij_generator(
        number_of_particles, 0)
    while i < number_of_particles:
        j = 0
        while j < number_of_particles:
            if i == j:
                assert generated_matrix[i, j] == 0.
            j += 1
        i += 1

# Test if the Jij coefficient matrix is composed of real numbers elsewhere.


def test_matrix_real_coefficient_Jij_generator():
    parameters = open("Simulation_parameters.txt", "r")
    lines = parameters.readlines()
    number_of_particles = int(lines[0])
    parameters.close()
    i = 0
    generated_matrix = Simulation.coefficient_Jij_generator(
        number_of_particles, 0)
    while i < number_of_particles:
        j = 0
        while j < number_of_particles:
            if i != j:
                assert isinstance(generated_matrix[i, j], float)
            j += 1
        i += 1

# Test if the Jij coefficient matrix is symmetric.


def test_symmetry_coefficient_Jij_generator():
    parameters = open("Simulation_parameters.txt", "r")
    lines = parameters.readlines()
    number_of_particles = int(lines[0])
    parameters.close()
    i = 0
    generated_matrix = Simulation.coefficient_Jij_generator(
        number_of_particles, 0)
    while i < number_of_particles:
        j = 0
        while j < number_of_particles:
            if i != j:
                assert generated_matrix[i, j] == generated_matrix[j, i]
            j += 1
        i += 1

# -----------------------------------------------------------------------------

# Test btest function, to see if it returns 1.0 or -1.0 when it should.


def test_btest_1():
    assert Simulation.btest(100, 2) == 1.0


def test_btest_2():
    assert Simulation.btest(99, 2) == -1, 0

# -----------------------------------------------------------------------------


def test_Hamiltonian_1():
    hi = numpy.array((1, 1))
    Jij = numpy.array(((0, 1), (1, 0)))
    expected_result = numpy.array(
        ((0., 0., 0., 0.), (0., 2., 0., 0.), (0., 0., 2., 0.), (0., 0., 0., 4.)))
    actual_result = Simulation.Hamiltonian_1(2, hi, Jij)
    assert actual_result.all() == expected_result.all()

# -----------------------------------------------------------------------------


def test_Hamiltonian_0():
    expected_result = numpy.array(
        ((0., 1., 1., 0.), (1., 0., 0., 1.), (1., 0., 0., 1.), (0., 1., 1., 0.)))
    actual_result = Simulation.Hamiltonian_0(2)
    assert actual_result.all() == expected_result.all()

# -----------------------------------------------------------------------------


def test_gap_two_smallest_in_array():
    array = [0, 2, 5, 9, 14]
    expected_result = 2
    actual_result = Simulation.gap_two_smallest_in_array(array)
    assert actual_result == expected_result

# -----------------------------------------------------------------------------


def test_position_smallest_in_array():
    array = [2, 5, 9, 0, 14]
    expected_result = 3
    actual_result = Simulation.position_smallest_in_array(array)
    assert actual_result == expected_result

# -----------------------------------------------------------------------------


def test_gap_quantum_simulation():
    hi = numpy.array((1, 1))
    Jij = numpy.array(((0, 1), (1, 0)))
    H0 = Simulation.Hamiltonian_0(2)
    H1 = Simulation.Hamiltonian_1(2, hi, Jij)
    expected_result = numpy.array((2., 0.4811943, 2.))
    actual_result, _ = Simulation.quantum_simulation(2, 2, H0, H1)
    assert actual_result.all() == expected_result.all()

# -----------------------------------------------------------------------------


def test_eigenvectors_quantum_simulation():
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
    hi = numpy.array((1, 2))
    Jij = numpy.array(((0, 1), (1, 0)))
    H0 = Simulation.Hamiltonian_0(2)
    H1 = Simulation.Hamiltonian_1(2, hi, Jij)
    expected_result = 0.6571
    actual_result = Simulation.results_of_simulation(100, 2, H0, H1)
    assert actual_result == expected_result
