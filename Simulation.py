import numpy as np
import scipy.linalg as la
import scipy.integrate as integrate


def coefficient_hi_reader(number_of_particles, hi_coefficients):
    """This function reads a set of coefficients (hi) used to generate
       Hamiltonian H1. They represent the local magnetic fields.

    Parameters:
        number_of_particles : number of particles present in the system of
                              particles.
        hi_coefficients : string that contains the values of the coefficients,
                          separated by commas.

    Returns:
        The hi coefficients in form of a vector of dimension number_of_particles.

    Raise:
        ValueError if :
            number_of_particles is not between 2 and 8.
            The dimension of the vector is greater than number_of_particles.
    """
    if abs(number_of_particles) > 8 or abs(number_of_particles) < 2:
        raise ValueError(
            "You must insert as number of particles a positive integer lower or equal to 8 and higher than 1.")
    else:
        hi = np.zeros(number_of_particles)
        numbers = hi_coefficients.split(",")
        if len(numbers) != number_of_particles:
            raise ValueError(
                "The dimension of the vector must be equal to the number of particles.")
        i = 0
        for number in numbers:
            num = float(number)
            hi[i] = num
            i += 1
    return(hi)


def coefficient_Jij_reader(number_of_particles, Jij_coefficients):
    """This function reads a set of coefficients (Jij) used to generate
       Hamiltonian H1. They represent the matching force between spins.

    Parameters:
        number_of_particles : number of particles present in the system of
                              particles.
        Jij_coefficients : string that contains the values of the coefficients,
                           separated by commas and slashes.

    Returns:
        The Jij coefficients in form of a matrix of dimension
        number_of_particles x number_of_particles.

    Raise:
        ValueError if:
            The dimension of the matrix is not consistent with number_of_particles.
            The diagonal of the matrix is not comprised of 0s.
            The matrix is not symmetric.
    """
    Jij = np.zeros((number_of_particles, number_of_particles))
    rows = Jij_coefficients.split("/")
    if len(rows) != number_of_particles:
        raise ValueError(
            "The number of rows must be consistent with the number of particles.")
    i = 0
    while i < number_of_particles:
        numbers = rows[i].split(",")
        if len(numbers) != number_of_particles:
            raise ValueError(
                "The number of columns must be consistent with the number of particles.")
        j = 0
        for number in numbers:
            num = float(number)
            Jij[i, j] = num
            if i == j:
                if Jij[i, i] != 0.:
                    raise ValueError(
                        "The diagonal of the J_ij matrix must contain 0s. Check the parameters file.")
            elif i > j:
                if Jij[i, j] != Jij[j, i]:
                    raise ValueError(
                        "The matrix of J_ij coefficients should be symmetric. Check the parameters file.")
            j += 1
        i += 1
    return(Jij)


def random_coefficient_hi_generator(number_of_particles):
    """This function randomly generates a set of coefficients (hi) used to
       generate Hamiltonian H1. They represent the local magnetic fields.

    Parameters:
        number_of_particles : number of particles present in the system of
                              particles.

    Returns:
        The hi coefficients in form of a vector of dimension number_of_particles.

    Raise:
        ValueError if :
            number_of_particles is not between 2 and 8.
    """
    if abs(number_of_particles) > 8 or abs(number_of_particles) < 2:
        raise ValueError(
            "You must insert as number of particles a positive integer lower or equal to 8 and higher than 1.")
    else:
        # Random values of the hi coefficients.
        hi = np.random.rand(number_of_particles)
        return(hi)


def random_coefficient_Jij_generator(number_of_particles):
    """This function randomly generates a set of coefficients (Jij) used to
       generate Hamiltonian H1. They represent the matching force between spins.

    Parameters:
        number_of_particles : number of particles present in the system of
                              particles.

    Returns:
        The Jij coefficients in form of a matrix of dimension
        number_of_particles x number_of_particles.
    """
    # Random values of the Jij coefficients.
    Jij = np.random.rand(number_of_particles, number_of_particles)
    i = 0
    while i < number_of_particles:
        j = 0
        while j < number_of_particles:
            if i == j:
                # The coupling coefficient is 0 if the index is equal.
                Jij[i, j] = 0.
            if i < j:
                # The coupling coefficient Jji is equal to Jij.
                Jij[j, i] = Jij[i, j]
            j += 1
        i += 1
    return(Jij)


def btest(i, n):
    """This function determines if the n-th bit of a number is 0 or 1. It is
       used for the generation of both Hamiltonians.

    Parameters:
        i : number.
        n : n-th bit of the number i.

    Returns:
        +1 if the n-th bit of i is 1.
        -1 if the n-th bit of i is 0.
    """
    if (i & (1 << n)):
        return(1.)
        # n-th bit is set (1).

    else:
        return(-1.)
        # n-th bit is not set (0).


def Hamiltonian_1(number_of_particles, hi, Jij):
    """This function creates the target Hamiltonian towards which the system
       will evolve.

    Parameters:
        number_of_particles : number of particles present in the system of
                              particles.
        hi : the hi coefficients representing the local magnetic fields.
        Jij : the Jij coefficients representing the matching force between spins.

    Returns:
        A diagonal matrix of dimension 2^number_of_particles x 2^number_of_particles
        representing the target Hamiltonian H1.
    """
    # Initialise H1 with 0s to give value later.
    H1 = np.zeros((2**number_of_particles, 2**number_of_particles))
    i = 0
    while i < 2**(number_of_particles):
        j = 0
        summation = 0.
        while j < number_of_particles:
            k = 0
            while k < number_of_particles:
                # The Sz operator acts multipliying by 1 if the qubit is 1 and
                # by -1 if it is 0.
                summation = summation + btest(i, j) * btest(i, k) * Jij[j, k]
                k = k + 1
            summation = summation + btest(i, j) * hi[j]
            j = j + 1
        H1[i, i] = summation
        i = i + 1
    return(H1)


def Hamiltonian_0(number_of_particles):
    """This function creates the initial Hamiltonian from which the system will
       evolve.

    Parameters:
        number_of_particles : number of particles present in the system of
                              particles.

    Returns:
        A symmetric matrix of dimension 2^number_of_particles x 2^number_of_particles
        composed of 0s and -1s, representing the initial Hamiltonian.

    """
    # Initialise H0 with 0s to give value later.
    H0 = np.zeros((2**number_of_particles, 2**number_of_particles))
    i = 0
    while i < 2**number_of_particles:
        j = 0
        while j < number_of_particles:
            # The Sx operator acts changing the qubit value 0 <--> 1.
            if btest(i, j) == 1.:
                H0[i, i - 2**j] = 1
                H0[i - 2**j, i] = 1
            else:
                H0[i, i + 2**j] = 1
                H0[i + 2**j, i] = 1
            j = j + 1
        i = i + 1
    H0 = -H0
    return(H0)


def gap_two_smallest_in_array(array):
    """This function calculates the difference between the two smallest numbers
       in an array. It is used to obtain the Energy Gap.

    Parameters:
        array : array of numbers.

    Returns:
        The difference between the two smallest numbers in array.
    """
    array_size = len(array)
    first = second = float("inf")
    for i in range(0, array_size):
        # If current element is smaller than first then
        # update both first and second.
        if array[i] < first:
            second = first
            first = array[i]
        # If arr[i] is in between first and second then
        # update second.
        elif (array[i] < second and array[i] != first):
            second = array[i]
    gap = second - first
    return(gap)


def position_smallest_in_array(array):
    """This function calculates the position of the smallest number in an array.
       It is used to locate the ground state.

    Parameters:
        array : array of numbers.

    Returns:
        The position of the smallest number in array.
    """

    array_size = len(array)
    smallest = float("inf")
    for i in range(0, array_size):
        # If current element is smaller than smallest, update smallest.
        if array[i] < smallest:
            smallest = array[i]
            position = i
    return(position)


def quantum_simulation(num_steps, number_of_particles, H0, H1):
    """This function performs the simulation, evolving from the initial to the
       target Hamiltonian, calculating in each step the Energy Gap and the
       probability of the states.

    Parameters:
        num_steps : number of steps when evolving from H0 to H1.
        number_of_particles : number of particles present in the system of
                              particles.
        H0 : Initial Hamiltonian.
        H1 : Target Hamiltonian, towards which the system evolves.

    Returns:
        A vector containing the value of the Gap in each step.
        A matrix containing the probability of each state for every step.
    """
    step = 1.0 / (float(num_steps))
    s = 0.0  # Adiabatic parameter that evolves from 0 to 1.
    i = 0
    # Initialise gap and eigenvectors with 0s to give values later.
    gap = np.zeros(num_steps + 1)
    eigenvectors = np.zeros((num_steps + 1, 2**number_of_particles))
    while i <= num_steps:
        # H_s evolves from H0 to H1.
        H_s = s * H1 + (1.0 - s) * H0
        # Obtaining the eigenvalues and eigenvectors the ground state is known.
        eigenvals, eigenvecs = la.eig(H_s)
        eigenvals = eigenvals.real
        gap[i] = gap_two_smallest_in_array(eigenvals)
        ground_state_position = position_smallest_in_array(eigenvals)
        eigenvectors[i] = np.array(eigenvecs[:, ground_state_position])
        s = s + step
        i = i + 1
    probability = eigenvectors**2
    return(gap, probability)


def results_of_simulation(num_steps, number_of_particles, H0, H1):
    """ This function calculates the time required for the AQC and saves the data
        concerning the evolution of the Gap, the speed and the probability of
        the states.

    Parameters:
        num_steps : number of steps when evolving from H0 to H1.
        number_of_particles : number of particles present in the system of
                              particles.
        H0 : Initial Hamiltonian.
        H1 : Target Hamiltonian, towards which the system evolves.

    Returns:
        The time required for the Adiabatic Quantum Computation.
    """
    # Calculating the time, integrating the inverse of the speed along the
    # adiabatic parameter s.
    gap, probability = quantum_simulation(
        num_steps, number_of_particles, H0, H1)
    xaxis = np.linspace(0., 1., 101)
    # The speed is proportional to the square of the gap.
    speed = gap**2
    speed_inverse = 1.0 / (speed)
    time = integrate.simps(speed_inverse, xaxis)
    time = round(time, 4)

    # Saving data of the Gap.
    gap_file = open("Data/Gap.txt", "w")
    np.savetxt(gap_file, gap)
    gap_file.close()

    # Saving data of the Speed.
    speed_file = open("Data/Speed.txt", "w")
    np.savetxt(speed_file, speed)
    speed_file.close()

    # Saving data of the evolution of probability.
    states_file = open("Data/States.txt", "w")
    for row in probability:
        np.savetxt(states_file, row)
    states_file.close()

    return(time)
