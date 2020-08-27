import numpy as np
import scipy.linalg as la
import scipy.integrate as integrate
import matplotlib.pyplot as plt


def coefficient_hi_generator(number_of_particles, random_or_not):
    """This function generates a set of coefficients (hi) used later by the
       Hamiltonian. They represent the local magnetic fields.

    Parameters:
        number_of_particles : number of particles present in the system of
                              particles.
        random_or_not : integer that determines if the coefficients are 
                        generated randomly or are given by the user 
                        (1=random, 0=not random).

    Returns:
        The hi coefficients in form of a vector of dimension number_of_particles.

    Raise:
        ValueError if :
            number_of_particles is not between 2 and 8.
            random_or_not is not 0 or 1.
            The dimension of the vector is greater than number_of_particles.
    """
    if abs(number_of_particles) > 8 or abs(number_of_particles) < 2:
        raise ValueError(
            "You must insert as number of particles a positive integer lower or equal to 8 and higher than 1.")
    elif random_or_not not in [0, 1]:
        raise ValueError(
            "The second line of the file containing the parameters must be 0 or 1.")
    else:
        if random_or_not == 0:
            with open("Simulation_parameters.txt", "r") as parameters:
                hi = np.zeros(number_of_particles)
                lines = parameters.read().splitlines()
                line_marker = 1
                for line in lines:
                    if line_marker == 3:
                        numbers = line.split(",")
                        if len(numbers) != number_of_particles:
                            raise ValueError(
                                "The dimension of the vector must be equal to the number of particles.")
                        k = 0
                        for number in numbers:
                            num = float(number)
                            hi[k] = num
                            k += 1
                    line_marker += 1
            return(hi)
        elif random_or_not == 1:
            hi = np.random.rand(number_of_particles)
            return(hi)


def coefficient_Jij_generator(number_of_particles, random_or_not):
    """This function generates a set of coefficients (Jij) used later by the
       Hamiltonian. They represent the matching force between spins.
       
    Parameters:
        number_of_particles : number of particles present in the system of
                              particles.
        random_or_not : integer that determines if the coefficients are
                        generated randomly or are given by the user
                        (1=random, 0=not random).
    
    Returns:
        The Jij coefficients in form of a matrix of dimension
        number_of_particles x number_of_particles.
    
    Raise:
        ValueError if:
            The dimension of the matrix is not consistent with number_of_particles.
            The diagonal of the matrix is not comprised of 0s.
            The matrix is not symmetric
    """
    if random_or_not == 0:
        with open("Simulation_parameters.txt", "r") as parameters:
            Jij = np.zeros((number_of_particles, number_of_particles))
            lines = parameters.read().splitlines()
            i = 0
            line_marker = 1
            for line in lines:
                if line_marker > 3:
                    numbers = line.split(",")
                    if len(numbers) != number_of_particles:
                        raise ValueError(
                            "The dimension of the matrix must be consistent with the number of particles.")
                    j = 0
                    for number in numbers:
                        num = float(number)
                        Jij[i, j] = num
                        if i == j:
                            if Jij[i, i] != 0.:
                                raise ValueError(
                                    "The diagonal of the J_ij matrix must contain 0s. Check the parameters file.")
                        if i > j:
                            if Jij[i, j] != Jij[j, i]:
                                raise ValueError(
                                    "The matrix of J_ij coefficients should be symmetric. Check the parameters file.")
                        j += 1
                    i += 1
                line_marker += 1
        return(Jij)
    elif random_or_not == 1:
        Jij = np.random.rand(number_of_particles, number_of_particles)
        i = 0
        while i < number_of_particles:
            j = 0
            while j < number_of_particles:
                if i == j:  
                    # The coupling coefficient is 0 if the index is equal
                    Jij[i, j] = 0.
                if i < j:
                    Jij[j, i] = Jij[i, j]
                j = j + 1
            i = i + 1
        return(Jij)


def btest(i, n):
    """This function determines if the n-th bit of a number is 0 or 1.
    
    Parameters: 
        i : number.
        n : n-th bit of the number i.
        
    Returns: 
        +1 if the n-th bit of i is 1.
        -1 if the n-th bit of i is 0.
    """
    if (i & (1 << n)):
        return(1.)
        # n-th bit is set (1)

    else:
        return(-1.)
        # n-th bit is not set (0)


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
    ham1 = np.zeros((2**number_of_particles, 2**number_of_particles))
    i = 0
    while i < 2**(number_of_particles):
        j = 0
        suma = 0.
        while j < number_of_particles:
            k = 0
            while k < number_of_particles:
                suma = suma + btest(i, j) * btest(i, k) * Jij[j, k]
                k = k + 1
            suma = suma + btest(i, j) * hi[j]
            j = j + 1
        ham1[i, i] = suma
        i = i + 1
    return(ham1)

# Create the initial Hamiltonian representing the transversal field


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
    H0 = np.zeros((2**number_of_particles, 2**number_of_particles))
    i = 0
    while i < 2**number_of_particles:
        j = 0
        while j < number_of_particles:
            if btest(i, j) == 1.:
                H0[i, i - 2**j] = 1
                H0[i - 2**j, i] = 1
            else:
                H0[i, i + 2**j] = 1
                H0[i + 2**j, i] = 1
            j = j + 1
        i = i + 1
    H0=-H0
    return(H0)


def gap_two_smallest_in_array(array):
    """This function calculates the difference between the two smallest numbers
       in an array.
    
    Parameters:
        array : array of numbers.
        
    Returns:
        The difference between the two smallest numbers in array.
    """
    array_size = len(array)
    first = second = float("inf")
    for i in range(0, array_size):
        # If current element is smaller than first then
        # update both first and second
        if array[i] < first:
            second = first
            first = array[i]
        # If arr[i] is in between first and second then
        # update second
        elif (array[i] < second and array[i] != first):
            second = array[i]
    gap = second - first
    return(gap)


def position_smallest_in_array(array):
    """This function calculates the position of the smallest number in an array.
    
    Parameters:
        array : array of numbers.
        
    Returns:
        The position of the smallest number in array.
    """
    
    array_size = len(array)
    smallest = float("inf")
    for i in range(0, array_size):
        # If current element is smaller than smallest, update smallest
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
    step = 1.0 / (float(num_steps))  # Size of the step
    s = 0.0 # Adiabatic parameter that evolves from 0 to 1.
    i = 0
    gap = np.zeros(num_steps + 1)  # array that will contain energy gap
    eigenvectors = np.zeros((num_steps + 1, 2**number_of_particles))
    while i <= num_steps:
        # We construct the matrix H(s) which evolves from H0 to ham1
        H_s = s * H1 + (1.0 - s) * H0
        # We diagonalize H_s and obtain its eigenvalues and eigenvectors
        eigenvals, eigenvecs = la.eig(H_s)
        eigenvals = eigenvals.real  # converts eigenvalues to real numbers
        # We obtain the two lowesr eigenvalues, which will give the ground and
        # the first excited state
        gap[i] = gap_two_smallest_in_array(eigenvals)
        position = position_smallest_in_array(eigenvals)
        # The eigenvectors are the columns of eigenvecs
        eigenvectors[i] = np.array(eigenvecs[:, position])
        s = s + step
        i = i + 1
    probability = eigenvectors**2
    return(gap, probability)

# Let's integrate to obtain the time required for the AQC


def results_of_simulation(num_steps, number_of_particles, H0, H1):
    """ This function calculates the time required for the AQC and plots the 
        evolution of the Gap, the speed and the probability of the states.
        
    Parameters: 
        num_steps : number of steps when evolving from H0 to H1.
        number_of_particles : number of particles present in the system of
                              particles.
        H0 : Initial Hamiltonian.
        H1 : Target Hamiltonian, towards which the system evolves.
        
    Returns:
        The time required for the Adiabatic Quantum Computation.
    """
    #Calculating the time
    gap, probability = quantum_simulation(num_steps, number_of_particles, H0, H1)
    xaxis = np.linspace(0., 1., 101)
    speed = gap**2
    speed_inverse = 1.0 / (speed)
    time = integrate.simps(speed_inverse, xaxis)
    time = round(time, 4)
    
    #Plotting of the Gap and Speed
    figure = plt.figure()
    yaxis = np.array([gap, speed])
    label = np.array(["Gap", "Speed"])
    j = 0
    for i in yaxis:
        plt.plot(xaxis, i, label=label[j])
        j = j + 1
    plt.ylabel("Gap, Speed")
    plt.xlabel("Adiabatic parameter (s)")
    plt.legend()
    figure.savefig("Graphs/Gap_and_Speed.png")
    
    #Plotting of the probability of the states
    figure = plt.figure()
    i = 0
    for i in range(2**number_of_particles):
        plt.plot(xaxis, probability[:, i])
    plt.yscale('log')
    plt.gca().set_ylim([0.0001, 1.])
    plt.ylabel("Probability")
    plt.xlabel("Adiabatic parameter (s)")
    figure.savefig('Graphs/States.png')
    
    return(time)

# Once all functions are written initialize the variables and call the
# functions


# Open the file to read the parameters
parameters = open("Simulation_parameters.txt", "r")
all_lines = parameters.readlines()
parameters.close()

# Assign the values to the functions' input variables
number_of_particles = int(all_lines[0])
random_or_not = int(all_lines[1])
hi = coefficient_hi_generator(number_of_particles, random_or_not)
Jij = coefficient_Jij_generator(number_of_particles, random_or_not)
H0 = Hamiltonian_0(int(all_lines[0]))
H1 = Hamiltonian_1(int(all_lines[0]), hi, Jij)

# Call the functions to return the desired plots and results
print("The time required for the Adiabatic Quantum Computing is: ",
      results_of_simulation(100, number_of_particles, H0, H1))
