import Simulation
import configparser


# Open the file to read the parameters.
initial_parameters = configparser.ConfigParser()
initial_parameters.read("Simulation_parameters.ini")


# Assign the values to the functions' input variables.
number_of_particles = int(
    initial_parameters.get(
        "Parameters",
        "Number of Particles"))
random_or_not = int(initial_parameters.get("Parameters", "Random selection"))
if random_or_not == 0:
    # The coeffcicients are the ones given by the user
    hi_coefficients = initial_parameters.get("Coefficients", "hi_coefficients")
    Jij_coefficients = initial_parameters.get(
        "Coefficients", "Jij_coefficients")
    hi = Simulation.coefficient_hi_reader(number_of_particles, hi_coefficients)
    Jij = Simulation.coefficient_Jij_reader(
        number_of_particles, Jij_coefficients)
elif random_or_not == 1:
    # The coefficients are randomly generated
    hi = Simulation.random_coefficient_hi_generator(number_of_particles)
    Jij = Simulation.random_coefficient_Jij_generator(number_of_particles)
H0 = Simulation.Hamiltonian_0(number_of_particles)
H1 = Simulation.Hamiltonian_1(number_of_particles, hi, Jij)

# Call the functions to return the desired plots and results.
print("The time required for the Adiabatic Quantum Computing is: ",
      Simulation.results_of_simulation(100, number_of_particles, H0, H1))
