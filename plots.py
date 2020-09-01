import numpy as np
import matplotlib.pyplot as plt
import configparser


def plot_gap_and_speed():
    """ This function plots the evolution of the Gap and Speed during the
        simulation.

    Parameters:
        gap : evolution of the gap.
        speed : evolution of the speed.
        xaxis : evolution of the adiabatic parameter, which will be the x axis
                of the plot.
        yaxis : the gap and the speed as y axis of the plot.
        figure : Figure to plot
    """
    gap = np.loadtxt("Data/Gap.txt").reshape(101,)
    speed = np.loadtxt("Data/Speed.txt").reshape(101,)
    xaxis = np.linspace(0., 1., 101)
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


def plot_probability(number_of_particles):
    """ This function plots the evolution of the probability of the different
        states during the simulation.

    Parameters:
        probability : evolution of the probability.
        xaxis : evolution of the adiabatic parameter, which will be the x axis
                of the plot.
        figure : Figure to plot
    """
    probability = np.loadtxt(
        "Data/States.txt").reshape(101, 2**number_of_particles)
    xaxis = np.linspace(0., 1., 101)
    figure = plt.figure()
    i = 0
    for i in range(2**number_of_particles):
        plt.plot(xaxis, probability[:, i])
    plt.yscale('log')
    plt.gca().set_ylim([0.0001, 1.])
    plt.ylabel("Probability")
    plt.xlabel("Adiabatic parameter (s)")
    figure.savefig('Graphs/States.png')


# Read the parameters contained in the parameters file
initial_parameters = configparser.ConfigParser()
initial_parameters.read("Simulation_parameters.ini")

# Assign value to the number of particles
number_of_particles = number_of_particles = int(
    initial_parameters.get(
        "Parameters",
        "Number of Particles"))

# Call the plots
plot_gap_and_speed()
plot_probability(number_of_particles)
