import random
import numpy as np
import scipy.linalg as la
import scipy.integrate as integrate
import matplotlib.pyplot as plt


parameters=open("Simulation_parameters.txt","r")
all_lines=parameters.readlines()
parameters.close()
#nspin determines the number of particles used in the simulation
nspin=int(all_lines[0])
