import random
import sys
import numpy as np
import scipy.linalg as la
import scipy.integrate as integrate
import matplotlib.pyplot as plt

parameters=open("Simulation_parameters.txt","r")
all_lines=parameters.readlines()
parameters.close()

#This function generates hi coefficients in order to generate the later Hamiltonians
def coefficient_hi_generator(number_of_particles,random_or_not):
    if abs(number_of_particles)>8 or abs(number_of_particles)==0:
        sys.exit("You must insert as number of particles a positive integer lower than 8 and higher than 0.")
    elif random_or_not not in [0,1]:
        sys.exit("The second line of the file containing the parameters must be 0 or 1.")
    else:
        if random_or_not==0:
            #Create hi array with 0s to give it value later
            hi=np.zeros(number_of_particles)
            i=0
            while i<number_of_particles:
                hi[i]=float(input("Introduce the element h{}: ".format(i+1)))
                i=i+1
            i=0
            return(hi)
        elif random_or_not==1:
            hi=np.random.rand(number_of_particles)
            return(hi)  

#This function generates Jij coefficients in order to generate the later Hamiltonians       
def coefficient_Jij_generator(number_of_particles,random_or_not):
    if abs(number_of_particles)>8 or abs(number_of_particles)==0:
        sys.exit("You must insert as number of particles a positive integer lower than 8 and higher than 0.")
    elif random_or_not not in [0,1]:
        sys.exit("The second line of the file containing the parameters must be 0 or 1.")
    else:
        if random_or_not==0:
            #Create Jij array with 0s to give it value later
            Jij=np.zeros((number_of_particles,number_of_particles))
            i=0
            while i<number_of_particles:
                j=0
                while j<number_of_particles:
                    if i<j:
                        Jij[i,j]=float(input("Introduce the element J{}{}: ".format(i+1,j+1)))
                        Jij[j,i]=Jij[i,j]
                    j=j+1
                i=i+1
            return(Jij)
        elif random_or_not==1:
            Jij=np.random.rand(number_of_particles,number_of_particles)
            i=0
            while i<number_of_particles:
                j=0
                while j<number_of_particles:
                    if i==j:#The coupling coefficient is 0 if the index is equal
                        Jij[i,j]=0.
                    if i<j:
                        Jij[j,i]=Jij[i,j]
                    j=j+1 
                i=i+1
            return(Jij)    
        
#Determines if the n-th bit of a number is 0 or 1 and returns +1 or -1. Useful for the construction of the Hamiltonians.
def btest(i,n):
    if (i & (1<<n)) :
        return(1.)
        ## n-th bit is set (1)

    else:
        return(-1.)
        ## n-th bit is not set (0)
    
#Create the target Hamiltonian towards which the system will evolve
def Hamiltonian_1(number_of_particles,hi=coefficient_hi_generator(int(all_lines[0]), int(all_lines[1])),Jij=coefficient_Jij_generator(int(all_lines[0]), int(all_lines[1]))):
    ham1=np.zeros((2**number_of_particles,2**number_of_particles))
    i=0
    while i<2**(number_of_particles):
        j=0
        suma=0.
        while j<number_of_particles:
            k=0
            while k<number_of_particles:
                suma=suma+btest(i,j)*btest(i,k)*Jij[j,k]
                k=k+1
            suma=suma+btest(i,j)*hi[j]
            j=j+1
        ham1[i,i]=suma
        i=i+1
    return(ham1)

