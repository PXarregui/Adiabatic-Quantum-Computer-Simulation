import sys
import numpy as np
import scipy.linalg as la
import scipy.integrate as integrate
import matplotlib.pyplot as plt

#This function generates hi coefficients in order to generate the later Hamiltonians
def coefficient_hi_generator(number_of_particles,random_or_not):
    if abs(number_of_particles)>8 or abs(number_of_particles)<2:
        sys.exit("You must insert as number of particles a positive integer lower or equal to 8 and higher than 1.")
    elif random_or_not not in [0,1]:
        sys.exit("The second line of the file containing the parameters must be 0 or 1.")
    else:
        if random_or_not==0:
            with open("Simulation_parameters.txt","r") as parameters:
                hi=np.zeros(number_of_particles)
                lines=parameters.read().splitlines()
                line_marker=1
                for line in lines:
                    if line_marker==3:
                        numbers=line.split(",")
                        if len(numbers)!=number_of_particles:
                            sys.exit("The dimension of the vector must be equal to the number of particles.")
                        k=0
                        for number in numbers:
                            num=float(number)
                            hi[k]=num
                            k+=1
                    line_marker+=1
            return(hi)
        elif random_or_not==1:
            hi=np.random.rand(number_of_particles)
            return(hi)  

#This function generates Jij coefficients in order to generate the later Hamiltonians       
def coefficient_Jij_generator(number_of_particles,random_or_not):
    if abs(number_of_particles)>8 or abs(number_of_particles)<2:
        sys.exit("You must insert as number of particles a positive integer lower or equal to 8 and higher than 1.")
    elif random_or_not not in [0,1]:
        sys.exit("The second line of the file containing the parameters must be 0 or 1.")
    else:
        if random_or_not==0:
            with open("Simulation_parameters.txt","r") as parameters:
                Jij=np.zeros((number_of_particles,number_of_particles))
                lines=parameters.read().splitlines()
                i=0
                line_marker=1
                for line in lines:
                    if line_marker>3:
                        numbers=line.split(",")
                        if len(numbers)!=number_of_particles:
                            sys.exit("The dimension of the matrix must be consistent with the number of particles.")
                        j=0
                        for number in numbers:
                            num=float(number)
                            Jij[i,j]=num
                            if i==j:
                                if Jij[i,i]!=0.:
                                    sys.exit("The diagonal of the J_ij matrix must contain 0s. Check the parameters file.")
                            if i>j:
                                if Jij[i,j]!=Jij[j,i]:
                                    sys.exit("The matrix of J_ij coefficients should be symmetric. Check the parameters file.")
                            j+=1
                        i+=1
                    line_marker+=1 
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
def Hamiltonian_1(number_of_particles,hi,Jij):
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

#Create the initial Hamiltonian representing the transversal field
def Hamiltonian_0(number_of_particles):
    H0=np.zeros((2**number_of_particles,2**number_of_particles))
    i=0
    while i<2**number_of_particles:
        j=0
        while j<number_of_particles:
            if btest(i,j)==1.:
                H0[i,i-2**j]=1
                H0[i-2**j,i]=1
            else:
                H0[i,i+2**j]=1
                H0[i+2**j,i]=1
            j=j+1
        i=i+1
    return(H0)

def gap_two_smallest_in_array(array):
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
            second = array[i]; 
    gap=second-first
    return(gap) 

def position_smallest_in_array(array):
    array_size = len(array) 
    smallest= float("inf")
    for i in range(0, array_size): 
        # If current element is smaller than smallest, update smallest
        if array[i] < smallest: 
            smallest = array[i] 
            position=i
    return(position)

def quantum_simulation(num_steps,number_of_particles,H0,H1):
    step=1.0/(float(num_steps)) #Size of the step
    s=0.0
    i=0
    gap=np.zeros(num_steps+1)#array that will contain energy gap
    eigenvectors=np.zeros((num_steps+1,2**number_of_particles))
    while i<=num_steps:
        #We construct the matrix H(s) which evolves from H0 to ham1
        H_s=s*H1+(1.0-s)*H0
        #We diagonalize H_s and obtain its eigenvalues and eigenvectors
        eigenvals,eigenvecs=la.eig(H_s)
        eigenvals=eigenvals.real #converts eigenvalues to real numbers
        #We obtain the two lowesr eigenvalues, which will give the ground and the first excited state
        gap[i]=gap_two_smallest_in_array(eigenvals)
        position=position_smallest_in_array(eigenvals)
        eigenvectors[i]=np.array(eigenvecs[:,position])#The eigenvectors are the columns of eigenvecs
        s=s+step
        i=i+1
    probability=eigenvectors**2
    return(gap,probability)

#We can represent the evolution of the gap and speed   
def plot_gap(num_steps,number_of_particles,H0,H1):
    gap,_=quantum_simulation(num_steps,number_of_particles,H0,H1)
    speed=gap**2
    xaxis=np.linspace(0.,1.,101)
    yaxis=np.array([gap,speed])
    plt.figure()
    label=np.array(["Gap","Speed"])
    j=0
    for i in yaxis:
        plt.plot(xaxis,i,label=label[j])
        j=j+1
    plt.ylabel("Gap, Speed")
    plt.xlabel("Adiabatic parameter (s)")
    plt.legend()
    plt.show()
    
#We also can represent the evolution of the probability of each state 
def plot_states(num_steps,number_of_particles,H0,H1):
    _,probability=quantum_simulation(num_steps,number_of_particles,H0,H1)
    xaxis=np.linspace(0.,1.,101)
    plt.figure()
    i=0
    for i in range(2**number_of_particles):
        plt.plot(xaxis,probability[:,i])
    plt.yscale('log')
    plt.gca().set_ylim([0.0001,1.])
    plt.ylabel("Probability")
    plt.xlabel("Adiabatic parameter (s)")
    plt.show()

#Let's integrate to obtain the time required for the AQC
def computation_time(num_steps,number_of_particles,H0,H1):
    gap,_=quantum_simulation(num_steps,number_of_particles,H0,H1)
    xaxis=np.linspace(0.,1.,101)
    speed_inverse=1.0/(gap**2)
    time=integrate.simps(speed_inverse,xaxis)
    time=round(time,4)
    return(time)

# Once all functions are written initialize the variables and call the functions
    
#Open the file to read the parameters
parameters=open("Simulation_parameters.txt","r")
all_lines=parameters.readlines()
parameters.close()

#Assign the values to the functions' input variables
number_of_particles=int(all_lines[0])
random_or_not=int(all_lines[1])
hi=coefficient_hi_generator(number_of_particles,random_or_not)
Jij=coefficient_Jij_generator(number_of_particles,random_or_not)
H0=Hamiltonian_0(int(all_lines[0]))
H1=Hamiltonian_1(int(all_lines[0]),hi,Jij)

#Call the functions to return the desired plots and results
plot_gap(100,number_of_particles,H0,H1)
plot_states(100,number_of_particles,H0,H1)
print("The time required for the Adiabatic quantum computing is: ", computation_time(100,number_of_particles,H0,H1))