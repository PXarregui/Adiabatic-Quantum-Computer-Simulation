import Simulation
import numpy
import pytest

# Tests concerning coefficient_hi_generator()

#Test the hi coefficient generator function to see if the return value is a properly shaped numpy array.
def test_shape_coefficient_hi_generator():
    assert numpy.shape(Simulation.coefficient_hi_generator(3, 1))==(3,)
 
#Test the hi coefficient generator function to see if the return value is an array of real numbers.
def test_real_coefficient_hi_generator():
    i=0
    while i<3:
        assert isinstance(Simulation.coefficient_hi_generator(3, 1)[i], float)
        i+=1

# -----------------------------------------------------------------------------
        
# Tests concerning coefficient_Jij_generator()

#Test the Jij generator function to see if the matrix has the apropriate dimension. 
def test_shape_coefficient_Jij_generator():
    assert numpy.shape(Simulation.coefficient_Jij_generator(3, 1))==(3,3)

#Test if the Jij coefficient matrix has 0s in the diagonal.
def test_matrix_diagonal_coefficient_Jij_generator():
    i=0
    generated_matrix=Simulation.coefficient_Jij_generator(3, 1)
    while i<3:
        j=0
        while j<3:
            if i==j:
                assert generated_matrix[i,j]==0.
            j+=1
        i+=1
         
#Test if the Jij coefficient matrix is composed of real numbers elsewhere.
def test_matrix_real_coefficient_Jij_generator():
    i=0
    generated_matrix=Simulation.coefficient_Jij_generator(3, 1)
    while i<3:
        j=0
        while j<3:
            if i!=j:
                assert isinstance(generated_matrix[i,j], float)
            j+=1
        i+=1
 
#Test if the Jij coefficient matrix is symmetric.
def test_symmetry_coefficient_Jij_generator():
    i=0
    generated_matrix=Simulation.coefficient_Jij_generator(3, 1)
    while i<3:
        j=0
        while j<3:
            if i!=j:
                assert generated_matrix[i,j]==generated_matrix[j,i]
            j+=1
        i+=1      
        
# -----------------------------------------------------------------------------        
  
#Test btest function, to see if it returns 1.0 or -1.0 when it should.
def test_btest_1():
    assert Simulation.btest(100,2)==1.0

def test_btest_2():
    assert Simulation.btest(99,2)==-1,0

# -----------------------------------------------------------------------------
    
def test_Hamiltonian_1():
    hi=numpy.array((1,1))
    Jij=numpy.array(((0,1),(1,0)))
    expected_result=numpy.array(((0.,0.,0.,0.),(0.,2.,0.,0.),(0.,0.,2.,0.),(0.,0.,0.,4.)))
    actual_result=Simulation.Hamiltonian_1(2,hi,Jij)
    assert actual_result.all()==expected_result.all()
    
# -----------------------------------------------------------------------------
        
def test_Hamiltonian_0():
    expected_result=numpy.array(((0.,1.,1.,0.),(1.,0.,0.,1.),(1.,0.,0.,1.),(0.,1.,1.,0.)))
    actual_result=Simulation.Hamiltonian_0(2)
    assert actual_result.all()==expected_result.all()       