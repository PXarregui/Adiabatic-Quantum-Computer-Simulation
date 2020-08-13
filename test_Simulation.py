import Simulation
import numpy
import pytest

#Test the hi coefficient generator function to see if the return value is a properly shaped numpy array of real numbers.
def test_coefficient_hi_generator():
    assert numpy.shape(Simulation.coefficient_hi_generator(3, 1))==(3,)
    i=0
    while i<3:
        assert isinstance(Simulation.coefficient_hi_generator(3, 1)[i], float)
        i+=1

#Test if the Jij coefficient matrix is symmetric with 0s in the diagonal and real numbers elsewhere.
def test_matrix_coefficient_Jij_generator():
    i=0
    generated_matrix=Simulation.coefficient_Jij_generator(3, 1)
    while i<3:
        j=0
        while j<3:
            if i==j:
                assert generated_matrix[i,j]==0.
            else:
                assert isinstance(generated_matrix[i,j], float)
                assert generated_matrix[i,j]==generated_matrix[j,i]
            j+=1
        i+=1
        