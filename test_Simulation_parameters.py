import Simulation
import numpy
import pytest

#Test the hi coefficient generator function to see if the return value is a properly shaped numpy array of real numbers.
def test_coefficient_hi_generator():
    assert isinstance(Simulation.coefficient_hi_generator(3, 1), numpy.ndarray)
    assert numpy.shape(Simulation.coefficient_hi_generator(3, 1))==(3,)
    i=0
    while i<3:
        assert isinstance(Simulation.coefficient_hi_generator(3, 1)[i], float)
        i+=1