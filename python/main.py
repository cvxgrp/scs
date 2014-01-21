import scs_direct as p
import scs_indirect as q
from sys import getrefcount
import numpy

#from guppy import hpy

Ax = numpy.matrix([1.,1.])
Ai = numpy.matrix([0,1])
Ap = numpy.matrix([0,1,2])
b = numpy.matrix([1.,1.])
c = numpy.matrix([1.,1.])

solution = p.solve(Ax, Ai, Ap, b, c)
print solution['x']
print solution['y']
print solution['status']

print getrefcount(solution['x'])

# h = hpy()
# print h.heap()