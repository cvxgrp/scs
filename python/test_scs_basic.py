# nost test suite copied initially from ECOS project
from __future__ import print_function
import platform

def import_error(msg):
  print()
  print("## IMPORT ERROR:" + msg)
  print()

try:
  from nose.tools import assert_raises, assert_almost_equals
except ImportError:
  import_error("Please install nose to run tests.")
  raise

try:
  import scs
except ImportError:
  import_error("You must install the scs module before running tests.")
  raise

try:
  import numpy as np
except ImportError:
  import_error("Please install numpy.")
  raise

try:
  import scipy.sparse as sp
except ImportError:
  import_error("Please install scipy.")
  raise

# global data structures for problem
c = np.array([-1.])
b = np.array([1., -0.])
A = sp.csc_matrix([1., -1.]).T.tocsc()
data = {'A':A, 'b':b, 'c':c}
cone = {'q': [], 'l': 2}

FAIL = 'Failure' # scs code for failure

def check_solution(solution, expected):
  assert_almost_equals(solution, expected, places=2)

def check_failure(sol):
  assert sol['info']['status'] == FAIL

def test_problems():
  sol = scs.solve(data, cone, use_indirect=False)
  yield check_solution, sol['x'][0], 1

  new_cone = {'q':[2], 'l': 0}
  sol = scs.solve(data, new_cone, use_indirect=False)
  yield check_solution, sol['x'][0], 0.5

  sol = scs.solve(data, cone, use_indirect = True )
  yield check_solution, sol['x'][0], 1

  sol = scs.solve(data, new_cone, use_indirect = True )
  yield check_solution, sol['x'][0], 0.5


if platform.python_version_tuple() < ('3','0','0'):
  def test_problems_with_longs():
    new_cone = {'q': [], 'l': long(2)}
    sol = scs.solve(data, new_cone, use_indirect=False)
    yield check_solution, sol['x'][0], 1
    sol = scs.solve(data, new_cone, use_indirect=True )
    yield check_solution, sol['x'][0], 1

    new_cone = {'q':[long(2)], 'l': 0}
    sol = scs.solve(data, new_cone, use_indirect=False)
    yield check_solution, sol['x'][0], 0.5
    sol = scs.solve(data, new_cone, use_indirect=True )
    yield check_solution, sol['x'][0], 0.5

def check_keyword(error_type, keyword, value):
  assert_raises(error_type, scs.solve, data, cone, **{keyword: value})

def test_failures():
  yield assert_raises, TypeError, scs.solve
  yield assert_raises, ValueError, scs.solve, data, {'q':[4], 'l':-2}
  yield check_keyword, ValueError, 'max_iters', -1
  # python 2.6 and before just cast float to int
  if platform.python_version_tuple() >= ('2', '7', '0'):
    yield check_keyword, TypeError, 'max_iters', 1.1

  yield check_failure, scs.solve( data, {'q':[1], 'l': 0} )

