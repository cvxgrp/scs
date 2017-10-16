from __future__ import print_function
import platform
## import utilities to generate random cone probs:
import sys
sys.path.insert(0, '../examples/python')
from genRandomScsConeProb import *


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

def check_solution(solution, expected):
  assert_almost_equals(solution, expected, places=2)

def assert_(str1, str2):
  if (str1 != str2):
    print("assert failure: %s != %s" % (str1, str2))
  assert str1 == str2

def check_infeasible(sol):
  assert_(sol['info']['status'], 'Infeasible')

def check_unbounded(sol):
  assert_(sol['info']['status'], 'Unbounded')

random.seed(0)
num_feas = 50
num_unb = 10
num_infeas = 10

opts={'max_iters':100000,'eps':1e-5} # better accuracy than default to ensure test pass
K = {'f':10, 'l':25, 'q':[5, 10, 0 ,1], 's':[2, 1, 2, 0, 1, 10, 8], 'ep':0, 'ed':0, 'p':[0.25, -0.75, 0.33, -0.33, 0.2]}
m = getScsConeDims(K)

def test_feasible():
    for i in range(num_feas):
        data, p_star = genFeasible(K, n = m // 3, density = 0.1)
        
        sol = scs.solve(data, K, use_indirect=False, **opts)
        yield check_solution, dot(data['c'],sol['x']), p_star
        yield check_solution, dot(-data['b'],sol['y']), p_star

        sol = scs.solve(data, K, use_indirect=True, **opts)
        yield check_solution, dot(data['c'],sol['x']), p_star
        yield check_solution, dot(-data['b'],sol['y']), p_star

def test_infeasible():
    for i in range(num_infeas):
        data = genInfeasible(K, n = m // 3)
        
        yield check_infeasible, scs.solve(data, K, use_indirect=False, **opts)
        yield check_infeasible, scs.solve(data, K, use_indirect=True, **opts)

def test_unbounded():
    for i in range(num_unb):
        data = genUnbounded(K, n = m // 2)
        
        yield check_unbounded, scs.solve(data, K, use_indirect=False, **opts)
        yield check_unbounded, scs.solve(data, K, use_indirect=True, **opts)

