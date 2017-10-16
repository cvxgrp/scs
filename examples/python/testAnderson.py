from __future__ import print_function, division
import scs
from numpy import *
from scipy import sparse, randn
from genRandomScsConeProb import *

#############################################
#  Uses scs to solve a random cone problem  #
#############################################


def main():
  #random.seed(0)
  solveFeasible()


def solveFeasible():
  # cone:
  K = {
      'f': 10,
      'l': 15,
      'q': [5, 10, 0, 1],
      's': [3, 4, 0, 0, 1],
      'ep': 10,
      'ed': 10,
      'p': [-0.25, 0.5, 0.75, -0.33]
  }
  m = getScsConeDims(K)
  data, p_star = genFeasible(K, n=m // 3, density=0.01)
  params = {
      'normalize': True,
      'scale': 5,
      'cg_rate': 2,
      'acceleration_lookback': 10
  }

  sol_i = scs.solve(data, K, use_indirect=True, **params)
  xi = sol_i['x']
  yi = sol_i['y']
  print('p*  = ', p_star)
  print('pri error = ', (dot(data['c'], xi) - p_star) / p_star)
  print('dual error = ', (-dot(data['b'], yi) - p_star) / p_star)
  # direct:
  sol_d = scs.solve(data, K, use_indirect=False, **params)
  xd = sol_d['x']
  yd = sol_d['y']
  print('p*  = ', p_star)
  print('pri error = ', (dot(data['c'], xd) - p_star) / p_star)
  print('dual error = ', (-dot(data['b'], yd) - p_star) / p_star)

  params['acceleration_lookback'] = 30
  sol_i = scs.solve(data, K, use_indirect=True, **params)
  xi = sol_i['x']
  yi = sol_i['y']
  print('p*  = ', p_star)
  print('pri error = ', (dot(data['c'], xi) - p_star) / p_star)
  print('dual error = ', (-dot(data['b'], yi) - p_star) / p_star)

  params['acceleration_lookback'] = 0
  sol_i = scs.solve(data, K, use_indirect=True, **params)
  xi = sol_i['x']
  yi = sol_i['y']
  print('p*  = ', p_star)
  print('pri error = ', (dot(data['c'], xi) - p_star) / p_star)
  print('dual error = ', (-dot(data['b'], yi) - p_star) / p_star)


if __name__ == '__main__':
  main()
