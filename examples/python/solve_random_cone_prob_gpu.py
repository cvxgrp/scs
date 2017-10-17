from __future__ import print_function, division
import scs
from numpy import *
from scipy import sparse, randn
from gen_random_cone_prob import *

#############################################
#  Uses scs to solve a random cone problem  #
#############################################


def main():
  random.seed(1)
  solve_feasible()
  solve_infeasible()
  solve_unbounded()


def solve_feasible():
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
  m = get_scs_cone_dims(K)
  data, p_star = gen_feasible(K, n=m // 3, density=0.01)
  params = {'eps': 1e-3, 'normalize': True, 'scale': 5, 'cg_rate': 2}

  sol_i = scs.solve(data, K, gpu=True, **params)
  xi = sol_i['x']
  yi = sol_i['y']
  print('p*  = ', p_star)
  print('pri error = ', (dot(data['c'], xi) - p_star) / p_star)
  print('dual error = ', (-dot(data['b'], yi) - p_star) / p_star)


def solve_infeasible():
  K = {
      'f': 10,
      'l': 15,
      'q': [5, 10, 0, 1],
      's': [3, 4, 0, 0, 1],
      'ep': 10,
      'ed': 10,
      'p': [-0.25, 0.5, 0.75, -0.33]
  }
  m = get_scs_cone_dims(K)
  data = gen_infeasible(K, n=m // 3)
  params = {'eps': 1e-4, 'normalize': True, 'scale': 0.5, 'cg_rate': 2}
  sol_i = scs.solve(data, K, gpu=True, **params)


def solve_unbounded():
  K = {
      'f': 10,
      'l': 15,
      'q': [5, 10, 0, 1],
      's': [3, 4, 0, 0, 1],
      'ep': 10,
      'ed': 10,
      'p': [-0.25, 0.5, 0.75, -0.33]
  }
  m = get_scs_cone_dims(K)
  data = gen_unbounded(K, n=m // 3)
  params = {'eps': 1e-4, 'normalize': True, 'scale': 0.5, 'cg_rate': 2}
  sol_i = scs.solve(data, K, gpu=True, **params)


if __name__ == '__main__':
  main()
