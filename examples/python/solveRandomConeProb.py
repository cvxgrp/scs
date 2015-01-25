from __future__ import print_function
import scs
from numpy import *
from scipy import sparse, randn
from genRandomConeProb import *

#############################################
#  Uses scs to solve a random cone problem  #
#############################################

def main():
    #random.seed(0)
    solveFeasible()
    solveInfeasible()
    solveUnbounded()
    
def solveFeasible():
    # cone:
    K = {'f':10, 'l':15, 'q':[5, 10, 0 ,1], 's':[3, 4, 0, 0, 1], 'ep':10, 'ed':10, 'p':[-0.25, 0.5, 0.75, -0.33]}
    m = getConeDims(K) 
    data, p_star = genFeasible(K, n = m/3, density = 0.01)
    params = {'eps':1e-3, 'normalize':True, 'scale':5, 'cg_rate':2}
    
    sol_i = scs.solve(data, K, use_indirect=True, **params)
    xi = sol_i['x']
    yi = sol_i['y']
    print('p*  = ', p_star)
    print('pri error = ', (dot(data['c'], xi) - p_star) / p_star)
    print('dual error = ', (-dot(data['b'], yi) - p_star) / p_star)
    # direct:
    sol_d = scs.solve(data, K, **params)
    xd = sol_d['x']
    yd = sol_d['y']
    print('p*  = ', p_star)
    print('pri error = ', (dot(data['c'], xd) - p_star) / p_star)
    print('dual error = ', (-dot(data['b'], yd) - p_star) / p_star)
    
def solveInfeasible():
    K = {'f':10, 'l':15, 'q':[5, 10, 0 ,1], 's':[3, 4, 0, 0, 1], 'ep':10, 'ed':10, 'p':[-0.25, 0.5, 0.75, -0.33]}
    m = getConeDims(K)
    data = genInfeasible(K, n = m/3)
    params = {'eps':1e-4, 'normalize':True, 'scale':0.5, 'cg_rate':2}
    sol_i = scs.solve(data, K, use_indirect=True, **params)
    sol_d = scs.solve(data, K, **params)

def solveUnbounded():
    K = {'f':10, 'l':15, 'q':[5, 10, 0 ,1], 's':[3, 4, 0, 0, 1], 'ep':10, 'ed':10, 'p':[-0.25, 0.5, 0.75, -0.33]}
    m = getConeDims(K)
    data = genUnbounded(K, n = m/3)
    params = {'eps':1e-4, 'normalize':True, 'scale':0.5, 'cg_rate':2}
    sol_i = scs.solve(data, K, use_indirect=True, **params)
    sol_d = scs.solve(data, K, **params)

if __name__ == "__main__":
   main()
