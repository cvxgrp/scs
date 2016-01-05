from __future__ import print_function, division
import scs
from numpy import *
from scipy import sparse, randn
from genRandomConeProb import *

#############################################
#  Uses scs to solve a random cone problem  #
#############################################

def main():
    random.seed(1)
    solveFeasible()
    solveInfeasible()
    solveUnbounded()
    
def solveFeasible():
    # cone:
    K = {'f':10, 'l':15, 'q':[5, 10, 0 ,1], 's':[3, 4, 0, 0, 1], 'ep':10, 'ed':10, 'p':[-0.25, 0.5, 0.75, -0.33]}
    m = getConeDims(K) 
    data, p_star = genFeasible(K, n = m//3, density = 0.01)
    params = {'eps':1e-3, 'normalize':True, 'scale':5, 'cg_rate':2}
    
    sol_i = scs.solve(data, K, gpu=True, **params)
    xi = sol_i['x']
    yi = sol_i['y']
    print('p*  = ', p_star)
    print('pri error = ', (dot(data['c'], xi) - p_star) / p_star)
    print('dual error = ', (-dot(data['b'], yi) - p_star) / p_star)
    
def solveInfeasible():
    K = {'f':10, 'l':15, 'q':[5, 10, 0 ,1], 's':[3, 4, 0, 0, 1], 'ep':10, 'ed':10, 'p':[-0.25, 0.5, 0.75, -0.33]}
    m = getConeDims(K)
    data = genInfeasible(K, n = m//3)
    params = {'eps':1e-4, 'normalize':True, 'scale':0.5, 'cg_rate':2}
    sol_i = scs.solve(data, K, gpu=True, **params)

def solveUnbounded():
    K = {'f':10, 'l':15, 'q':[5, 10, 0 ,1], 's':[3, 4, 0, 0, 1], 'ep':10, 'ed':10, 'p':[-0.25, 0.5, 0.75, -0.33]}
    m = getConeDims(K)
    data = genUnbounded(K, n = m//3)
    params = {'eps':1e-4, 'normalize':True, 'scale':0.5, 'cg_rate':2}
    sol_i = scs.solve(data, K, gpu=True, **params)

if __name__ == "__main__":
   main()
