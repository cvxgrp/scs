#!/usr/bin/env python -u
# encoding: utf-8

import numpy as np
from scipy import interpolate, ndimage
import gc
import os
import psutil
from cvxpy import *
import numpy

import objgraph

def sample_mem():
    gc.collect()
    process = psutil.Process(os.getpid())
    return process.get_memory_info()[0] / float(2 ** 20)

def cvx_solve():
    m = 128
    n = 10
    np.random.seed(1)
    A = np.random.randn(m, n)
    b = np.random.randn(m)

    # Construct the problem.
    x = Variable(n)
    objective = Minimize(norm(A*x - b))
    constraints = [0 <= x, x <= 1]
    prob = Problem(objective, constraints)
    prob.solve('SCS')
    # prob.solve('ECOS')
    # prob.solve('CVXOPT')

if __name__=='__main__':
    for _ in range(10000):
        cvx_solve()
        print(sample_mem())

