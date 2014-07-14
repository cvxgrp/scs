import scs
from numpy import *
from scipy import sparse, randn

#############################################
#      Generate random cone problems        #
#############################################

def genFeasible(K, n, density):
    m = getConeDims(K)
    
    z = randn(m,)
    z = symmetrizeSDP(z, K)  # for SD cones
    y = proj_dual_cone(z, K)  # y = s - z;
    s = y - z  # s = proj_cone(z,K)
    
    A = sparse.rand(m, n, density, format='csc')
    A.data = randn(A.nnz)
    x = randn(n)
    c = -transpose(A).dot(y)
    b = A.dot(x) + s
    
    data = {'A': A, 'b': b, 'c': c}
    return data, dot(c, x)


def genInfeasible(K, n):
    m = getConeDims(K)
    
    z = randn(m,)
    z = symmetrizeSDP(z, K)  # for SD cones
    y = proj_dual_cone(z, K)  # y = s - z;
    A = randn(m, n)
    A = A - outer(y, transpose(A).dot(y)) / linalg.norm(y) ** 2  # dense...
    
    b = randn(m);
    b = -b / dot(b, y);
    
    data = {'A':sparse.csc_matrix(A), 'b':b, 'c':randn(n)}
    return data

def genUnbounded(K, n):
    m = getConeDims(K)
    
    z = randn(m);
    z = symmetrizeSDP(z, K);  # for SD cones
    s = proj_cone(z, K);
    A = randn(m, n);
    x = randn(n);
    A = A - outer(s + A.dot(x), x) / linalg.norm(x) ** 2;  # dense...
    c = randn(n);
    c = -c / dot(c, x);
    
    data = {'A':sparse.csc_matrix(A), 'b':randn(m), 'c':c}
    return data    

def pos(x):
    return (x + abs(x)) / 2

def getConeDims(K):
    l = K['f'] + K['l']
    for i in range(0, len(K['q'])):
        l = l + K['q'][i];
    
    for i in range(0, len(K['s'])):
        l = l + K['s'][i] ** 2;

    l = l + K['ep'] * 3;
    l = l + K['ed'] * 3;
    return l

def proj_dual_cone(z, c):
    return z + proj_cone(-z, c)

def proj_cone(z, c):
    z = copy(z)
    free_len = c['f']
    lp_len = c['l']
    q = c['q']
    s = c['s']
    # free/zero cone
    z[0:free_len] = 0;
    # lp cone
    z[free_len:lp_len + free_len] = pos(z[free_len:lp_len + free_len])
    # SOCs
    idx = lp_len + free_len;
    for i in range(0, len(q)):
        z[idx:idx + q[i]] = proj_soc(z[idx:idx + q[i]])
        idx = idx + q[i]
    # SDCs
    for i in range(0, len(s)):
        z[idx:idx + s[i] ** 2] = proj_sdp(z[idx:idx + s[i] ** 2], s[i])
        idx = idx + s[i] ** 2
    # Exp primal
    for i in range(0, c['ep']):
        z[idx:idx + 3] = project_exp_bisection(z[idx:idx + 3])
        idx = idx + 3
    # Exp dual
    for i in range(0, c['ed']):
        z[idx:idx + 3] = z[idx:idx + 3] + project_exp_bisection(-z[idx:idx + 3])
        idx = idx + 3
    return z

def proj_soc(tt):
    tt = copy(tt)
    if len(tt) == 0:
        return
    elif len(tt) == 1:
        return pos(tt)
    
    v1 = tt[0]
    v2 = tt[1:];
    if linalg.norm(v2) <= -v1:
        v2 = zeros(len(v2))
        v1 = 0
    elif linalg.norm(v2) > abs(v1):
        v2 = 0.5 * (1 + v1 / linalg.norm(v2)) * v2
        v1 = linalg.norm(v2)
    tt[0] = v1
    tt[1:] = v2
    return tt

def proj_sdp(z, n):
    z = copy(z)
    if n == 0:
        return
    elif n == 1:
        return pos(z)
    z = reshape(z, (n, n), order='F')
    zs = (z + transpose(z)) / 2
    
    w, v = linalg.eig(zs)  # cols of v are eignvectors
    w = pos(w)
    z = dot(v, dot(diag(w), transpose(v)))
    return reshape(z, (n * n,))

def symmetrizeSDP(z, K):
    l = K['f'] + K['l']
    for i in range(0, len(K['q'])):
        l = l + K['q'][i]
    for i in range(0, len(K['s'])):
        n = K['s'][i]
        V = reshape(z[l:l + n * n], (n, n), order='F')
        V = (V + transpose(V)) / 2
        z[l:l + n * n] = reshape(V, (n * n,))
        l = l + n * n
    return z

def project_exp_bisection(v):
    v = copy(v)
    r = v[0]; s = v[1]; t = v[2]
    # v in cl(Kexp)
    if (s * exp(r / s) <= t and s > 0) or (r <= 0 and s == 0 and t >= 0):
        return v
    # -v in Kexp^*
    if (-r < 0 and r * exp(s / r) <= -exp(1) * t) or (-r == 0 and -s >= 0 and -t >= 0):
        return zeros(3,)
    # special case with analytical solution
    if r < 0 and s < 0:
        v[1] = 0
        v[2] = max(v[2], 0)
        return v

    x = copy(v)
    ub, lb = getRhoUb(v)
    for iter in range(0, 100):
        rho = (ub + lb) / 2
        g, x = calcGrad(v, rho, x)
        if g > 0:
            lb = rho
        else:
            ub = rho
        if (ub - lb < 1e-6):
            break
    return x

def getRhoUb(v):
    lb = 0
    rho = 2 ** (-3)
    g, z = calcGrad(v, rho, v)
    while g > 0:
        lb = rho
        rho = rho * 2
        g, z = calcGrad(v, rho, z)
    ub = rho
    return ub, lb

def calcGrad(v, rho, warm_start):
    x = solve_with_rho(v, rho, warm_start[1])
    if x[1] == 0:
        g = x[0]
    else:
        g = (x[0] + x[1] * log(x[1] / x[2]))
    return g, x

def solve_with_rho(v, rho, w):
    x = zeros(3)
    x[2] = newton_exp_onz(rho, v[1], v[2], w)
    x[1] = (1 / rho) * (x[2] - v[2]) * x[2]
    x[0] = v[0] - rho
    return x

def newton_exp_onz(rho, y_hat, z_hat, w):
    t = max(max(w - z_hat, -z_hat), 1e-6)
    for iter in range(0, 100):
        f = (1 / rho ** 2) * t * (t + z_hat) - y_hat / rho + log(t / rho) + 1
        fp = (1 / rho ** 2) * (2 * t + z_hat) + 1 / t
    
        t = t - f / fp
        if t <= -z_hat:
            t = -z_hat
            break
        elif t <= 0:
            t = 0
            break
        elif abs(f) < 1e-6:
            break
    return t + z_hat

if __name__ == "__main__":
   main()
