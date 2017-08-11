import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def compute_error(u_exact, u_sph):
    n = len(u_exact)
    u_exact = np.asarray(u_exact)
    u_sph = np.asarray(u_sph)
    l2 = np.sqrt((sum((u_exact - u_sph)**2)) / len(u_exact))
    l1 = sum(abs(u_exact - u_sph)) / n
    linf = max(abs(u_exact - u_sph))
    return l1, l2, linf

def exact_solution(n, total_time):
    x = np.linspace(0, 1, n)
    f = np.sin(2*np.pi*x)
    dt = 0.0001
    p = int(total_time/dt)
    for t in range(p):
        for i in range(n):
            x[i] = x[i] + f[i]*dt
    return x, f

def exact_newton(x, tf):
    n = len(x)
    u1 = np.zeros(n)
    u1 = u1 + 10.0
    for i in range(n):
        f = lambda u : u - np.sin(2*np.pi*(x[i] - u*tf))
        initial = 0.0
        while (abs(check_root(u1[i], x[i], tf)) > 1e-15):
            u1[i] = fsolve(f, initial)
            initial = initial + 0.1
    return u1

def check_root(u, x, tf):
    a = u - np.sin(2*np.pi*(x - u*tf))
    return a

tf = 0.15
N = [51, 101, 201, 501, 1001, 2001]
L1 = []
L2 = []
Linf = []
for i in range(len(N)):
    with open('xx' + str(N[i]) + '.txt') as f:
        xx = f.read().splitlines()
    with open('yy' + str(N[i]) + '.txt') as f:
        yy = f.read().splitlines()
    xx = [float(j) for j in xx]
    yy = [float(j) for j in yy]
    x = np.linspace(0, 1, N[i])
    # x_exact, u_exact = exact_solution(N[i], tf)
    u_exact = exact_newton(xx, tf)
    l1, l2, linf = compute_error(u_exact, yy)
    # l1, l2, linf = compute_error(x_exact, xx)
    L1.append(l1)
    L2.append(l2)
    Linf.append(linf)



def get_coefficients(N, L1, L2, Linf):
    logA = np.log10(N)
    logB = np.log10(L1)
    logC = np.log10(L2)
    logD = np.log10(Linf)
    c1 = np.polyfit(logA, logB, 1)
    c2 = np.polyfit(logA, logC, 1)
    cinf = np.polyfit(logA, logD, 1)
    c1 = [round(i,2) for i in c1]
    c2 = [round(i,2) for i in c2]
    cinf = [round(i,2) for i in cinf]
    return c1, c2, cinf


def plotting(N, l1, l2, linf):
    c1, c2, cinf = get_coefficients(N, l1, l2, linf)
    plt.figure(figsize=(9,7))
    N = np.asarray(N)
    plt.plot(N, 10**(c1[1])*pow(N, c1[0]), label='$N^{' + str(c1[0]) + '}$')
    plt.plot(N, 10**(c2[1])*pow(N, c2[0]), label='$N^{' + str(c2[0]) + '}$')
    plt.plot(N, 10**(cinf[1])*pow(N, cinf[0]), label='$N^{' + str(cinf[0]) + '}$')
    plt.loglog(N, l1, 'd', label='L1')
    plt.loglog(N, linf, 'v', label='Linf')
    plt.loglog(N, l2, '*', label='L2')
    plt.xlabel('N')
    plt.ylabel('Error Norm')
    plt.legend()
    plt.show()

plotting(N, L1, L2, Linf)
