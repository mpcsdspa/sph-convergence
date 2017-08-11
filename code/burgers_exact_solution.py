import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

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

### Plotting Functions ###


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
    plt.loglog(N, l1, label='L1')
    plt.loglog(N, linf, label='Linf')
    plt.loglog(N, l2, label='L2')
    N = np.asarray(N)
    plt.plot(N, 10**(c1[1])*pow(N, c1[0]), 'd', label='$N^{' + str(c1[0]) + '}$')
    plt.plot(N, 10**(c2[1])*pow(N, c2[0]), 'v', label='$N^{' + str(c2[0]) + '}$')
    plt.plot(N, 10**(cinf[1])*pow(N, cinf[0]), '*', label='$N^{' + str(cinf[0]) + '}$')
    plt.ylabel('N')
    plt.xlabel('Error Norm')
    plt.legend()
    plt.show()

def compute_coefficients(N, a1, a2, a3, a4):
    logA = np.log10(N)
    logB = np.log10(a1)
    logC = np.log10(a2)
    logD = np.log10(a3)
    logE = np.log10(a4)
    c1 = np.polyfit(logA, logB, 1)
    c2 = np.polyfit(logA, logC, 1)
    c3 = np.polyfit(logA, logD, 1)
    c4 = np.polyfit(logA, logE, 1)
    c1 = [round(i,2) for i in c1]
    c2 = [round(i,2) for i in c2]
    c3 = [round(i,2) for i in c3]
    c4 = [round(i,2) for i in c4]
    a = (c1[0] + c2[0] + c3[0] + c4[0])*0.25
    b = (c1[1] + c2[1] + c3[1] + c4[1])*0.25
    return round(a,2), round(b,2)

def misc_plotting(N, kernel_factor=1.0, norm='L2'):
    l1_spline, linf_spline, l2_spline = compute_error(N, kernel_factor, kernel=spline)
    l1_gauss, linf_gauss, l2_gauss = compute_error(N, kernel_factor, kernel=gauss)
    l1_quintic, linf_quintic, l2_quintic = compute_error(N, kernel_factor, kernel=quintic_spline)
    l1_wc2, linf_wc2, l2_wc2 = compute_error(N, kernel_factor, kernel=wendland_c2)
    l1_wc4, linf_wc4, l2_wc4 = compute_error(N, kernel_factor, kernel=wendland_c4)
    l1_supergauss, linf_supergauss, l2_supergauss = compute_error(N, kernel_factor, kernel=super_gauss)
    N = np.asarray(N)
    if norm=='L2':
        a,b = compute_coefficients(N, l2_spline, l2_gauss, l2_supergauss, l2_quintic)
        print a, b
        plt.figure(figsize=(9,7))
        plt.loglog(N, l2_gauss, label='Gauss')
        plt.loglog(N, l2_spline, label='Cubic-Spline')
        plt.loglog(N, l2_quintic, label='Quintic-Spline')
        plt.loglog(N, l2_supergauss, label='Super-Gaussian')
        plt.plot(N, 10**(b)*pow(N, a), 'v', label='$N^{' + str(a) + '}$')
        plt.legend()
        plt.show()


n = 101
x, f = exact_solution(n, 0.1)
u1 = exact_newton(x, 0.1)
print np.asarray(u1) - np.asarray(f)
plt.plot(x, u1, label='newton')
plt.plot(x, f, label='integration')
plt.show()