import numpy as np
import matplotlib.pyplot as plt
import SPHA


def function(x):
    return np.sin(2*np.pi*x)

def sph_approx(x, f, h, kernel=SPHA.gauss):
    n = len(x)
    dx = (max(x) - min(x)) / (n - 1)
    y = []
    for i in range(0, n):
        v = 0.0
        for j in range(0, n):
            w = kernel(x[i], x[j], h)
            t = w * f[j] * dx
            v = v + t
        y.append(v)
    return y

def sph_approx_density(x, f, h, m, kernel=SPHA.gauss):
    n = len(x)
    dx = (max(x) - min(x)) / (n - 1)
    # print dx
    y = []
    for i in range(0, n):
        v = 0.0
        p = 0.0
        for j in range(0, n):
            w = kernel(x[i], x[j], h)
            p = p + (m * w)
        rho = p
        # print rho
        for j in range(0, n):
            w = kernel(x[i], x[j], h)
            t = w * f[j] * m / rho
            v = v + t
        y.append(v)
    return y

def sph_approx_density_interior(x, x_int, f, h, m, kernel=SPHA.gauss):
    n = len(x)
    n1 = len(x_int)
    dx = (max(x) - min(x)) / (n - 1)
    y = []
    for i in range(0, n1):
        v = 0.0
        p = 0.0
        for j in range(0, n):
            w = kernel(x_int[i], x[j], h)
            p = p + (m * w)
        rho = p
        # print rho
        for j in range(0, n):
            w = kernel(x_int[i], x[j], h)
            t = w * f[j] * m / rho
            v = v + t
        y.append(v)
    return y




def compute_approx(n, kernel_factor, kernel=SPHA.spline):
    x = np.linspace(-1, 1, n)
    dx = x[1] - x[0]
    x = modify_x(x, 0.5)  ### x = modify_x(x, 0.0) for no noise addition
    x_int = np.linspace(-1, 1, 5*(n-1)+1)
    func = function(x)
    func1 = function(x_int)
    h = kernel_factor * dx
    func_approx = sph_approx_density_interior(x, x_int, func, h, 1.0, kernel)
    return x_int, func_approx, func1

def compute_approx_new(n, kernel_factor, kernel=SPHA.spline):
    x = np.linspace(-1, 1, n)
    x = modify_x(x, 1.0)
    x_int = np.linspace(-1, 1, n)
    # x_int = modify_x(x_int, 0.5)
    func = function(x)
    func1 = function(x_int)
    dx = x[1] - x[0]
    h = kernel_factor * dx
    func_approx = sph_approx_interior(x, x_int, func, h, kernel)
    return func_approx, func1


def modify_x(x, c):
    n = len(x)
    dx = (max(x) - min(x)) / (n - 1)
    c = dx * c
    for i in range(1, n - 1):
        x[i] = x[i] + np.random.random_sample() * c
    return x


def compute_error(N, kernel_factor, kernel=SPHA.spline):
    l1_err = []
    linf_err = []
    l2_err = []
    for n in N:
        print n
        x, func_approx, func_exact = compute_approx(n, kernel_factor, kernel)
        # func_approx, func_exact = compute_approx_new(n, kernel_factor, kernel)
        func_approx = np.asarray(func_approx)
        err = abs(func_approx - func_exact)
        linf_err.append(max(err))
        l1_err.append(sum(err)/len(func_approx))
        l2_err.append(np.sqrt(sum(err**2))/len(func_approx))
    return l1_err, linf_err, l2_err


def plots(N, l1_err, linf_err, l2_err):
    plt.plot(N, l1_err, label='L1')
    plt.plot(N, linf_err, label='Linf')
    plt.plot(N, l2_err, label='L2')
    plt.legend()
    plt.show()

def log_plots(N, l1_err, linf_err, l2_err):
    plt.loglog(N, l1_err, label='L1')
    plt.loglog(N, linf_err, label='Linf')
    plt.loglog(N, l2_err, label='L2')
    N = np.asarray(N)
    plt.plot(N, 2.0*(pow(N, -3.0)), '--', label='$N^{2.0}$')
    plt.plot(N, 1.0*(pow(N, -1.0)), '--', label='$N^{1.0}$')
    plt.legend()     
    plt.show()


def quintic_spline(a, b, h):
    sigma = 1.0/120
    r = abs(a - b) / h
    if(r <= 1.0):
        w = sigma * ((3 - r)**5 - 6*(2 - r)**5 + 15*(1 - r)**5)
    elif(1.0 < r <= 2.0):
        w = sigma * ((3 - r)**5 - 6*(2 - r)**5)
    elif(2.0 < r <= 3.0):
        w = sigma * ((3 - r)**5)
    else:
        w = 0
    return (w / h)


def super_gauss(a, b, h):
    r = abs(a - b) / h
    if(r <= 3.0):
        w = (np.exp(-r**2)*(1.5 - r**2))
    else:
        w = 0
    return w / (h * (np.sqrt(np.pi)))

def wendland_c2(a, b, h):
    sigma = 5.0/8.0
    r = abs(a - b) / h
    if(r <= 2.0):
        w = sigma * ((1.0 - 0.5*r)**3) * (3.0*0.5*r + 1)
    else:
        w = 0
    return (w / h)


def wendland_c4(a, b, h):
    sigma = 0.75
    r = abs(a - b) / h
    if(r <= 2.0):
        w = sigma * ((1 - 0.5*r)**5) * (5*0.5*r + 1 + 2.0*(r**2))
    else:
        w = 0
    return (w / h)


def wendland_c6(a, b, h):
    sigma = 55.0/64.0
    r = abs(a - b) / h
    if(r <= 1.0):
        w = sigma * ((1 - 0.5*r)**7) * (7*0.5*r + 1 + 19.0*(r**2)*0.25 + 21*(r**3)*0.125)
    else:
        w = 0
    return (w / h)

def spline(a, b, h):
    r = abs(a - b) / h
    if(r < 1.0):
        w = ((2 / 3.0) - r**2 + r**3 / 2)
    elif(1.0 <= r < 2.0):
        w = ((1.0 / 6) * (2 - r)**3)
    else:
        w = 0
    return (w / h)


def gauss(a, b, h):
    r = abs(a - b) / h
    if(r < 3.0):
        w = (np.exp(-r**2))
    else:
        w = 0
    return w / (h * (np.sqrt(np.pi)))

def sph_approx(x, f, h, kernel=super_gauss):
    n = len(x)
    dx = (max(x) - min(x)) / (n - 1)
    y = []
    for i in range(0, n):
        v = 0.0
        for j in range(0, n):
            w = kernel(x[i], x[j], h)
            t = w * f[j] * dx
            v = v + t
        y.append(v)
    return y

def sph_approx_interior(x, x_int, f, h, kernel=super_gauss):
    n = len(x)
    n1 = len(x_int)
    dx = x[1] - x[0]
    y = []
    for i in range(n1):
        v = 0.0
        for j in range(0, n):
            w = kernel(x_int[i], x[j], h)
            t = w * f[j] * dx
            v = v + t
        y.append(v)
    return y
        
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


def plotting(N, kernel_factor=1.0, kernel=spline):
    l1, linf, l2 = compute_error(N, kernel_factor, kernel)
    c1, c2, cinf = get_coefficients(N, l1, l2, linf)
    plt.figure(figsize=(9,7))
    plt.loglog(N, l1, label='L1')
    plt.loglog(N, linf, label='Linf')
    plt.loglog(N, l2, label='L2')
    N = np.asarray(N)
    plt.plot(N, 10**(c1[1])*pow(N, c1[0]), 'd', label='$N^{' + str(c1[0]) + '}$')
    plt.plot(N, 10**(c2[1])*pow(N, c2[0]), 'v', label='$N^{' + str(c2[0]) + '}$')
    plt.plot(N, 10**(cinf[1])*pow(N, cinf[0]), '*', label='$N^{' + str(cinf[0]) + '}$')
    plt.ylabel('Error Norm')
    plt.xlabel('Number of particles')
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


N = [51, 101, 201, 501, 1001, 2001]
plotting(N, kernel=quintic_spline)