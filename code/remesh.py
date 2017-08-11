import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy import interpolate

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


def check_root(u, x, tf):
    a = u - np.sin(2*np.pi*(x - u*tf))
    return a


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


def berger_remesh(x, y, total_time, dt, e, h):
    dx = x[1] - x[0]
    n = len(x)
    iter = int(total_time / dt)
    for k in range(iter):
        dx = x[1] - x[0]
        h = dx
        for i in range(n):
            v = 0.0
            for j in range(0, n):
                w = wendland_c4(x[i], x[j], h)
                t = w * (y[j] - y[i]) * dx
                v = v + t
            a = y[i] + e * v
            x[i] = x[i] + a * dt
        if (k % 10 == 0):
            print k
            x, y = remesh(x, y, n)
    return x, y


def remesh(x, y, n):
    fint = interpolate.interp1d(x, y, kind='cubic')
    xnew = np.linspace(min(x), max(x), n)
    ynew = fint(xnew)
    return xnew, ynew


def compute_error(y, y_exact):
    y = np.asarray(y)
    y_exact = np.asarray(y_exact)
    err = abs(y_exact - y)
    l1_err = sum(err) / len(y)
    linf_err = max(err)
    l2_err = np.sqrt(sum(err**2)) / len(y)
    return l1_err, linf_err, l2_err



tf = 0.15
N = [51, 101, 201, 251, 501, 1001]
L1 = []
L2 = []
Linf = []
for i in range(len(N)):
    print N[i]
    x = np.linspace(0, 1, N[i])
    dx = x[1] - x[0]
    f = np.sin(2 * np.pi * x)
    xx, yy = berger_remesh(x, f, tf, 0.1*dx, 0.5, dx)

    thefile = open('xx' + str(N[i]) + '.txt', 'w')
    for item in xx:
        thefile.write("%s\n" % item)

    thefile = open('yy' + str(N[i]) + '.txt', 'w')
    for item in yy:
        thefile.write("%s\n" % item)