import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import fsolve
import csv

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

def berger_xsph1(x, y, total_time, dt, e, h):
    dx = x[1] - x[0]
    n = len(x)
    iter = int(total_time / dt)
    for k in range(iter):
        for i in range(n):
            v = 0.0
            for j in range(0, n):
                if i!=n-1:
                    h = x[i+1] - x[i]
                else:
                    h = x[i] - x[i-1]
                w = wendland_c6(x[i], x[j], h)
                t = w * (y[j] - y[i]) * dx
                v = v + t
            a = y[i] + e * v
            x[i] = x[i] + a * dt
    return x, y

def berger_xsph(x, y, total_time, dt, e, h):
    dx = x[1] - x[0]
    n = len(x)
    iter = int(total_time / dt)
    for k in range(iter):
        for i in range(n):
            v = 0.0
            for j in range(0, n):
                w = spline(x[i], x[j], h)
                t = w * (y[j] - y[i]) * dx
                v = v + t
            a = y[i] + e * v
            x[i] = x[i] + a * dt
    return x, y

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

def compute_error(y, y_exact):
    y = np.asarray(y)
    y_exact = np.asarray(y_exact)
    err = abs(y_exact - y)
    l1_err = sum(err) / len(y)
    linf_err = max(err)
    l2_err = np.sqrt(sum(err**2)) / len(y)
    return l1_err, linf_err, l2_err


def berger_xsph_ADKE(x, y, total_time, dt, e, h):
    n = len(x)
    dx = x[1] - x[0]
    h_list = [dx]*n
    iter = int(total_time / dt)
    for k in range(iter):
        for i in range(n):
            v = 0.0
            h_list = get_lambda(x, 1.0, h_list, 0.03, 0.5)*dx
            for j in range(0, n):
                w = gauss(x[i], x[j], h_list[i])
                t = w * (y[j] - y[i]) * dx
                v = v + t
            a = y[i] + e * v
            x[i] = x[i] + a * dt
    return x, y


def get_lambda(x, m, h, k, eps):
    m = 1.0
    rho = get_density(x, m, h)
    g = np.sum(np.log(rho))/len(rho)
    g = 10**(g)
    l = np.zeros(len(rho))
    for i in range(len(rho)):
        l[i] = k * ((g / rho[i])**eps)
    return l

def get_density(x, m, h):
    rho = []
    for i in range(0, len(x)):
        p = 0.0
        for j in range(0, len(x)):
            w = gauss(x[i], x[j], h[i])
            p = p + (m * w)
        rho.append(p)
    return rho


tf = 0.15
N = [51, 101, 151, 201, 251]
L1 = []
L2 = []
Linf = []
for i in range(len(N)):
    print N[i]
    x = np.linspace(0, 1, N[i])
    dx = x[1] - x[0]
    f = np.sin(2 * np.pi * x)
    xx, yy = berger_xsph_ADKE(x, f, tf, 0.1*dx, 0.5, dx)

    thefile = open('xx' + str(N[i]) + '.txt', 'w')
    for item in xx:
        thefile.write("%s\n" % item)

    thefile = open('yy' + str(N[i]) + '.txt', 'w')
    for item in yy:
        thefile.write("%s\n" % item)