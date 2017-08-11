"""
THIS PROGRAM IS USED AS MODULE. ALL THE FUNCTIONS ARE RUN THROUGH OUTPUT FUNCTION.
"""

import numpy as np
import matplotlib.pyplot as plt


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

def splined(a, b, h):
    r = abs(a - b) / h
    rij = abs(a - b)
    if (0 < r < 1.0):
        w = (1.5 * (r**2) - 2 * r)
        w = (a - b) * w / rij
    elif(1.0 <= r < 2.0):
        w = (-0.5 * ((2 - r)**2))
        w = (a - b) * w / rij
    else:
        w = 0
    return (w / h**2)


def gaussd(a, b, h):
    r = abs(a - b) / h
    rij = abs(a - b)
    if(0 < r < 3.0):
        w = -2 * r * (np.exp(-(r**2)))
        w = (a - b) * w / rij
    else:
        w = 0
    return (w / (h * h * (np.sqrt(np.pi))))


def sph_approx(x, f, h, kernel=gauss):
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



def sph_approx_density(x, f, h, m, kernel=gauss):
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


def sph_approx_density_derivative(x, f, h, m, kernel=gauss, kerneld=gaussd):
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
        print rho
        for j in range(0, n):
            w = kerneld(x[i], x[j], h)
            t = w * f[j] * m / rho
            v = v + t
        y.append(v)
    return y


def generatex(lower, upper, m, rho):
    dx = m / rho
    n = 1 + (upper - lower) / dx
    x = np.mgrid[lower:upper:n * 1j]
    return x



def solve_advection(x, y, total_time, a, dt):
    upper = max(x)
    lower = min(x)
    n = len(x)
    dx = (upper - lower) / (len(x) - 1)
    q = [min(x) - dx]
    x = np.concatenate([q, x])
    # print x
    r = [y[n - 2]]
    y = np.concatenate([r, y])
    iter = int(total_time / dt)
    for t in range(iter):
        for i in range(len(x)):
            x[i] = x[i] + dt * a
            # print 1
            if(x[i] > upper):
                c = [x[0] - dx]
                x = np.concatenate([c, x])
                x = np.delete(x, i + 1)
                b = [y[i - 2]]
                y = np.concatenate([b, y])
                y = np.delete(y, i + 1)
    # print x
    # print y
    return x, y


def berger_xsph(x, y, total_time, dt, e, h):
    upper = max(x)
    lower = min(x)
    n = len(x)
    dx = (upper - lower) / (n - 1)
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
        if(iter < 100):
            d = 10
        else:
            d = 100
        if(k % d == 0):
            print k
            #plt.plot(x, y, label=str(k*dt))
            #plt.ylim([0, 1])
    #plt.legend()
    plt.show()
    return x, y



if __name__ == '__main__':
    n = 101
    x = np.mgrid[0:1:n*1j]
    #u = [0] * n
    #b = 1.0/3.0
    dx = 1.0/(n-1)
    #rho = 1.0
    f = np.sin(2*np.pi*x)
    dx = (max(x) - min(x))/(n-1)
    #m = rho * dx
    total_time = 0.25
    dt = 0.1
    xx, yy = berger_xsph(x, f, total_time, 0.1*dx, 0.5, dx)
    plt.plot(xx, yy)
    plt.show()
