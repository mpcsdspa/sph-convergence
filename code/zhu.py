import numpy as np
# import matplotlib.pyplot as plt

def spline(a, b, h):
    r = abs(a - b) / h
    if(r < 1.0):
        w = ((2 / 3.0) - r**2 + r**3 / 2)
    elif(1.0 <= r < 2.0):
        w = ((1.0 / 6) * (2 - r)**3)
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


def wendland_c4(a, b, h):
    sigma = 0.75
    r = abs(a - b) / h
    if(r <= 2.0):
        w = sigma * ((1 - 0.5*r)**5) * (5*0.5*r + 1 + 2.0*(r**2))
    else:
        w = 0
    return (w / h)

def berger_zhu(x, xint, y, yint, total_time, dt, e, h):
    n = len(x)
    nint = len(xint)
    iter = int(total_time / dt)
    for k in range(iter):
        dx = xint[1] - xint[0]
        h = x[1] - x[0]
        for i in range(n):
            v = 0.0
            for j in range(0, nint):
                w = super_gauss(x[i], xint[j], h)
                t = w * (yint[j] - y[i]) * dx
                v = v + t
            a = y[i] + e * v
            x[i] = x[i] + a * dt
    return x, y







tf = 0.15
N = [64, 100, 225, 400, 625]
L1 = []
L2 = []
Linf = []
for i in range(len(N)):
    print N[i]
    x = np.linspace(0, 1, N[i])
    nnb = np.sqrt(N[i])*0.5*(N[i]-1) + N[i]
    xint = np.linspace(0, 1, nnb)
    dx = x[1] - x[0]
    f = np.sin(2 * np.pi * x)
    fint = np.sin(2*np.pi*xint)
    xx, yy = berger_zhu(x, xint, f, fint,  tf, 0.1*dx, 0.1, dx)

    thefile = open('xx' + str(N[i]) + '.txt', 'w')
    for item in xx:
        thefile.write("%s\n" % item)

    thefile = open('yy' + str(N[i]) + '.txt', 'w')
    for item in yy:
        thefile.write("%s\n" % item)