import sys
import numpy as np
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, h, c, A):
    return A/(1.0+(x/c)**h)

T = [37.0, 42.0, 47.0, 52.0, 57.0, 62.0] #, 72.0]
S = ['blue', 'red', 'blue', 'red', 'blue', 'red', 'blue', 'red', 'blue', 'red']
p0 = [10.0, 57.0, 100.0] #initial guess
b0 = [(0, 30.0, 95.0), (40.0, 90.0, 105.0)]
T0 = 27.0
TN = 92.0

data = []
for fn in sys.argv[1:]:
    d1 = []
    d2 = []
    d3 = []
    d4 = []
    d5 = []
    d6 = []
    lines = open(fn,'r').readlines()
    for l in lines[1:]:
        es = l.strip().split('\t')
        d1.append(float(es[1]))
        d2.append(float(es[2]))
        d3.append(float(es[3]))
        d4.append(float(es[4]))
        d5.append(float(es[5]))
        d6.append(float(es[6]))
    data.append([])
    data[-1].append(np.mean(d1))
    data[-1].append(np.mean(d2))
    data[-1].append(np.mean(d3))
    data[-1].append(np.mean(d4))
    data[-1].append(np.mean(d5))
    data[-1].append(np.mean(d6))

#for each TMT
topI = np.median([x[0] for x in data] + [x[1] for x in data])
x1 = []
y1 = []
for n, dat in enumerate(data):
    x1 = x1 + [T0] + T + [TN]
    y1 = y1 + [topI] + dat + [0.0]

x1 = np.array(x1)
y1 = np.array(y1)/topI*100
scale = np.ones(y1.shape)
nd = int(y1.shape[0] / len(data))
skip = []
for i in range(len(data)):
    skip.append(i*nd)
    skip.append((i+1)*nd-1)

tol = 1.0
max_step = 100
for i in range(max_step):
    #fitting
    ys = y1*scale #scaled
    popt, pcov = curve_fit(func, x1, ys, p0, bounds=b0)
    #fit curve
    x2 = np.linspace(30.0, 80.0, 201)
    y2 = [func(i, popt[0], popt[1], popt[2]) for i in x2]
    #plt.plot(x2, y2, 'black')
    #split dataset
    for n, dat in enumerate(data):
        xdata = np.array([T0] + T + [TN])
        ydata = np.array([topI] + dat + [0.0]) / topI * 100
        s = scale[nd*n:nd*(n+1)]
        #raw data
        #plt.scatter(xdata, ydata, s = 100, color="", edgecolors=S[n], marker='o', lw=2)
        #scaled data
        #plt.scatter(xdata, ydata*s, s = 20, color=S[n], marker='o')
    #plt.show()

    #check diff, change scale
    yfit = [func(i, popt[0], popt[1], popt[2]) for i in x1] #target
    delta = yfit - ys
    tmp = []
    for j, d2 in enumerate(delta * delta):
        tmp.append([d2, j])
    for d2, j in sorted(tmp, key = lambda x: x[0], reverse=True):
        if j not in skip:
            k = j
            break
    print(k, ys[k], yfit[k])
    if (ys[k]-yfit[k])**2 < tol:
        print("Done!")
        break
    new_scale = yfit[k] / y1[k]
    scale[k] = new_scale

for n, fn in enumerate(sys.argv[1:]):
    s = scale[nd*n:nd*(n+1)]
    print(fn, s[1:-1])
