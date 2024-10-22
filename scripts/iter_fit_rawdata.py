import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

def func(x, h, c, A):
    return A/(1.0+(x/c)**h)

def get_Ts(temp_str):
    Ts = temp_str.strip('"')
    Ts = re.split(r'[,\s\t]+', Ts)
    Ts = [float(t) for t in Ts]
    return Ts

#"37.0, 42.0, 47.0, 52.0, 57.0, 62.0"
T = get_Ts(sys.argv[1])
nT = len(T)

p0 = [10.0, 57.0, 100.0] #initial guess
b0 = [(0, 30.0, 95.0), (40.0, 90.0, 105.0)]
T0 = 27.0
TN = 97.0
bg = 2.0

data = []
for fn in sys.argv[2:]:
    #for each file
    ds = [[] for _ in range(nT)]
    lines = open(fn,'r').readlines()
    for l in lines[1:]:
        es = re.split(r'[,\s\t]+', l.strip())
        for i, e in enumerate(es[1:]):
            ds[i].append(float(e))
        
    fn_data = [np.mean(dat_lst) for dat_lst in ds]
    data.append(fn_data)

#for each TMT
norm_data = []
for dat in data:
    topI = np.max(dat[:2])
    dat = [ d/topI*100 for d in dat ]
    norm_data.append(dat)

x1 = []
y1 = []
for n, dat in enumerate(norm_data):
    x1 = x1 + [T0] + T + [TN]
    y1 = y1 + [100.0] + dat + [bg]

x1 = np.array(x1)
y1 = np.array(y1)
scale = np.ones(y1.shape)
nd = int(y1.shape[0] / len(data))
skip = []
for i in range(len(data)):
    skip.append(i*nd)
    skip.append((i+1)*nd-1)

print("DB:", [f"{x:.2f}" for x in x1] )
print("DB:", [f"{x:.2f}" for x in y1] )

tol = 5.0
max_step = 100
with PdfPages('fitting.pdf') as pdf:
    for i in range(max_step):
        plt.figure()
        #fitting
        ys = y1*scale #scaled
        popt, pcov = curve_fit(func, x1, ys, p0, bounds=b0)
        #fit curve
        x2 = np.linspace(30.0, 80.0, 201)
        y2 = [func(i, popt[0], popt[1], popt[2]) for i in x2]
        plt.plot(x2, y2, 'black')
        #split dataset
        S="b"
        for n, dat in enumerate(norm_data):
            if S=="b": S = "r"
            elif S=="r": S = "b"
            xdata = np.array([T0] + T + [TN])
            ydata = np.array([100.0] + dat + [bg])
            s = scale[nd*n:nd*(n+1)]
            #raw data
            plt.scatter(xdata, ydata, s = 100, color="w", edgecolors=S, marker='o', lw=2)
            #scaled data
            plt.scatter(xdata, ydata*s, s = 20, color=S, marker='o')
        #plt.show()
        pdf.savefig()
    
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
        print(k, f"{ys[k]:.2f}", f"{yfit[k]:.2f}")
        if (ys[k]-yfit[k])**2 < tol:
            print("Done!")
            break
        new_scale = yfit[k] / y1[k]
        scale[k] = new_scale

for n, fn in enumerate(sys.argv[2:]):
    s = scale[nd*n:nd*(n+1)]
    print(fn, " ".join([f"{x:.2f}" for x in s[1:-1]]))
