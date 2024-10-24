import sys
import re
import numpy as np
from collections import defaultdict
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

def func(x, h, c, A):
    return A/(1.0+(x/c)**h)

def get_Ts(temp_str):
    Ts = temp_str.strip('"')
    Ts = re.split(r'[,\s\t]+', Ts)
    Ts = [float(t) for t in Ts]
    return Ts

T = get_Ts(sys.argv[1])
print("Temperature", T, file=sys.stderr)

color_lst = ['r', 'b', 'g', 'y', 'p'] #max file number = 5
p0 = [10.0, 57.0, 100.0] #initial guess
b0 = [(0, 30.0, 95.0), (40.0, 90.0, 105.0)]
T0 = 27.0
TN = 97.0

#load file name and scale
databases = []
prot_count = defaultdict(int)
nf = 0
for l in open(sys.argv[2], 'r').readlines():
    #for each data file, typically two: ctrl, exp
    map_pro_dat = defaultdict(list)
    es = l.strip().split(" ")
    fn = es[0]
    scale = []
    for s in es[1:]:
        scale.append(float(s))
    print("Loading:", fn, scale, file=sys.stderr)
    nf += 1

    lines = open(fn, 'r').readlines()
    for l in lines[1:]:
        es = re.split(r'[,\s\t]+', l.strip())
        prot = es[0]
        for i, s in enumerate(scale):
            val = float(es[i+1])*s
            map_pro_dat[prot].append(val)
        assert(len(T) == len(map_pro_dat[prot]))
        prot_count[prot] += 1

    databases.append(map_pro_dat)

tol_sq = 2.0
tol_R2 = 0.95
with PdfPages('combined-results.pdf') as pdf:
    for p in tqdm(prot_count.keys(), desc="Processing proteins"):
        if prot_count[p] < nf: continue
        with_in_tol = True
        skip_output = False
        R2_lst = []
        ERR_lst = []
        p0_lst = []
        norm_lst = []
        for i in range(nf):
            #for each file
            data = databases[i][p]
            topI = np.max(data)
            norm_data = [ d/topI*100 for d in data ]
            norm_lst.append(norm_data)
            xdata = [T0] + T + [TN]
            ydata = [100.0] + norm_data + [0.0]
            popt, pcov = curve_fit(func, xdata, ydata, p0, bounds=b0)
            perr = np.sqrt(np.diag(pcov))
            Tm = popt[1]
            Hill = popt[0]
            p0_lst.append(popt)
            Tm_err = perr[1]

            ycalc0 = [func(i, popt[0], popt[1], popt[2]) for i in xdata]
            R2 = np.corrcoef(ydata, ycalc0)
            R2 = R2[0,1]
            R2_lst.append(R2)
            ERR_lst.append(Tm_err)
            if R2 < tol_R2 or Tm_err > tol_sq:
                with_in_tol = False
            if R2 < 0.8 or Tm_err > 10.0:
                skip_output = True
                break
        
        if skip_output: continue
        #output
        out_str = "PRO= " + p
        legend_lst = []
        for i in range(nf):
            legend_lst.append(f"Tm{i+1}= {p0_lst[i][1]:.2f} ({R2_lst[i]:.2f}, {ERR_lst[i]:.2f})")
            out_str += f" Tm{i+1}= {p0_lst[i][1]:.2f} R2= {R2_lst[i]:.2f} ERR={ERR_lst[i]:.2f}"
        print(out_str)

        if with_in_tol:
            #print(p, p0_lst, R2_lst)
            #draw
            plt.figure()
            for i in range(nf):
                x1 = T
                y1 = norm_lst[i]
                x2 = np.linspace(30.0, 80.0, 201)
                y2 = [func(x, p0_lst[i][0], p0_lst[i][1], p0_lst[i][2]) for x in x2]
                #fitting curve
                plt.plot(x2, y2, color=color_lst[i], label=legend_lst[i])
                #norm data
                plt.scatter(x1, y1, s = 100, color="w", edgecolors=color_lst[i], marker='o', lw=2)
                plt.title(p)
                plt.legend(loc='upper right')
            pdf.savefig()
            plt.close()
            #output stderr
