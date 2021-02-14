import sys
import numpy as np
from collections import defaultdict
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

map_pro_dat = defaultdict(list)

#load scale
scale = []
for s in sys.argv[2:]:
    scale.append(float(s))
#sys.exit()

#read data for each protein group
lines = open(sys.argv[1], 'r').readlines()
for l in lines[1:]:
    es = l.strip().split('\t')
    for i in [1,2,3,4,5,6]:
        val = float(es[i])*scale[i-1]
        map_pro_dat[es[0]].append(val)
#print map_pro_dat

def func(x, h, c, A):
    return A/(1.0+(x/c)**h)

x0 = [37.0, 42.0, 47.0, 52.0, 57.0, 62.0]
p0 = [10.0, 57.0, 100.0] #initial guess
b0 = [(0, 30.0, 95.0), (40.0, 90.0, 105.0)]
#p0 = [10.0, 57.0]
for p in map_pro_dat.keys(): #for p in ["sp|Q9C0B1|FTO_HUMAN"]:
  try:
    xdata = x0
    ydata = map_pro_dat[p]
    ydata = [y*100.0/np.max(ydata[:3]) for y in ydata]
    y0 = ydata
    popt, pcov = curve_fit(func, xdata, ydata, p0, bounds=b0)
    perr = np.sqrt(np.diag(pcov))
    Tm = popt[1]
    Hill = popt[0]
    Tm_err = perr[1]
    #print x0
    #print y0

    #fix
    if Tm_err>5.0:
        #if Tm<52:
        xdata = [27.0] + xdata
        ydata = [100.0] + ydata
        #if Tm>47 or np.min(ydata)>10.0:
        xdata = xdata + [92.0]
        ydata = ydata + [0.0]
        popt, pcov = curve_fit(func, xdata, ydata, p0, bounds=b0)
        perr = np.sqrt(np.diag(pcov))
        Tm = popt[1]
        Hill = popt[0]
        Tm_err = perr[1]

    #output
    y_calc0 = [func(i, popt[0], popt[1], popt[2]) for i in x0]
    #print x0, y_calc0, y0
    R2 = np.corrcoef(y0, y_calc0)
    y_calc = [func(i, popt[0], popt[1], popt[2]) for i in xdata]
    #print xdata, y_calc, ydata
    #pR2 = pearsonr(ydata, y_calc)
    #pR2 = pearsonr(y0, y_calc) #same as R2
    for pname in p.split(';'):
      print pname, "Tm=", Tm, "+/-", Tm_err, "Hill=", Hill, "R2=", R2[0,1] #,pR2[0]
    #debug
    continue
    #plt.scatter(ydata, y_calc, c='b')
    #plt.scatter(y0, y_calc0, c='r')
    plt.scatter(xdata, ydata, c='b')
    x2 = np.linspace(30.0, 70.0, 201)
    y2 = [func(i, popt[0], popt[1], popt[2]) for i in x2]
    plt.plot(x2, y2, 'r')
    plt.show()
  except:
    ydata = map_pro_dat[p]
    print "debug:", p, ydata

