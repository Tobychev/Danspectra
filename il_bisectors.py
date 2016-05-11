import matplotlib.pyplot as pl
import scipy.interpolate as si
import scipy.optimize as ro
import matplotlib.cm as cm
import numpy as np
regnames = ["5053","5215","5654","6405","6449"]
import pdb

def bisect(lmbd,mnspec,line):
    #find highest points to the left and right of linecenter
    left, = np.where(lmbd < line.cent) ; 
    left  = left[mnspec[ left ].argmax()]
    right,= np.where(lmbd > line.cent)
    right = right[mnspec[ right ].argmax()]
    # Limit bisector to the height of lowest edge
    ys = np.linspace(line.dept,min([mnspec[left],mnspec[right]]),129)
    mf    = si.UnivariateSpline(lmbd[::-1],mnspec[::-1],s=0)
    bisec = np.zeros(len(ys)) 
    for i,y in enumerate(ys[1:]):
#        print(i,mf(lmbd[left])-y,mf(lin.cent)-y,mf(lmbd[right])-y)
        try:
            lf = ro.brentq(lambda x : mf(x)-y,lmbd[left],lin.cent)
            rg = ro.brentq(lambda x : mf(x)-y,lmbd[right],lin.cent)
        except:
            lf = rg = lin.cent
        bisec[i+1] = (lf+rg)/2 - lin.cent
    return bisec,ys

col = cm.Oranges

if True:
    subs =  ['Myst    505.358','Myst    521.557','Myst    565.450','Myst    640.576','Myst    644.913',]
    j = 0
    for i,regname in enumerate(regnames):
        rg = __import__(regname)
        for lin in rg.qu1lines:
            print(lin.name)
            if lin.name in subs:
                xs,ys = bisect(rg.qu1.lmbd[lin.idx],rg.qu1m[lin.idx],lin)
                pl.plot(xs[:-1],ys[:-1],'-',label=lin.name)    
                j+=1
    pl.legend(loc="best")
    pl.ylabel("Relative intensity")
    pl.xlabel("Distance from line center [nm]")
    pl.xlim(-0.00158, 0.00288)
    pl.ylim( 0.81045, 0.98785)
    #pl.show(block=False)
    pl.savefig("../thesis/figures/MystBisectors.png")

if False:

    alla = ['Ni I    504.885','C I     505.214','C_2     505.262','Ti I    505.287','Fe I    505.464','Fe I    505.850','Fe I    520.988','Ti II   521.153','Co I    521.269','Cu I    521.820','Fe I    521.970','Fe I    565.147','Si I    565.492','V I     565.744','Sc II   565.790','Fe I    566.134','Fe I    640.032','Si I    640.729','Fe I    640.802','Si I    640.867','Ca I    644.982','Co I Si I 645.018','Si I    645.231','tel H_2O 645.414','Co I    645.500','tel H_2O 645.889']

    subs =  ['C_2     505.262','C I     505.214','Fe I    505.464','Co I    521.269','Si I    640.867','tel H_2O 645.889']
    j = 0
    for i,regname in enumerate(regnames):
        rg = __import__(regname)
        for lin in rg.qu1lines:
            print(lin.name)
            if lin.name in subs:
                xs,ys = bisect(rg.qu1.lmbd[lin.idx],rg.qu1m[lin.idx],lin)
                pl.plot(xs[:-1],ys[:-1],'-',label=lin.name)    
                j+=1
    pl.legend(loc="best")
    pl.ylabel("Relative intensity")
    pl.xlabel("Distance from line center [nm]")
#    pl.xlim(-0.00158, 0.00288)
#    pl.ylim( 0.81045, 0.98785)
    pl.show(block=False)
    pl.savefig("../thesis/figures/OtherBisectors.png")
