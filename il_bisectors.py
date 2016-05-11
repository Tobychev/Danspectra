import matplotlib.pyplot as pl
import scipy.interpolate as si
import scipy.optimize as ro
import numpy as np
regnames = ["5053","5215","5654","6405","6449"]

def bisect(lmbd,mnspec,line):
    #find highest points to the left and right of linecenter
    left = mnspec[ np.where(lmbd > line.cent)].argmax()
    right= mnspec[ np.where(lmbd < line.cent)].argmax()
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
            pl.close()
            pl.plot(lmbd,mnspec)
            pl.plot(lmbd,mf(lmbd))
            pl.plot([lmbd[0],lmbd[-1]],[y,y])
            pl.show(block=True)
        bisec[i+1] = (lf+rg)/2 - lin.cent
    return bisec,ys

for regname in regnames:
    rg = __import__(regname)
    for lin in rg.qu1lines:
        xs,ys = bisect(rg.qu1.lmbd[lin.idx],rg.qu1m[lin.idx],lin)
        pl.plot(xs,ys,'-o',label=lin.name)

pl.legend()
pl.show()
