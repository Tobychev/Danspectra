import numpy as np
import pickle as pic
import visualize as vis
import matplotlib.pyplot as pl
import scipy.interpolate as si


bot  = 0; vel  = 1; fwhm = 2; as12 = 3; fw13 = 4; as13 = 5; fw23 = 6; as23 = 7; err  = 8; ew   = 9; con  = 10;
reg = "6449"
lne = 'Co_6455'
rg = __import__(reg)
dat = np.load("bin/{}_qu1.npz".format(reg))
lines = pic.load(open("bin/{}_qu1.lin".format(reg),"rb"))

#vis.kde(dat[lne][:,vel])


upsel, = np.where(dat[lne][:,vel] > 10)
dwsel, = np.where(dat[lne][:,vel] < 10)

Co = rg.qu1Co
lmbd = rg.qu1.lmbd[rg.qu1Co.idx]
spec = rg.qu1[upsel[0],rg.qu1Co.idx]
spl  = si.UnivariateSpline(lmbd[::-1],spec[::-1],s=1)
vis.splineevalplot(spl,spec,lmbd)


if False:
    pl.plot(rg.qu1.lmbd,rg.qu1[upsel[np.random.randint(len(upsel))],:],'r')
    pl.plot(rg.qu1.lmbd,rg.qu1[dwsel[np.random.randint(len(dwsel))],:],'b')
    pl.show()

