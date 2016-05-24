import matplotlib.pyplot as pl
import scipy.optimize as opt
import numpy as np

def gaussian(x, amp, cen, wid):
    return 1-  amp * np.exp(-(x-cen)**2 /wid)

l1,l2 = 641,654
l3,l4 = 587,677

regnames = ["5053","5215","5654","6405","6449"]
rg = __import__("5215")

x = rg.qu1.lmbd[slice(l1,l2)]
y = rg.qu1m[slice(l1,l2)]

wd = 0.011999999/100
cen = 521.52
amp = 1- 0.21

(famp,fcen,fwd),_ = opt.curve_fit(gaussian,x,y,p0=[amp,cen,wd])

xx,yy  = rg.qu1.lmbd[slice(l3,l4)],rg.qu1m[slice(l3,l4)]

#pl.plot(x,y-gaussian(x,famp,fcen,fwd),'.-')
#pl.plot(x,y,'.-')
pl.plot(xx, yy- gaussian(xx,famp,fcen,fwd))
pl.show()

from astropy.modeling import fitting
from astropy.modeling.models import Lorentz1D

modL = Lorentz1D(amplitude=1-amp,x_0=cen,fwhm=wd*100)
res = fitter(modL,xx,1-yy)
xx,yy  = rg.qu1.lmbd[slice(l3,l4)],rg.qu1m[slice(l3,l4)]
pl.plot(xx,1-((1-yy)-res(xx)))
