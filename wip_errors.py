import matplotlib.pyplot as pl
import scipy.stats as st
import spectra as spc
import numpy as np
import errors as er
import copy 

from imp import reload

reload(spc)

# For moments measurement
#vel,bot,cont,err,ew,mn,var,ske,kur,wvar,wske,wkur = np.arange(0,12)
# For spline measurement
bot,cnt,fwhm,as12,fw13,as13,fw23,as23,err,ew,cont = np.arange(0,11)

wCN   = [289, 342];   cCN   = 640.86681253796701 
wFeI  = [376, 441];   cFeI  = 640.80263785584464
wSiFe = [467, 519];   cSiFe = 640.72889305314993
wMyst = [647, 750];   cMyst = 640.57558287431198
wCNq  = [1010, 1039]; cCNq  = 640.31162454661717

sf = spc.SpectraFactory("data/6405_aS1")
sf.frame_row_cut([0,799])
sf.contrast_cut(50)
sf.set_continua("segments")

sp = sf.make_spectra()
def kde(idx):
    rt =  st.gaussian_kde(mesFlat[:,idx])
    x  = np.linspace(mesFlat[:,idx].min(),mesFlat[:,idx].max(),121)
    pl.plot(x,rt(x))
    pl.show()

def mystmutant(x,dist):
    x[:,wFlat[0]:wFlat[1]] = x[:,wFlat[0]:wFlat[1]]*dist
    return x

def percentiles(measure):
    return ( np.percentile(measure,2.274),
             np.percentile(measure,15.87),
             measure.mean(),
             np.percentile(measure,84.13),
             np.percentile(measure,97.725))

def intervalplot(errors,valrange,prop):
    pl.plot(valrange,errors[:,prop,4],'r')
    pl.plot(valrange,errors[:,prop,3],'b')
    pl.plot(valrange,errors[:,prop,1],'b')
    pl.plot(valrange,errors[:,prop,0],'r')
    pl.show()


# Hand picked flat area
wFlat = [1179,1239]; cFlat = sp.lmbd[range(wFlat[0],wFlat[1]+1)].mean()

# Range of line bottoms.
botrange = np.linspace(0.25,0.9,11)
varrange = np.linspace(0.6,2.5,14)

errors = []
bot = 0.8
for var in varrange:
    norm = st.norm.pdf(np.linspace(-6,6,60),scale=var); norm = norm/norm.max()
    line = 1-bot*norm
    my = copy.deepcopy(sp) 
    my.modify(lambda x: mystmutant(x,line))
    print("Estimating errors for {}".format(bot))
    Flat = spc.splineline(wFlat,cFlat,my.meta)
    mesFlat =  Flat.measure(my)
    errors.append( er.err_spline_mes(mesFlat) )

errors = np.array(errors)
