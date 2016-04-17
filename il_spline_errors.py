import matplotlib.pyplot as pl
import scipy.stats as st
import spectra as spc
import numpy as np
import errors as er
import copy 

def intervalplot(errors,rng,prop,rngname="",yname="",title=""):
    pl.plot(rng,errors[:,prop,4],'r')
    pl.plot(rng,errors[:,prop,3],'b')
    pl.plot(rng,errors[:,prop,1],'b')
    pl.plot(rng,errors[:,prop,0],'r')

    pl.title(title)
    pl.ylabel(yname)
    pl.xlabel(rngname)
    pl.show()

def mystmutant(x,dist):
    x[:,wFlat[0]:wFlat[1]] = x[:,wFlat[0]:wFlat[1]]*dist
    return x



# For spline measurement
bot,cnt,fwhm,as12,fw13,as13,fw23,as23,err,ew,cont = np.arange(0,11)

sf = spc.SpectraFactory("data/6405_aS1")
sf.frame_row_cut([0,799])
sf.contrast_cut(50)
sf.set_continua("segments")
sp = sf.make_spectra()

# Hand picked flat area
wFlat = [1179,1239]; cFlat = sp.lmbd[range(wFlat[0],wFlat[1]+1)].mean()

def mystmutant(x,dist):
    x[:,wFlat[0]:wFlat[1]] = x[:,wFlat[0]:wFlat[1]]*dist
    return x

# Range of line bottoms and variance.
botrange = np.linspace(0.1,0.9,14)
varrange = np.linspace(0.6,2.5,14)

errors = []
#bot = 0.8
var = 1.5
for bot in botrange:
    norm = st.norm.pdf(np.linspace(-6,6,60),scale=var); norm = norm/norm.max()
    line = 1-bot*norm
    my = copy.deepcopy(sp) 
    my.modify(lambda x: mystmutant(x,line))
    print("Estimating errors for {}".format(bot))
    Flat = spc.splineline(wFlat,cFlat,my.meta)
    mesFlat =  Flat.measure(my)
    errors.append( er.err_spline_mes(mesFlat) )

errors = np.array(errors)

intervalplot(errors,1-botrange,fwhm,"Line bottom","Equivalent width")
