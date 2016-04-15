import matplotlib.pyplot as pl
import scipy.stats as st
import spectra as spc
import numpy as np
import copy 

from imp import reload

reload(spc)

vel = 0; bot = 1; cont = 2; err = 3; ew  = 4; mn  = 5; var = 6; ske = 7; kur = 8; wvar = 9; wske = 10; wkur = 11;

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

# Hand picked flat area
wFlat = [1179,1239]; cFlat = sp.lmbd[range(wFlat[0],wFlat[1]+1)].mean()


norm = st.norm.pdf(np.linspace(-6,6,60))

simFeI  = 1 - 0.62*norm/norm.max() # Line values eye estimated from atlas
simSiFe = 1 - 0.19*norm/norm.max() # Line values eye estimated from atlas
simMyst = 1 - 0.09*norm/norm.max() # Line values eye estimated from atlas

def mystmutant(x):
    x[:,wFlat[0]:wFlat[1]] = x[:,wFlat[0]:wFlat[1]]*simSiFe
    return x

my = copy.deepcopy(sp) 
my.modify(mystmutant)

Flat = spc.splineline(wFlat,cFlat,my.meta)
mesFlat =  Flat.measure(my)
sigm = np.zeros((11,5))
for i in range(0,11):
    # The percentile ranges are one and two sigmas, calculated by hand from a probability table on wikipedia
    sigm[i,:]= (
        np.percentile(mesFlat[:,i],2.274),
        np.percentile(mesFlat[:,i],15.87),
        mesFlat[:,i].mean(),
        np.percentile(mesFlat[:,i],84.13),
        np.percentile(mesFlat[:,i],97.725))
    m2s,m1s,mu,p1s,p2s = sigm[i,:]
    print("-2s: {: 8.6f}, -1s {: 8.6f}, mu {: 8.6f}, +1s {: 8.6f}, +2s {: 8.6f}".format( m2s/mu,m1s/mu,mu,p1s/mu,p2s/mu ))

def kde(idx):
    rt =  st.gaussian_kde(mesFlat[:,idx])
    x  = np.linspace(mesFlat[:,idx].min(),mesFlat[:,idx].max(),121)
    pl.plot(x,rt(x))
    pl.show()
