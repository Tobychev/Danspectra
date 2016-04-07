import matplotlib.pyplot as pl
import interactive as intr
import visualize as vis
import danframe as dan
import kontin as con
import lines as lin
import numpy as np
import imp
import astropy.io.fits as ft
import numpy.random as rnd
import binner as bn
import scipy.stats as st

imp.reload(lin); imp.reload(vis);imp.reload(dan);imp.reload(lin);imp.reload(intr);imp.reload(bn)

s6405_t5p = dan.frameseries("data/6405_aS1","top 5%")
count_name = str(s6405_t5p.wave)+"_"+str(s6405_t5p.series)+"__meanspec.fits"
counts,  = ft.open("data/"+count_name)

# Load continuum contrast 
cont, =  ft.open("data/6405_aS1__concont.fits")

# Group cont with order value, sort by cont in decreasing order, then save only the top twelve frames
qual = [(x,i) for i,x in enumerate(cont.data)] ; qual.sort(key=lambda x: x[0],reverse=True)
s6405_t5p.frames = [s6405_t5p.frames[x[1]] for x in qual[:12]]
s6405_t5p.normalize()

# Loading line definitions and creating moment and spline lines
mCN,mFeI,mSiFe,mmyst,mCNq = lin.make_lines_from_wins(s6405_t5p,s6405_t5p.pkwindows)
CN,FeI,SiFe,myst,CNq      = lin.make_splines_from_wins(s6405_t5p,s6405_t5p.pkwindows)

# Making a big block of spectral data across all good frames
block = s6405_t5p.frames[0].data
for frm in s6405_t5p.frames[1:]:
    block = np.vstack((block,frm.data))

# Fitting a continuum to the meanspectra
top20   = np.argsort(counts.data)[-20:]
con_pol = np.poly1d(np.polyfit(s6405_t5p.lmbd[top20],counts.data[top20],deg=1))

# Convenience variables
lmMyst = s6405_t5p.lmbd[myst.idx]
reps    = 10000 # Number of poisson numbers

# Simulating spectra by drawing poisson numbers of apropriate mean for each pixel
# Fails to make convincing synthetic spectra, errors seem far too big

#cnMyst  = counts.data[myst.idx].astype(int)
#mcblock = rnd.poisson(cnMyst,(reps,len(cnMyst)))
#con     = np.ones(reps)*con_pol(myst.cent)

# Spline measuring on the simulated data
#mes    = myst.measure_on_block(s6405_t5p,mcblock/con_pol(lmMyst),con)

mesmyst = myst.measure(s6405_t5p)
binMyst = bn.binned_framegroup(myst,s6405_t5p,mesmyst[:,10])
bincent = np.convolve(binMyst.bins,np.array([0.5,0.5]),mode="valid")
sel, = np.where(binMyst.counts > 30)
sorting = np.digitize(mesmyst[:,10],binMyst.bins[:-1]) 

var = np.zeros( (len(sel),len(myst.idx)) )
men = np.zeros( (len(sel),len(myst.idx)) )
kur = np.zeros( (len(sel),len(myst.idx)) )

for i,idx in enumerate(sel):
    subs = block[ sorting == (idx+1), : ]
    subs = subs[:,myst.idx]
    print(subs.shape)
    for j,row in enumerate(subs.T):
        men[i,j] = row.mean()
        var[i,j] = (row-men[i,j]).std()/np.sqrt(subs.shape[0])
        kur[i,j] = st.kurtosis(row-men[i,j])

con     = np.ones(reps)#*con_pol(myst.cent)
for binNr,bn in enumerate(sel):
    mcblock = rnd.normal( loc=men[binNr,:] , scale=np.sqrt(var[binNr,:]) ,size=(reps,len(myst.idx)) )
    mes     = myst.measure_on_block(s6405_t5p,mcblock,con)
    print("")
    print("For bin {:.4f} with {} spectra".format(bincent[binNr],binMyst.counts[bn]))
    for i in range(0,10):
        # The percentile ranges are one and two sigmas, calculated by hand from a probability table on wikipedia
        print("-2s: {: 8.6f}, -1s {: 8.6f}, mu {: 8.6f}, +1s {: 8.6f}, +2s {: 8.6f}".format(
            np.percentile(mes[:,i],2.274),
            np.percentile(mes[:,i],15.87),
            mes[:,i].mean(),
            np.percentile(mes[:,i],84.13),
            np.percentile(mes[:,i],97.725)))
