import matplotlib.pyplot as pl
import visualize as vis
import danframe as dan
import kontin as con
import lines as lin
import numpy as np
import astropy.io.fits as fts
import astropy.stats as ast
import imp
imp.reload(lin); imp.reload(dan); imp.reload(con); imp.reload(vis)

def moments(x,lines):
    rows = len(lines)    
    mu, mu2, mu3, mu4 = np.zeros(rows),np.zeros(rows),np.zeros(rows),np.zeros(rows)
    for i in range(0,rows):
        line = lines[i]
        dpdf = (1-line/line.max())
        dpdf = dpdf/dpdf.sum()
        mu[i] = np.sum(dpdf*x) # Reshaping enables broadcasting
        mu2[i]= np.sum(dpdf*(x-mu[i])**2)
        mu3[i]= np.sum(dpdf*(x-mu[i])**3)
        mu4[i]= np.sum(dpdf*(x-mu[i])**4)

    skew = mu3/mu2**(3/2) 
    kurt = (mu4/mu2**2 - 3)
    return mu,mu2,skew,kurt

s6405_t5p = dan.frameseries("data/6405_aS1","top 5%")

# Load continuum contrast 
cont, =  fts.open("data/6405_aS1__concont.fits")

# Group cont with order value, sort by cont in decreasing order, then save only the top twelve frames
qual = [(x,i) for i,x in enumerate(cont.data)] ; qual.sort(key=lambda x: x[0],reverse=True)
s6405_t5p.frames = [s6405_t5p.frames[x[1]] for x in qual[:12]]
s6405_t5p.normalize()

block = s6405_t5p.frames[0].data

for frm in s6405_t5p.frames[1:]:
    block = np.vstack((block,frm.data))

CN, FeI, SiFe, myst, CNq = lin.make_lines_from_wins(s6405_t5p,s6405_t5p.pkwindows)

#Helpful constants
vel = 0; bot = 1; con = 2; err = 3; ew  = 4; mn  = 5; var = 6; ske = 7; kur = 8
mesFeI  = FeI.measure(s6405_t5p)
mesSiFe = SiFe.measure(s6405_t5p)
mesmyst = myst.measure(s6405_t5p)
mesCN   = CN.measure(s6405_t5p)
mesCNq  = CNq.measure(s6405_t5p)

# Set line to test on
if False:
    line = FeI
    step = 6 #FeI number
    cuts  = mesFeI[err] < np.percentile(mesFeI[err],95)    
    quant = mesFeI[con].reshape(-1)
    binFeI    = lin.binned_framegroup(FeI,s6405_t5p,quant,cuts)
    mesBinFeI = binFeI.measure()
    quant  = quant[cuts.reshape(-1)].reshape(-1)
    sort   = binFeI.partition_data(quant)
    binned = binFeI
    

if True:
    line  = myst    
    step  = 8 #myst number
    cuts  = mesmyst[err] < np.percentile(mesmyst[err],95)
    quant = mesmyst[con].reshape(-1)
    binMyst     = lin.binned_framegroup(myst,s6405_t5p,quant,cuts)
    mesBinMyst = binMyst.measure()
    quant  = quant[cuts.reshape(-1)].reshape(-1)
    sort   = binMyst.partition_data(quant)
    binned = binMyst


nr = 1
lmbd  = s6405_t5p.lmbd[line.idx]
block = block[cuts.reshape(-1),:]
binNr = block[(sort==nr),:]
