import matplotlib.pyplot as pl
import visualize as vis
import danframe as dan
import kontin as con
import lines as lin
import numpy as np
import astropy.io.fits as fts

pl.rcParams["figure.figsize"] = (10,6) # Bigger figures
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
mesFeI   = FeI.measure(s6405_t5p)
mesSiFe  = SiFe.measure(s6405_t5p)
mesmyst  = myst.measure(s6405_t5p)
mesCN    = CN.measure(s6405_t5p)
mesCNq   = CNq.measure(s6405_t5p)

quant = mesFeI[con].reshape(-1)
cuts  = mesFeI[err] < np.percentile(mesFeI[err],89)
binFeI = lin.binned_framegroup(FeI,s6405_t5p,quant,cuts)
mesBinnFeI = binFeI.measure()

quant = quant[cuts.reshape(-1)].reshape(-1)

sort = binFeI.partition_data(quant)
binnr = np.unique(sort)

step = 6
i    = 0
lmbd = s6405_t5p.lmbd[FeI.idx]

for nr in binnr[:]:
    print("Events: ",
            binFeI.counts[nr-1],
            (sort == nr).sum(),
            ((binFeI.bins[i] <= quant ) & (quant < binFeI.bins[i+1] ) ).sum()   )
    i+=1
