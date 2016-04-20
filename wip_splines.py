import matplotlib.pyplot as pl
import visualize as vis
import danframe as dan
import kontin as con
import lines as lin
import binner as bn
import numpy as np
import astropy.io.fits as fts
import astropy.stats as ast
import scipy.interpolate as si
import imp
imp.reload(lin); imp.reload(dan); imp.reload(con); imp.reload(vis)

s6405_t5p = dan.frameseries("data/6405_aS1","top 5%")

# Load continuum contrast 
cont, =  fts.open("data/6405_aS1__concont.fits")

# Group cont with order value, sort by cont in decreasing order, then save only the top twelve frames
qual = [(x,i) for i,x in enumerate(cont.data)] ; qual.sort(key=lambda x: x[0],reverse=True)
s6405_t5p.frames = [s6405_t5p.frames[x[1]] for x in qual[:12]]
s6405_t5p.normalize()
lmbd = s6405_t5p.lmbd

block = s6405_t5p.frames[0].data

for frm in s6405_t5p.frames[1:]:
    block = np.vstack((block,frm.data))

CN, FeI, SiFe, myst, CNq = lin.make_lines_from_wins(s6405_t5p,s6405_t5p.pkwindows)

Mys, = lin.make_lines_from_wins(s6405_t5p,[ [646, 722]])

#Helpful constants
vel = 0; bot = 1; con = 2; err = 3; ew  = 4; mn  = 5; var = 6; ske = 7; kur = 8
#mesFeI  = FeI.measure(s6405_t5p)
#mesSiFe = SiFe.measure(s6405_t5p)
#mesmyst = myst.measure(s6405_t5p)
#mesCN   = CN.measure(s6405_t5p)
#mesCNq  = CNq.measure(s6405_t5p)

def measure(x,tck,dl=1e-7):
    smallstep = 1e-7
    numsmallstep = 1e3
    lmbd = np.linspace(x[0],x[-1],int( (x[-1]-x[0])/dl ))
    #Do tow rounds to get better acc
    icnt = lmbd[si.splev(lmbd,tck).argmin()]
    botl = np.linspace(icnt*(1-smallstep),icnt*(1+smallstep),int(numsmallstep))
    bot  = si.splev(botl,tck).min()
    cnt  = botl[si.splev(botl,tck).argmin()]
    bo12 = (1 +   bot)/2
    bo13 = (1 + 2*bot)/3
    bo23 = (2 +   bot)/3
    tck13 = (tck[0],tck[1]-bo13,tck[2])
    tck12 = (tck[0],tck[1]-bo12,tck[2])
    tck23 = (tck[0],tck[1]-bo23,tck[2])
    dl13 = si.sproot(tck13)
    dl12 = si.sproot(tck12)
    dl23 = si.sproot(tck23)
    wd13  = dl13[1] - dl13[0]
    wd12  = dl12[1] - dl12[0]
    wd23  = dl23[1] - dl23[0]
    print(wd13,wd12,wd23)
#    assm  = cnt  - (lmbd1 + lmbd2)/2



for i in range(1,2):
    row = np.random.randint(0,len(block))
    x = lmbd[Mys.idx][::-1]
    y = block[row,Mys.idx][::-1]
    d =  1-y/y.max()
    w = d/d.sum()
    #Spl = si.UnivariateSpline(x,y,w=w,s=8e-3)
    tck  = si.splrep(x,y,s=8e-3)
    measure(x,tck)
#   Spl  = lambda x: si.splev(x,tck)

#   pl.step(lmbd[Mys.idx],block[row,Mys.idx],'.',where="mid")#,label="Residual: {:.4e}".format( () ))
#   pl.plot(x,Spl(x))
#   pl.legend()
#   pl.show()
