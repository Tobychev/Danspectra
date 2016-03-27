import matplotlib.pyplot as pl
import visualize as vis
import danframe as dan
import kontin as con
import lines as lin
import numpy as np
import astropy.io.fits as fts
import astropy.stats as ast
import scipy.interpolate as si
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


def makespline(spec,line,group,kns=6):
    lmbd  = group.lmbd[line.idx]
    _,kno = np.histogram(lmbd,kns+2)
    kno   = kno[1:-2]
    return si.LSQUnivariateSpline(lmbd[::-1],spec[::-1],kno)

def measure_spline(spl,line,group,dl=2e-5):
    lmbd = group.lmbd[line.idx]
    lmbd = np.linspace(lmbd[0],lmbd[-1],int( (lmbd[0]-lmbd[-1])/dl )) 
    bot  = spl(lmbd).min()
    icnt = spl(lmbd).argmin()
    cnt  = lmbd[icnt]
    bo12 = (1 +   bot)/2
    bo13 = (1 + 2*bot)/3
    bo23 = (2 +   bot)/3
    fwhm,as12 = __width_assym(spl,lmbd,bo12,cnt)
    fw13,as13 = __width_assym(spl,lmbd,bo13,cnt)
    fw23,as23 = __width_assym(spl,lmbd,bo23,cnt)
    cnt = 299792.458*(cnt-line.cent)/line.cent

    return bot,cnt,fwhm,as12,fw13,as13,fw23,as23,spl.get_residual()

def __width_assym(spl,lmbd,lev,cnt):
    ilev,= np.where(spl(lmbd) <= lev)
    wdth = lmbd[ilev[0]] - lmbd[ilev[-1]]
    assm = cnt  - (lmbd[ilev[0]] + lmbd[ilev[-1]])/2
    return wdth,assm


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
#block = block[cuts.reshape(-1),:]
binNr = block[(sort==nr),:]

splmes = np.zeros(block[:,1:10].shape)

for i,row in enumerate(block):
    mf = makespline(row[line.idx],line,s6405_t5p,9)
    splmes[i,:] = measure_spline(mf,line,s6405_t5p) 

pl.plot(mesmyst[con].reshape(-1),splmes[:,1],'o')
pl.show()
