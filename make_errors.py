import numpy as np
import spectra as spc
import scipy.interpolate as si 
import matplotlib.pyplot as pl
import numpy.polynomial.chebyshev as ch

def errlims(measure):
    return np.array([ np.percentile(measure,2.274), np.percentile(measure,15.87),
                      np.percentile(measure,84.13), np.percentile(measure,97.725)])-measure.mean()

# For spline measurement
s_bot,s_cnt,s_fwhm,s_as12,s_fw13,s_as13,s_fw23,s_as23,s_err,s_ew,s_cont = np.arange(0,11)

reg = __import__("6405")
flat1 = list(range(169,184))
flat2 = list(range(913,932))
flat3 = list(range(979,995))
flat4 = list(range(1147,1239))

flt1 = reg.qu1m[flat1]/reg.qu1m[flat1].mean() 
flt2 = reg.qu1m[flat2]/reg.qu1m[flat2].mean() 
flt3 = reg.qu1m[flat3]/reg.qu1m[flat3].mean() 
flt4 = reg.qu1m[flat4]/reg.qu1m[flat4].mean() 
flt = np.hstack((flt1,flt2,flt3,flt4))
pos  = np.arange(0,len(flt))
chb  = ch.chebfit(pos,flt,5) # Use for normalization later

bak1 = reg.qu1[:,flat1]/reg.qu1[:,flat1].mean(axis=1).reshape(-1,1)
bak2 = reg.qu1[:,flat2]/reg.qu1[:,flat2].mean(axis=1).reshape(-1,1)
bak3 = reg.qu1[:,flat3]/reg.qu1[:,flat3].mean(axis=1).reshape(-1,1)
bak4 = reg.qu1[:,flat4]/reg.qu1[:,flat4].mean(axis=1).reshape(-1,1)

bak = np.hstack((bak1,bak2,bak3,bak4))/ch.chebval(pos,chb).T

width = 1420 # Magic number found by trial
bak = bak.reshape(-1,width)
blmbd = reg.qu1.lmbd[:width]
bmet  = reg.qu1.meta; bmet.lmbd = blmbd; bmet.ref = None
bmet.cont = ( 0 , np.ones(bak.shape[0]) )


k=0
newlines = []
regnames = ["5053","5215","5654","6405","6449"]
for regname in regnames:
    rg = __import__(regname)
    for lin in rg.qu1lines:
        rawlin = rg.qu1m[lin.idx]
        leng   = len(rawlin)
        sli    = slice(k,leng+k)
        bak[:,sli] = bak[:,sli]*rawlin
        win,mlin,mspec = lin.deconstruct()
        mspec.lmbd = blmbd
        win[0],win[1] = sli.start,sli.stop
        newlines.append(spc.splineline(win,mlin,mspec))
        k += leng+1
bmean = bak.mean(axis=0)
bspec = spc.Spectra("Background+lines",blmbd,bak,bmet)

res = []
for line in newlines:
    line.recenter(bmean)
    res.append(line.measure(bspec))

dump  = {line.name: res[i] for i,line in enumerate(newlines)}
np.savez_compressed("bin/errorfile",**dump)

errs = {}
for key in dump.keys():
    rowerr = [np.abs(errlims(row)) for row in dump[key].T]
    errs[key[:-2]] = rowerr

np.savez_compressed("bin/errlimfile",**errs)

