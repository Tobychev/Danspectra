import matplotlib.pyplot as pl
import scipy.interpolate as si
import spectra as spc
import numpy as np

wNiI  = [1064,1094] 
wCI   = [731,761];  
wC2   = [690,701];  
wTiI  = [666,676];  
wMyst = [580,620];  
wFeI  = [480,504];  
wFeI2 = [87,114];   

C2test = [335,1341,1616,1670]
MystTest = [3065,1140,2719,1536,2828] 
TiTest = [1996,2645,2004,419]
FeTest = [1004,1602,1881,1480]

reg = __import__("5053")
C2   = spc.testspline(wC2,   reg.C2,   reg.qu1.meta)
Myst   = spc.testspline(wMyst,   reg.Myst,   reg.qu1.meta)

def spline_multiplot(smooths,spec,line,testrows=None):
    cols = 4
    rows = len(smooths)
    fig,ax = pl.subplots(rows,cols)
    nrow = spec[:,:].shape[0]
    lmbd = np.linspace(spec.lmbd[line.idx][0],spec.lmbd[line.idx][-1],1e4)
    we   = np.ones(len(spec.lmbd[line.idx])); we[len(we)*2/5:len(we)*3/5] = 1.2
    print(we)
    bins = len(line.idx)
    i = 0
    for smooth in smooths:
        if testrows is None:
            rnds = np.random.randint(nrow,size=4) 
        else:
            rnds = testrows
        spl1 = si.UnivariateSpline(spec.lmbd[line.idx][::-1],spec[rnds[0],line.idx][::-1],s=smooth,w=we)
        spl2 = si.UnivariateSpline(spec.lmbd[line.idx][::-1],spec[rnds[1],line.idx][::-1],s=smooth,w=we)
        spl3 = si.UnivariateSpline(spec.lmbd[line.idx][::-1],spec[rnds[2],line.idx][::-1],s=smooth,w=we)
        spl4 = si.UnivariateSpline(spec.lmbd[line.idx][::-1],spec[rnds[3],line.idx][::-1],s=smooth,w=we)     

        fig.axes[i+0].plot(lmbd,spl1(lmbd));fig.axes[i+0].plot(spec.lmbd[line.idx],spec[rnds[0],line.idx],'r')
        fig.axes[i+0].set_title("Row {}, Smooth {:4e},\n err {:.4e}".format(rnds[0],smooth,spl1.get_residual()/bins))
        fig.axes[i+1].plot(lmbd,spl2(lmbd));fig.axes[i+1].plot(spec.lmbd[line.idx],spec[rnds[1],line.idx],'r')
        fig.axes[i+1].set_title("Row {}, Smooth {:4e},\n err {:.4e}".format(rnds[1],smooth,spl1.get_residual()/bins))
        fig.axes[i+2].plot(lmbd,spl3(lmbd));fig.axes[i+2].plot(spec.lmbd[line.idx],spec[rnds[2],line.idx],'r')
        fig.axes[i+2].set_title("Row {}, Smooth {:4e},\n err {:.4e}".format(rnds[2],smooth,spl1.get_residual()/bins))
        fig.axes[i+3].plot(lmbd,spl4(lmbd));fig.axes[i+3].plot(spec.lmbd[line.idx],spec[rnds[3],line.idx],'r')
        fig.axes[i+3].set_title("Row {}, Smooth {:4e},\n err {:.4e}".format(rnds[3],smooth,spl1.get_residual()/bins))

        i += 4
    fig.tight_layout()
    return fig

def splineevalplot(spline,spec,lmbd):
    pl.plot(lmbd,spec,'r')
    lmbds = np.linspace(lmbd[0],lmbd[-1],1e4)
    pl.plot(lmbds,spline(lmbds),'k')
    pl.show()

aspl = si.UnivariateSpline(reg.qu1.lmbd[C2.idx][::-1],reg.qu1[0,C2.idx][::-1],s=1)
mspl = si.UnivariateSpline(reg.qu1.lmbd[Myst.idx][::-1],reg.qu1[0,Myst.idx][::-1],s=1)
#C2.measure(reg.qu1)

