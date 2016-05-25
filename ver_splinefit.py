import matplotlib.ticker as tck
import matplotlib.pyplot as pl
import scipy.interpolate as si
import spectra as spc
import numpy as np

def spline_multiplot(smooths,spec,line,testrows=None):
    cols = 3
    rows = len(smooths)
    fig,ax = pl.subplots(rows,cols,sharex=True,sharey=True)
    nrow = spec[:,:].shape[0]
    lmbd = np.linspace(spec.lmbd[line.idx][0],spec.lmbd[line.idx][-1],1e4)
    we   = np.ones(len(spec.lmbd[line.idx])); we[len(we)*2/5:len(we)*3/5] = 1.5
    bins = len(line.idx)
    if testrows is None:
        rnds = np.random.randint(nrow,size=4) 
    else:
        rnds = testrows   
    i = 0
    for smooth in smooths:

        spl1 = si.UnivariateSpline(spec.lmbd[line.idx][::-1],spec[rnds[0],line.idx][::-1],s=smooth,w=we)
        spl2 = si.UnivariateSpline(spec.lmbd[line.idx][::-1],spec[rnds[1],line.idx][::-1],s=smooth,w=we)
        spl3 = si.UnivariateSpline(spec.lmbd[line.idx][::-1],spec[rnds[2],line.idx][::-1],s=smooth,w=we)
#        spl4 = si.UnivariateSpline(spec.lmbd[line.idx][::-1],spec[rnds[3],line.idx][::-1],s=smooth,w=we)     

        fig.axes[i+0].plot(lmbd,spl1(lmbd));fig.axes[i+0].plot(spec.lmbd[line.idx],spec[rnds[0],line.idx],'.r')
        fig.axes[i+0].set_title("Row {},\n Smooth {:.2e} \n err {:.4e}".format(rnds[0],smooth,spl1.get_residual()/bins))
        fig.axes[i+1].plot(lmbd,spl2(lmbd));fig.axes[i+1].plot(spec.lmbd[line.idx],spec[rnds[1],line.idx],'.r')
        fig.axes[i+1].set_title("Row {},\n Smooth {:.2e} \n err {:.4e}".format(rnds[1],smooth,spl2.get_residual()/bins))
        fig.axes[i+2].plot(lmbd,spl3(lmbd));fig.axes[i+2].plot(spec.lmbd[line.idx],spec[rnds[2],line.idx],'.r')
        fig.axes[i+2].set_title("Row {},\n Smooth {:.2e} \n err {:.4e}".format(rnds[2],smooth,spl3.get_residual()/bins))
#        fig.axes[i+3].plot(lmbd,spl4(lmbd));fig.axes[i+3].plot(spec.lmbd[line.idx],spec[rnds[3],line.idx],'.r')
#        fig.axes[i+3].set_title("Row {}, Smooth {:4e},\n err {:.4e}".format(rnds[3],smooth,spl1.get_residual()/bins))

#        i += 4
        i += 3
    fig.tight_layout()
    return fig

def splineevalplot(spline,spec,lmbd):
    pl.plot(lmbd,spec,'r')
    lmbds = np.linspace(lmbd[0],lmbd[-1],1e4)
    pl.plot(lmbds,spline(lmbds),'k')
    pl.show()


def spline_errplot(smooths,spec,line1,line2,testrows=None):
    cols = 2
    rows = len(smooths)
    fig,ax = pl.subplots(rows,cols)#,sharex=,sharey=True)
    nrow = spec[:,:].shape[0]
    lmbd1 = np.linspace(spec.lmbd[line1.idx][0],spec.lmbd[line1.idx][-1],1e4)
    lmbd2 = np.linspace(spec.lmbd[line2.idx][0],spec.lmbd[line2.idx][-1],1e4)
    we1  = np.ones(len(spec.lmbd[line1.idx])); we1[int(len(we1)*2/5):int(len(we1)*3/5)] = 1.5
    we2  = np.ones(len(spec.lmbd[line2.idx])); we2[int(len(we2)*2/5):int(len(we2)*3/5)] = 1.5
    bins1 = len(line1.idx)
    bins2 = len(line2.idx)
    if testrows is None:
        rnds = np.random.randint(nrow,size=4) 
    else:
        rnds = testrows   
    i = 0
    for smooth in smooths:
        sm1,sm2 = smooth*len(line1.idx),smooth*len(line2.idx)
        spl1 = si.UnivariateSpline(spec.lmbd[line1.idx][::-1],spec[rnds[0],line1.idx][::-1],s=sm1,w=we1)
        spl2 = si.UnivariateSpline(spec.lmbd[line2.idx][::-1],spec[rnds[1],line2.idx][::-1],s=sm1,w=we2)

        fig.axes[i+0].plot(lmbd1,spl1(lmbd1));fig.axes[i+0].plot(spec.lmbd[line1.idx],spec[rnds[0],line1.idx],'.r')
        fig.axes[i+0].set_title("Row {},\n Smooth {:.2e} ".format(rnds[0],smooth,spl1.get_residual()/bins1))
        fig.axes[i+0].xaxis.set_major_formatter(tck.StrMethodFormatter("{x:6.2f}"))
        fig.axes[i+0].set_xlabel("Wavelength")
        fig.axes[i+0].set_ylabel("Relative intensity")

        fig.axes[i+1].plot(lmbd2,spl2(lmbd2));fig.axes[i+1].plot(spec.lmbd[line2.idx],spec[rnds[1],line2.idx],'.r')
        fig.axes[i+1].set_title("Row {},\n Smooth {:.2e} ".format(rnds[1],smooth,spl2.get_residual()/bins2))
        fig.axes[i+1].xaxis.set_major_formatter(tck.StrMethodFormatter("{x:6.2f}"))
        fig.axes[i+1].set_xlabel("Wavelength")
        fig.axes[i+1].set_ylabel("Relative intensity")
        i += 2
    
    return fig

reg = __import__("6405")
smooths = np.array([1e-5, 1.1e-4, 1e-3])

fg = spline_errplot(smooths,reg.qu1,reg.qu1Myst,reg.qu1FeI)
fg.set_size_inches( 10.5   ,  8.1875)
fg.tight_layout()
#pl.show()
fg.savefig("../thesis/figures/splinerr.png")
