import matplotlib.pyplot as pl
import spectra as spc
import numpy as np

r5053 = __import__("5053")
r5215 = __import__("5215")
r5654 = __import__("5654")
r6405 = __import__("6405")
r6449 = __import__("6449")

def spectraplot(savename,spectra,lines,myst):
    fig,ax = pl.subplots(1)
    ax.plot(spectra.lmbd,spectra[:,:].mean(axis=0),'k')
    for line in lines:
        ax.plot(np.ones(2)*line.cent,[1.01,1.04],'k',linewidth=1.)
    
    ax.plot(np.ones(2)*myst.cent,[1.01,1.04],'r',linewidth=1.3)
    ax.set_xlabel("Wavelenght [nm]")
    ax.set_ylabel("Relative intensity")
    fig.set_size_inches([12.2875, 6.575])
    fig.savefig(savename)

spectraplot("../../thesis/figures/ex5053.png",r5053.as1,r5053.as1lines,r5053.as1Myst)
spectraplot("../../thesis/figures/ex5215.png",r5215.as1,r5215.as1lines,r5215.as1Myst)
spectraplot("../../thesis/figures/ex5654.png",r5654.as2,r5654.as2lines,r5654.as2Myst)
spectraplot("../../thesis/figures/ex6405.png",r6405.as1,r6405.as1lines,r6405.as1Myst)
spectraplot("../../thesis/figures/ex6449.png",r6449.as1,r6449.as1lines,r6449.as1Myst)
