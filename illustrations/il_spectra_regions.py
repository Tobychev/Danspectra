import matplotlib.pyplot as pl
import matplotlib.cm as cm
import spectra as spc
import numpy as np

r5053 = __import__("5053")
r5215 = __import__("5215")
r5654 = __import__("5654")
r6405 = __import__("6405")
r6449 = __import__("6449")

col = cm.Oranges

def spectraplot(savename,spectra,lines,myst,lims=None):
    fig,ax = pl.subplots(1)
    ax.plot(spectra.lmbd,spectra[:,:].mean(axis=0),'k')
    for line in lines:
        ax.plot(np.ones(2)*line.cent,[1.01,1.04],'k',linewidth=1.)
        ax.plot(spectra.lmbd[line.idx],spectra[:,:].mean(axis=0)[line.idx],'x',color=col(0.45),ms=3.3)
    
    ax.plot(np.ones(1)*myst.cent,1.01,'rv',linewidth=1.3)
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Relative intensity")
    fig.set_size_inches([12.2875, 6.575])
    if lims is not None:
        ax.set_xlim(lims[0])
        ax.set_ylim(lims[1])
   
#    fig.show() 
    fig.savefig(savename)

spectraplot("../../thesis/figures/ex5053.png",r5053.qu1,r5053.qu1lines,r5053.qu1Myst,(
(504.80896484680829,505.89230575757881),(0.051964413375518298, 1.1137618356424892) ))
spectraplot("../../thesis/figures/ex5215.png",r5215.qu1,r5215.qu1lines,r5215.qu1Myst,(
(520.93987222301951, 522.01819558641614),(0.095599448368017897, 1.1728240073047951)))
spectraplot("../../thesis/figures/ex5654.png",r5654.qu1,r5654.qu1lines,r5654.qu2Myst,(
(565.08022182712523, 566.18386748782734),(0.14545504053124736, 1.2099881562449895) ))
spectraplot("../../thesis/figures/ex6405.png",r6405.qu1,r6405.qu1lines,r6405.qu1Myst,(
(639.97118032433366, 640.90599346765202),(0.21664216776069278, 1.1837892132498784)))
spectraplot("../../thesis/figures/ex6449.png",r6449.qu1,r6449.qu1lines,r6449.qu1Myst,(
(644.83444403715976, 645.91897246151677),(0.27096874499276535, 1.1084835818294947)))
