import matplotlib.pyplot as pl
import matplotlib.ticker as tck
import matplotlib.cm as cm
import visualize as vis
import numpy as np

regname = "6449"

def spotspec(regname):

    region = __import__(regname)
    spt = region.spt
    con = region.sptcon
    um, = np.where(con < region.um)
    wl, = np.where((region.um <= con ) & (con < region.wl))
    pn, = np.where((region.wl <= con ) & (con < region.pn))
    qs, = np.where(region.pn <= con)

    colr = cm.bwr
    fig,ax = pl.subplots(1)
    ax.plot(np.ones(2)*region.sptMyst.cent,[1.01,1.03],'k',lw=1)
    ax.plot(spt.lmbd,spt[um,:].mean(axis=0),color=colr(0.15),label="Umbra")
    ax.plot(spt.lmbd,spt[wl,:].mean(axis=0),color=colr(0.4), label="Wall")
    ax.plot(spt.lmbd,spt[pn,:].mean(axis=0),color=colr(0.6), label="Penumbra")
    ax.plot(spt.lmbd,spt[qs,:].mean(axis=0),color=colr(0.75),label="Quiet sun")

    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Relative intensity")
    ax.set_xlim(region.xspotlims); 
    ax.set_ylim(region.yspotlims); 
    ax.xaxis.set_major_formatter(tck.StrMethodFormatter("{x:6.2f}"))
    pl.legend(loc="lower left")
    fig.savefig("../thesis/figures/spot{}.png".format(regname))
#    fig.show()

spotspec("5053")
spotspec("5215")
spotspec("5654")
spotspec("6405")
spotspec("6449")
