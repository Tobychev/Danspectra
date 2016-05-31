import matplotlib.ticker as ticker
import matplotlib.pyplot as pl
import matplotlib.cm as cm
import numpy as np

regnames =  ["5053","5215","5654","6405","6449"]
for regname in regnames:
    fg,ax     = pl.subplots(1)
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:3.3f}"))
    
    if regname == "5053":
        ax.set_xlim(505.344, 505.361)
        ax.set_ylim(0.885, 0.899)
        num = 3
    if regname == "5215":
        ax.set_xlim(521.549, 521.572)
        ax.set_ylim(0.79289, 0.86898)
        num = 4
    if regname == "5654":
        ax.set_xlim(565.436, 565.461)
        ax.set_ylim(0.86491, 0.90258)
        num = 3
    if regname == "6405":
        ax.set_xlim(640.565, 640.591)
        ax.set_ylim(0.93524, 0.95434)
        num = 3
    if regname == "6449":
        ax.set_xlim(644.907, 644.919)
        ax.set_ylim(0.89524, 0.90507)
        num = 3
    xtic = ax.get_xticks()
    ax.set_xticks(xtic[1:-1])
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Relative intensity")
    
    rg = __import__(regname)
    lmbd,cont = rg.qu1.lmbd,rg.qu1con
    cn, ed    = np.histogram(cont,num)
    scaling   = np.linspace(0.30,0.87,num)[::-1]
    for i in range(0,num):
        sel, = np.where(np.logical_and(ed[i] < cont , cont <= ed[i+1]))
        binmn,contm = rg.qu1[sel,:].mean(axis=0),cont[sel].mean()
        colscale = scaling[i]*contm/ed[-1]
        ax.plot(lmbd,binmn,color=cm.Oranges(colscale),label="{:2.3f}".format(contm))
        print(contm,len(sel))
    ax.legend(loc="best")
    pl.savefig("../thesis/figures/Linecore_{}.png".format(regname))
#    pl.show(block=True)


