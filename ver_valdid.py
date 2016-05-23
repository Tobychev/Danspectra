import numpy as np
import matplotlib.pyplot as pl
import matplotlib.ticker as ticker

regnames = ["5053","5215","5654","6405","6449"]

lines = {"CaI":505.35035,"CI":505.35150,"TiI":505.40742,"CoI":521.54823,"NdII":521.565,"Ti I":565.47666,"SiI":640.56873,"Ca I":644.98080}

for regname in regnames:
    rg = __import__(regname)
    fg,ax     = pl.subplots(1)
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:3.3f}"))
   
    ax.plot(rg.qu1.lmbd,rg.qu1m,'.-')
 
    if regname == "5053":
        ax.set_xlim(505.344, 505.361)
        ax.set_ylim(0.885, 0.899)
    if regname == "5215":
        ax.set_xlim(521.549, 521.572)
        ax.set_ylim(0.79289, 0.86898)
    if regname == "5654":
        ax.set_xlim(565.436, 565.461)
        ax.set_ylim(0.86491, 0.90258)
    if regname == "6405":
        ax.set_xlim(640.565, 640.591)
        ax.set_ylim(0.93524, 0.95434)
    if regname == "6449":
        ax.set_xlim(644.907, 644.919)
        ax.set_ylim(0.89524, 0.90507)
    xtic = ax.get_xticks()
    ax.set_xticks(xtic[1:-1])
    ax.set_xlabel("Wavelenght [nm]")
    ax.set_ylabel("Relative intensity")

    
    y = min(ax.get_ylim()); dy = ax.get_ylim()[1] - ax.get_ylim()[0]
    y = y+ dy*np.array([0.9,0.95])
    for lin in lines.keys():

        ax.plot(np.ones(2)*lines[lin],y,'k')

    pl.show(block=True)
