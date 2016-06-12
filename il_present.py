import numpy as np
import spectra as spc
import pickle as pic
import matplotlib.pyplot as pl
import matplotlib.ticker as ticker

d =  __import__("6405")

if False:
    fg,ax = pl.subplots(1)
    idx = range(655,745)
    x, y = d.qu1.lmbd[idx],d.qu1m[idx]
    ax.step(x,y)
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:3.1f}"))   
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Relative intensity")
    fg.savefig("../thesis/presentation/clean_spectra.png")
    
    fg.show()

if False:
    fg,ax = pl.subplots(1)
    idx = range(655,745)
    rndrow = 177 # Chosen by looking for nasty examples
    x, y = d.qu1.lmbd[idx],d.qu1[rndrow,idx]
    ax.step(x,y)
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:3.1f}"))   
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Relative intensity")
    fg.savefig("../thesis/presentation/noise_spectra.png")
    fg.show()

d =  __import__("6405")
rndrow = 177
con =  spc.continua(d.qu1m,d.qu1.lmbd,"top N")
scon =  spc.continua(d.qu1m,d.qu1.lmbd,"segments")

pl.plot(d.qu1.lmbd,d.qu1m)
pl.plot(d.qu1.lmbd[con.idx],d.qu1m[con.idx],'rx')
pl.plot(d.qu1.lmbd[scon.idx],d.qu1m[scon.idx],'ko')
pl.show()

n_con =  spc.continua(d.qu1[rndrow,:],d.qu1.lmbd,"top N")
n_scon =  spc.continua(d.qu1[rndrow,:],d.qu1.lmbd,"segments")

pl.plot(d.qu1.lmbd,d.qu1[rndrow,:])
pl.plot(d.qu1.lmbd[n_con.idx],d.qu1[rndrow,n_con.idx],'rx')
pl.plot(d.qu1.lmbd[n_scon.idx],d.qu1[rndrow,n_scon.idx],'ko')
pl.show()
