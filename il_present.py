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
