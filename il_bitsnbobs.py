import matplotlib.pyplot as pl
import numpy as np
import os

## For making captions to the line tables
if False:
    caption = "\caption{{Lines in the {reg} region, that extends from {start:.2f} to {finish:.2f} nm.}}\label{{tab:{reg}}}"

    regnames =  ["5053","5215","5654","6405","6449"]
    capdata  = []

    for regname in regnames:
        reg = __import__(regname)
        capdata.append((regname,reg.qu1.lmbd[-1],reg.qu1.lmbd[0]))

    for i in capdata:
        print(caption.format(reg=i[0],start=i[1],finish=i[2]))
