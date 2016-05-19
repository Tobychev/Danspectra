import matplotlib.pyplot as pl
import numpy as np
import os

def atlasplot(atlas,region):
    reg = {"5053":19750,"5215":19150,"5654":17650,"6405":15600,"6449":15500}
    qs = atlas["qs{}".format(reg[region])]
    sp = atlas["sp{}".format(reg[region])]
    
    qsnorm = qs[np.where(qs[:,1]> np.percentile(qs[:,1],90) ),1].mean()
    spnorm = sp[np.where(sp[:,1]> np.percentile(sp[:,1],96) ),1].mean()
    pl.plot(1e7/qs[::-1,0],qs[::-1,1]/qsnorm)
    pl.plot(1e7/sp[::-1,0],sp[::-1,1]/spnorm)
    pl.show()

if False:
    caption = "\caption{{Lines in the {reg} region, that extends from {start:.2f} to {finish:.2f} nm.}}\label{{tab:{reg}}}"

    regnames =  ["5053","5215","5654","6405","6449"]
    capdata  = []

    for regname in regnames:
        reg = __import__(regname)
        capdata.append((regname,reg.qu1.lmbd[-1],reg.qu1.lmbd[0]))

    for i in capdata:
        print(caption.format(reg=i[0],start=i[1],finish=i[2]))

datadir = "local_data/"
spotatlas = list(filter(lambda x: x[-4:]=="spot",
                   os.listdir(datadir)) )
spotatlas.sort()

quietatlas = list(filter(lambda x: x[-2:]=="qs",
                   os.listdir(datadir)) )
quietatlas.sort()

atlas = {}

for spot in spotatlas:
    atlas[spot[:-5]] = np.loadtxt(datadir+spot)
    print(spot)

for quiet in quietatlas:
    tmp = np.loadtxt(datadir+quiet)
    atlas["qs"+quiet[2:-3]] = tmp[:,[0,3]]
    print(quiet)
