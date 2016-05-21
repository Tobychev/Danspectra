import matplotlib.pyplot as pl
import matplotlib.ticker as tck
import matplotlib.cm as cmm
import numpy as np
import os

def atlasplot(atlas,regname):
    reg = {"5053":19750,"5215":19150,"5654":17650,"6405":15600,"6449":15500}
    qs = atlas["qs{}".format(reg[regname])]
    sp = atlas["sp{}".format(reg[regname])]
    lm = atlas["lm{}".format(reg[regname])]
    region = __import__(regname)
    
    colr = cmm.Oranges
    fig,ax = pl.subplots(1)
    ax.plot(1e7/qs[::-1,0],qs[::-1,1],color=colr(0.75),label="Quiet sun")
    ax.plot(lm[:,0]/10,lm[:,1],color=colr(0.55), label="Limb")
    ax.plot(1e7/sp[::-1,0],sp[::-1,1],color=colr(0.15),label="Spot")

    ax.set_xlabel("Wavelenght [nm]")
    ax.set_ylabel("Relative intensity")
    ax.set_xlim(region.xspotlims); 
    ax.set_ylim(region.yspotlims); 
    ax.xaxis.set_major_formatter(tck.StrMethodFormatter("{x:6.2f}"))
    pl.legend(loc="lower left")
#    fig.savefig("../thesis/figures/spot{}.png".format(regname))
    fig.show()


datadir = "local_data/"
spotatlas = list(filter(lambda x: x[-4:]=="spot",
                   os.listdir(datadir)) )
spotatlas.sort()

quietatlas = list(filter(lambda x: x[-2:]=="qs",
                   os.listdir(datadir)) )
quietatlas.sort()

limbatlas = np.loadtxt(datadir+"FTS_atlas_limb.txt",skiprows=2)[:,[0,2]]
limblines = {"5053":[19750,193091,195955,0.1394],"5215":[19150,224647,227693,0.1435],"5654":[17650,312894,316479,0.1548],"6405":[15600,460851,465439,0.1776],"6449":[15500,469066,473713,0.177]}


atlas = {}

for spot in spotatlas:
    atlas[spot[:-5]] = np.loadtxt(datadir+spot)
    print(spot)

for quiet in quietatlas:
    tmp = np.loadtxt(datadir+quiet)
    atlas["qs"+quiet[2:-3]] = tmp[:,[0,3]]
    print(quiet)

for key in limblines.keys():
    cm,beg,end,cor = limblines[key]
    tmp = limbatlas[range(beg,end+1),:]
    tmp[:,0] = tmp[:,0]+cor
    atlas["lm"+str(cm)] = tmp

