import matplotlib.pyplot as pl
import matplotlib.ticker as tck
import matplotlib.cm as cmm
import numpy as np
import os

def atlasplot(atlas,regname):
    reg = {"5053":19750,"5215":19150,"5654":17650,"6405":15600,"6449":15500}
    xlim = {"5053":(505.39283866772666, 505.57376665634223),
            "5215":(521.54930858395733, 521.81535307328977),
            "5654":(565.36080224893726, 565.92956831584172),
            "6405":(640.56549492715919, 640.8788013787721 ),
            "6449":(644.92785094602743, 645.34074495254845)}
    ylim = {"5053":(0.63975783972631528, 1.0820081717200596),
            "5215":(0.63457553733797023, 1.1250485137224484),
            "5654":(0.75923171274393061, 1.0552637808945176),
            "6405":(0.86458333333333337, 1.0437500000000002),
            "6449":(0.67393994279930325, 1.0507659350912264)}
    qs = atlas["qs{}".format(reg[regname])]
#    sp = atlas["sp{}".format(reg[regname])]
    lm = atlas["lm{}".format(reg[regname])]
#    region = __import__(regname)
    
    colr = cmm.RdYlGn
    fig,ax = pl.subplots(1)
    ax.plot(1e7/qs[::-1,0],qs[::-1,1],color=colr(0.95),label="Quiet sun")
    ax.plot(lm[:,0]/10,lm[:,1],color=colr(0.70), label="Limb")
#    ax.plot(1e7/sp[::-1,0],sp[::-1,1],color=colr(0.20),label="Spot")

    ax.set_xlabel("Wavelenght [nm]")
    ax.set_ylabel("Relative intensity")
    ax.set_xlim(xlim[regname]); 
    ax.set_ylim(ylim[regname]); 
    ax.xaxis.set_major_formatter(tck.StrMethodFormatter("{x:6.2f}"))
    pl.legend(loc="lower left")
    fig.savefig("../thesis/figures/atlascmp{}.png".format(regname))
    fig.show()


datadir = "local_data/"
spotatlas = list(filter(lambda x: x[-4:]=="spot",
                   os.listdir(datadir)) )
spotatlas.sort()

quietatlas = list(filter(lambda x: x[-2:]=="qs",
                   os.listdir(datadir)) )
quietatlas.sort()

limbatlas = np.loadtxt(datadir+"FTS_atlas_limb.txt",skiprows=2)[:,[0,2]]
limblines = {"5053":[19750,193091,195955,1.394],"5215":[19150,224647,227693,1.435],"5654":[17650,312894,316479,1.548],"6405":[15600,460851,465439,1.776],"6449":[15500,469066,473713,1.77]}


atlas = {}

#for spot in spotatlas:
#    atlas[spot[:-5]] = np.loadtxt(datadir+spot)
#    print(spot)

for quiet in quietatlas:
    tmp = np.loadtxt(datadir+quiet)
    atlas["qs"+quiet[2:-3]] = tmp[:,[0,3]]
    print(quiet)

for key in limblines.keys():
    cm,beg,end,cor = limblines[key]
    tmp = limbatlas[range(beg,end+1),:]
    tmp[:,0] = tmp[:,0]+cor
    atlas["lm"+str(cm)] = tmp

atlasplot(atlas,"5053")
atlasplot(atlas,"5215")
atlasplot(atlas,"5654")
atlasplot(atlas,"6405")
atlasplot(atlas,"6449")
