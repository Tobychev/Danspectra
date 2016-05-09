import spectra as spc
import matplotlib.pyplot as pl
import visualize as vis
import pickle as pic
import numpy as np
import traceback

wMyst = [1211,1327]; cMyst = 644.9127
wCa   = [1150,1198]; cCa   = 644.9820
wBlnd = [1088,1151]; cBlnd = 645.0179
wSiI  = [834,889];   cSiI  = 645.2315
wH2O  = [613,637];   cH2O  = 645.4139
wCo   = [449,534];   cCo   = 645.5001
wH2O2 = [8,34];      cH2O2 = 645.8892

sf_qu1 = spc.SpectraFactory("data/6449_aS1",framerows=802,framecols=1514)
sf_qu1.frame_row_cut([0]+list(range(668,674))+[801])
sf_qu1.frame_col_cut([0,1513])
sf_qu1.contrast_cut(50)
sf_qu1.set_continua("segments")

qu1 = sf_qu1.make_spectra()
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1] # Define continua for qu1 series

qu1Myst = spc.testspline(wMyst, cMyst, qu1.meta,"Myst"   ) # Define lines for qu1 series...
qu1Ca   = spc.testspline(wCa,   cCa,   qu1.meta,"CaI"    )
qu1Blnd = spc.testspline(wBlnd, cBlnd, qu1.meta,"Blend" )
qu1SiI  = spc.testspline(wSiI,  cSiI,  qu1.meta,"SiI" )
qu1H2O  = spc.testspline(wH2O,  cH2O,  qu1.meta,"tel H2O")
qu1Co   = spc.testspline(wCo,   cCo,   qu1.meta,"CoI"    )
qu1H2O2 = spc.testspline(wH2O2, cH2O2, qu1.meta,"tel H2O")
qu1lines = [qu1Myst,qu1Ca,qu1Blnd,qu1SiI,qu1H2O,qu1Co,qu1H2O2]

regname = "6449"
res = {} 
for i,line in enumerate(qu1lines):
    name ="{}_{}".format(line.name.split()[0],int(line.cent*10))
    try:
        res[name] = line.measure(qu1)
    except Exception as err:
        print("Line {} failed to measure:".format(name))
        print(err)
        traceback.print_exc()
np.savez_compressed("bin/test_{}_qu1".format(regname),**res)
pic.dump(qu1lines,open("bin/test_{}_qu1.lin".format(regname),"wb"))

dat = np.load("bin/test_6449_qu1.npz")
lines = pic.load(open("bin/test_6449_qu1.lin","rb"))
lims = vis.make_linemap_lims([dat[key] for key in dat.keys()])

pl.close("all")

for i,line in enumerate(lines):
    name ="{}_{}".format(line.name.split()[0],int(line.cent*10))
    fig = vis.spline_linemap(dat[name],lines[i],lims=lims)
    fig.set_size_inches([ 14.2625,   10.])
#   fig.savefig("../thesis/figures/{}_{}_{}.png".format(regname,typname,name))
    pl.show(block=False)
