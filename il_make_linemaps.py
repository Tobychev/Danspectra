import matplotlib.pyplot as pl
import visualize as vis
import pickle as pic
import numpy as np


dat = np.load("bin/6449_qu2.npz")
lines = pic.load(open("bin/6449_qu2.lin","rb"))
lims = vis.make_linemap_lims([dat[key] for key in dat.keys()])

pl.close("all")

for i,line in enumerate(lines):
    name ="{}_{}".format(line.name.split()[0],int(line.cent*10))
    fig = vis.spline_linemap(dat[name],lines[i],lims=lims)
    fig.set_size_inches([ 14.2625,   10.])
#   fig.savefig("../thesis/figures/{}_{}_{}.png".format(regname,typname,name))
    pl.show(block=False)
