import visualize as vis
import pickle as pic
import numpy as np

dat = np.load("bin/5053_qu1.npz")
lines=[lin[key] for key in lin.keys()]
lims = vis.make_linemap_lims([dat[key] for key in dat.keys()])


for i,line in enumerate(lines):
    name = line.name.split(" ")[0]
    fig = vis.spline_linemap(res[name],lines[i],lims=lims)
    fig.set_size_inches([ 14.2625,   10.])
#   fig.savefig("../thesis/figures/{}_{}_{}.png".format(regname,typname,name))
    pl.show(block=False)
