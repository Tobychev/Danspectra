import matplotlib.pyplot as pl
import visualize as vis
import numpy as np

regname = "6405"
typname = "spt"

region = __import__(regname)

spec  = region.spt
lines = region.sptlines
lims  = region.sptlims

res = {} 
for i,line in enumerate(lines):
    name = line.name.split(" ")[0]
    res[name] = line.measure(spec)

np.savez("{}-{}".format(regname,typname),res)

for i,line in enumerate(lines):
    name = line.name.split(" ")[0]
    fig = vis.spline_linemap(res[name],lines[i],lims=lims)
    fig.set_size_inches([ 14.2625,   10.])
    fig.savefig("../thesis/figures/{}_{}_{}.png".format(regname,typname,name))
    pl.show(block=False)


