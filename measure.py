import matplotlib.pyplot as pl
import visualize as vis
region = __import__("6449")

regname = "6449"
typname = "qu1"
spec  = region.qu1
lines = region.qu1lines
lims = region.qu1lims

res = {} 
for i,line in enumerate(lines):
    name = line.name.split(" ")[0]
    res[name] = line.measure(spec)

for i,line in enumerate(lines):
    name = line.name.split(" ")[0]
    fig = vis.spline_linemap(res[name],lines[i],lims=lims)
    fig.set_size_inches([ 14.2625,   10.])
#    fig.savefig("../thesis/figures/{}_{}_{}.png".format(regname,typname,name))
    pl.show(block=False)


