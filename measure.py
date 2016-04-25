import matplotlib.pyplot as pl
import visualize as vis
region = __import__("6449")

lims = {}
lims["ewlim"]   = ( -0.7 , 3.6  )
lims["vellim"]  = (-8 , 10  )
lims["rellim"]  = ( -1 , 6  )

lims["fw13lim"] = ( - 0.1, 1.1  )
lims["fwhmlim"] = (  -0.2 , 1.1 )
lims["fw23lim"] = (-0.2 , 1.1    )

lims["as13lim"] = (-1.4 , 0.4  )
lims["as12lim"] = (-1.5 , 0.5  )
lims["as23lim"] = (-1.3 , 0.4  )

res = {} 
regname = "6449"
typname = "qu"
spec  = region.qu1
lines = region.qu1lines
lims  = region.qu1lims


for i,line in enumerate(lines):
    name = line.name.split(" ")[0]
    res[name] = line.measure(spec)

for i,line in enumerate(lines):
    fig = vis.spline_linemap(res[name],lines[i],lims=lims)
    fig.set_size_inches([ 14.2625,   9.7625])
    pl.show(block=False)
#    fig.savefig("../thesis/figures/{}_{}_{}.png".format(regname,typname,name))


