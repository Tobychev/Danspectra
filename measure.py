import matplotlib.pyplot as pl
import visualize as vis
region = __import__("5053")

lims = {}
lims["ewlim"]   = ( 0.5 , 1.5  )
lims["vellim"]  = (-6 , 7  )
lims["rellim"]  = ( 0.2 , 1.3  )
lims["fw13lim"] = ( -0.1 , 0.9  )
lims["fwhmlim"] = ( -0.1 , 1.3  )
lims["fw23lim"] = (-0.2 , 1.6  )
lims["as13lim"] = (-0.6 , 0.4  )
lims["as12lim"] = (-0.7 , 0.6  )
lims["as23lim"] = (-0.9 , 0.6  )

as2res = {} 
for line in region.as2lines:
    name = line.name.split(" ")[0]
    print(name)
    as2res[name] = line.measure(region.as2)





