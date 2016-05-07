import matplotlib.pyplot as pl
import visualize as vis
import pickle as pic
import numpy as np

regnames = ["5053","5215","5654","6405","6449"]


for region in regnames:
    dat = np.load("bin/{}_qu1.npz".format(region))
    lfil = open("bin/{}_qu1.lin".format(region),"rb")
    lines = pic.load(lfil); lfil.close()
    lims = vis.make_linemap_lims([dat[key] for key in dat.keys()])

    pl.close("all")
    print("Plotting region {}".format(region))
    for i,line in enumerate(lines):
        name ="{}_{}".format(line.name.split()[0],int(line.cent*10))
        print("Plotting line {}".format(name))
        try:
            fig = vis.spline_linemap(dat[name],lines[i],lims=lims)
            fig.set_size_inches([  11.075,   8.975])
            fig.subplots_adjust(top=0.94,right=0.96,bottom=0.06,left=0.08,wspace=0.33,hspace=0.43)
            fig.savefig("../thesis/figures/{}_{}_{}.png".format(region,name,"qu1"))
        except KeyError as err:
            print("Line {} not measured".format(line))
            print(err)

pl.close("all")

