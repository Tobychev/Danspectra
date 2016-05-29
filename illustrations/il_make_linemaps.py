import matplotlib.pyplot as pl
import visualize as vis
import pickle as pic
import numpy as np

regnames = ["5053","5215","5654","6405","6449"]

errs = np.load("bin/errlimfile.npz")

for region in regnames:
    for typ in ["qu1","qu2","spt"]:
        dat = np.load("bin/{}_{}.npz".format(region,typ))
        lfil = open("bin/{}_{}.lin".format(region,typ),"rb")
        lines = pic.load(lfil); lfil.close()
        lims = vis.make_linemap_lims([dat[key] for key in dat.keys()])

        pl.close("all")
        print("Plotting region {}".format(region))
        for i,line in enumerate(lines):
            name ="{}_{}".format(line.name.split()[0],int(line.cent*10))
            print("Plotting line {}".format(name))        
            try:
                fig = vis.spline_linemap(dat[name],lines[i],lims=lims)
                fig = vis.add_errs_linemap(fig,errs[line.name[:-2]],dat[name][:,9].mean())
                fig.set_size_inches([12.2375,   8.7125])
                fig.subplots_adjust(top=0.94,right=0.96,bottom=0.06,left=0.08,wspace=0.33,hspace=0.26)
                fig.savefig("../thesis/figures/{}_{}_{}.png".format(region,name,typ))
            except KeyError as err:
                print("Line {} not plotted:".format(line))
                print(err)

pl.close("all")

