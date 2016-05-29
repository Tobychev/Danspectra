import matplotlib.pyplot as pl
import visualize as vis
import pickle as pic    
import numpy as np

regnames = ["5053","5215","5654","6405","6449"]


fg,ax = pl.subplots(1)

for i,regname in enumerate(regnames):
    dat = np.load("bin/{}_{}.npz".format(regname,"qu1"))
    lines = pic.load(open("bin/{}_{}.lin".format(regname,"qu1"),"rb"))  
    sel = lines[ np.array([lin.dept for lin in lines]).argmin() ]
    sel = "{}_{:4.0f}".format(sel.name.split()[0],int(sel.cent*10))
    vis.kde(dat[sel][:,10],axis=ax)
    ax.lines[-1].set_label(sel)


pl.legend(loc="best")
pl.show()
