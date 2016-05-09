import matplotlib.pyplot as pl
import visualize as vis
import pickle as pic
import numpy as np

def kde_multiplot(datlist,order=None):
    cols = 3
    rows = int(np.ceil(len(datlist)/cols))
    fig,ax = pl.subplots(rows,cols)
    if order is None:
        order = range(0,len(datlist))
    for el in list(order):       
        fig.axes[el] = vis.kde(datlist[el],fig.axes[el])
    return fig

def cut_measure(measure,sel):
    return

names = ["bot","vel","fwhm","as12","fw13","as13","fw23","as23","err","ew","con"]
bot,vel,fwhm,as12,fw13,as13,fw23,as23,err,ew,con = list(range(0,11))

regnames = ["5053"] # ["5053","5215","5654","6405","6449"]

for region in regnames:
    raw  = __import__(region)
    raw  = raw.qu1
    dat  = np.load("bin/{}_qu1.npz".format(region))
    lfil = open("bin/{}_qu1.lin".format(region),"rb")
    lines = pic.load(lfil); lfil.close()



lin = dat.keys()[1]
for lin in dat.keys():
    datlist = [dat[lin][:,el] for el in range(0,11)]
    fg = kde_multiplot(datlist)
    fg.suptitle(lin)
    fg.show()

