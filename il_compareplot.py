import matplotlib.ticker as ticker
import matplotlib.pyplot as pl
import visualize as vis
import pickle as pic
import stats as st
import numpy as np

bot  = 0; vel  = 1; fwhm = 2; as12 = 3; fw13 = 4; as13 = 5; fw23 = 6; as23 = 7; err  = 8; ew   = 9; con  = 10;

def comp_plot(datlist,xpos,ypos,normalize=True):
    cols = 4
    rows = int(np.ceil(len(datlist)/cols))
    fig,ax = pl.subplots(rows,cols,sharex=True,sharey=True)

    for i,(x,y,name,_) in enumerate(datlist):
        if normalize:
            norm = y.mean()
            fig.axes[i].text(xpos,ypos,"{} {:5.1f}".format(name[:-5],norm))
        else:
            norm = 1
            fig.axes[i].text(xpos,ypos,"{}".format(name))

        regx,regy = st.kern_reg(x,y,bins=73)
        fig.axes[i].plot(x,y/norm,'bo',alpha=0.1)
        fig.axes[i].plot(regx,regy/norm,'w',linewidth=2.1)
        fig.axes[i].plot(regx,regy/norm,'r',linewidth=1.2)

    fig.axes[i].set_xlim(0.8,1.27)
    xtic,ytic = fig.axes[i].get_xticks(),fig.axes[i].get_yticks()
    fig.axes[i].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:3.1f}"))
    fig.axes[i].set_xticks(xtic[1:-1])
    fig.axes[i].set_yticks(ytic[1:-1])

    return fig

regnames = ["5053","5215","5654","6405","6449"]

ewdata = []
lcdata = []
for region in regnames:
    dat = np.load("bin/{}_qu1.npz".format(region))
    lines = pic.load(open("bin/{}_qu1.lin".format(region),"rb"))
    def dictgen(i):
        yield lines[i].name.split()[0]+"_"+str(int(lines[i].cent*10))
        yield lines[i]
    lines = dict(map(dictgen,range(len(lines))))

    for line in dat.keys():
        ewdata.append((dat[line][:,con],dat[line][:,ew],line,lines[line]))
        lcdata.append((dat[line][:,con],dat[line][:,vel],line,lines[line]))

        

#fig = comp_plot(sorted(ewdata,key=lambda itm: itm[3].dept),0.83,1.62)
fig = comp_plot(sorted(lcdata,key=lambda itm: itm[3].dept),0.83,1.62,False)
fig.set_size_inches([ 10.125 ,   9.7625])
fig.tight_layout(h_pad=0.0, w_pad=0.0)
fig.show()

