import matplotlib.ticker as ticker
import matplotlib.pyplot as pl
import visualize as vis
import pickle as pic
import stats as st
import numpy as np
bot  = 0; vel  = 1; fwhm = 2; as12 = 3; fw13 = 4; as13 = 5; fw23 = 6; as23 = 7; err  = 8; ew   = 9; con  = 10;

regnames = ["5053","5215","5654","6405","6449"]

def dat_plot(datlist,xpos,ypos,normalize=True,epos=None):
    cols = 2
    rows = int(np.ceil(len(datlist)/cols))
    fig,ax = pl.subplots(rows,cols,sharex=True,sharey=True)

    for i,(x,y,name,_,errs) in enumerate(datlist):
        if normalize:
            norm = y.mean()
            fig.axes[i].text(xpos,ypos,"{} {:5.1f}".format(name,norm))
        else:
            norm = 1
            fig.axes[i].text(xpos,ypos,"{}".format(name))

        s2err = np.array([errs[3],errs[0]]).reshape(2,1)
        s1err = np.array([errs[2],errs[1]]).reshape(2,1)

        regx,regy = st.kern_reg(x,y,bins=73)
        fig.axes[i].plot(x,y/norm,'b.',alpha=0.3)
        fig.axes[i].plot(regx,regy/norm,'w',linewidth=3.1)
        fig.axes[i].plot(regx,regy/norm,'r',linewidth=1.2)
        if epos is not None:
            ex,ey = epos
            fig.axes[i].errorbar( ex,ey,yerr=s1err/norm,fmt='b',linewidth=1.5 )
            fig.axes[i].errorbar( ex,ey,yerr=s2err/norm,fmt='r',alpha=0.3 )


    fig.axes[i].set_xlim(0.8,1.27)
    xtic,ytic = fig.axes[i].get_xticks(),fig.axes[i].get_yticks()
    fig.axes[i].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:3.1f}"))
#    fig.axes[i].set_xticks(xtic[1:-1])
#    fig.axes[i].set_yticks(ytic[1:-1])

    return fig


errs = np.load("bin/errlimfile.npz")
errs = {key: errs[key]  for key in errs.keys()}
errs['Co_6454'] = errs['Co_6450'] # Quick fix, here to stay
ewdata = []
lcdata = []
lddata = []
asdata = []
fwdata = []
for region in regnames:
    if region == "6449":
        typ = "qu2"
    else:
        typ = "qu1"
    dat = np.load("bin/{}_{}.npz".format(region,typ))
    lines = pic.load(open("bin/{}_{}.lin".format(region,typ),"rb"))
    def dictgen(i):
        yield lines[i].name.split()[0]+"_"+str(int(lines[i].cent*10))
        yield lines[i]
    lines = dict(map(dictgen,range(len(lines))))

    for line in dat.keys():
        errs[line]
        ewdata.append((dat[line][:,con],dat[line][:,ew]  ,line,lines[line],errs[line][ew,:]))
        lcdata.append((dat[line][:,con],dat[line][:,vel] ,line,lines[line],errs[line][vel,:]))
        lddata.append((dat[line][:,con],dat[line][:,bot] ,line,lines[line],errs[line][bot,:]))
        fwdata.append((dat[line][:,con],dat[line][:,fw13],line,lines[line],errs[line][fw13,:]))
        asdata.append((dat[line][:,con],dat[line][:,as13],line,lines[line],errs[line][as13,:]))

ewdata.sort(key=lambda itm: 1-itm[3].dept)
lcdata.sort(key=lambda itm: 1-itm[3].dept)
lddata.sort(key=lambda itm: 1-itm[3].dept)
fwdata.sort(key=lambda itm: 1-itm[3].dept)
asdata.sort(key=lambda itm: 1-itm[3].dept)



if False:
    data = [ewdata[14],ewdata[19]]
    fig  = dat_plot(sorted(data,key=lambda itm: 1-itm[3].dept,),0.84,1.19,epos=(1.23,1.1))
    fig.set_size_inches([15.925,   6. ])
    fig.tight_layout(h_pad=0.0, w_pad=0.0)
    fig.axes[-1].set_xlim(0.79907330635612739, 1.2690733063561275)
    fig.axes[-1].set_ylim(0.665093813383692, 1.3141211342188408)
    pl.savefig("../thesis/presentation/Mystdouble.png")
    fig.show()

if False:
    data = [ewdata[14],ewdata[15]]
    fig  = dat_plot(sorted(data,key=lambda itm: 1-itm[3].dept,),0.84,1.19,epos=(1.23,1.1))
    fig.set_size_inches([15.925,   6. ])
    fig.tight_layout(h_pad=0.0, w_pad=0.0)
    fig.axes[-1].set_xlim(0.79907330635612739, 1.2690733063561275)
    fig.axes[-1].set_ylim(0.665093813383692, 1.3141211342188408)
    pl.savefig("../thesis/presentation/MystSidouble.png")
    fig.show()

if True:
    data = [ewdata[22],ewdata[26]]
    fig  = dat_plot(sorted(data,key=lambda itm: 1-itm[3].dept,),0.84,1.19,epos=(1.23,1.1))
    fig.set_size_inches([15.925,   6. ])
    fig.tight_layout(h_pad=0.0, w_pad=0.0)
    fig.axes[-1].set_xlim(0.79907330635612739, 1.2690733063561275)
    fig.axes[-1].set_ylim(0.665093813383692, 1.3141211342188408)
    pl.savefig("../thesis/presentation/MystSidouble.png")
    fig.show()
