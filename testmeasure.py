import spectra as spc
import matplotlib.pyplot as pl
import visualize as vis
import pickle as pic
import numpy as np
import traceback

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

bot,vel,fwhm,as12,fw13,as13,fw23,as23,err,ew,con,smu,sskew,skurt,evar,eske,ekur = range(0,17)



NiI  = {"El":3.8474,"gf":-0.380,"lam":504.8847,"dep":0.580,"name":"Ni I"}
CI   = {"El":7.6848,"gf":-1.303,"lam":505.2144,"dep":0.141,"name":"C I"}
C2   = {"El":0.4891,"gf": 0.158,"lam":505.2620,"dep":0.020,"name":"C_2 "}
TiI  = {"El":2.1747,"gf":-0.270,"lam":505.2868,"dep":0.277,"name":"Ti I"}
Myst = {"El":    -1,"gf":     0,"lam":505.3577,"dep":0.111,"name":"Myst"}
FeI  = {"El":3.6398,"gf":-1.921,"lam":505.4642,"dep":0.476,"name":"Fe I"}
FeI2 = {"El":3.6417,"gf":-2.830,"lam":505.8496,"dep":0.130,"name":"Fe I"}

wNiI  = [1064,1094]
wCI   = [731,761]
wC2   = [690,701]
wTiI  = [666,676]
wMyst = [580,620]
wFeI  = [480,504]
wFeI2 = [87,114]

sf_qu1 = spc.SpectraFactory("data/5053_aS1",framerows=792,framecols=1466)
sf_qu1.frame_row_cut([0]+list(range(666,679))+[791])
sf_qu1.frame_col_cut([0,1,1465])
sf_qu1.contrast_cut(80)
sf_qu1.set_continua("segments")
lowNiIcut = -0.511135 # Read off from Ni I vel
qu1sel = np.arange(0,3107) # Givin pretty pure cut

qu1 = sf_qu1.make_spectra()
qu1.modify(lambda x : x[qu1sel,:]) # Cutting
qu1.meta.cont = (qu1.meta.cont[0][qu1sel], qu1.meta.cont[1][qu1sel]) # Updating continua

qu1m = qu1[:,:].mean(axis=0) 
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1] # Define continua for qu1 series


qu1NiI  = spc.splineline(wNiI,  NiI,  qu1.meta);qu1NiI.recenter(qu1m)
qu1CI   = spc.splineline(wCI,   CI,   qu1.meta);qu1CI.recenter(qu1m)
qu1C2   = spc.splineline(wC2,   C2,   qu1.meta);qu1C2.recenter(qu1m)
qu1TiI  = spc.splineline(wTiI,  TiI,  qu1.meta);qu1TiI.recenter(qu1m)
qu1Myst = spc.splineline(wMyst, Myst, qu1.meta);qu1Myst.recenter(qu1m)
qu1FeI  = spc.splineline(wFeI,  FeI,  qu1.meta);qu1FeI.recenter(qu1m)
qu1FeI2 = spc.splineline(wFeI2, FeI2, qu1.meta);qu1FeI2.recenter(qu1m)
qu1lines = [qu1NiI,qu1CI,qu1C2,qu1TiI,qu1Myst,qu1FeI,qu1FeI2]

regname = "6449"
res = {} 
for i,line in enumerate(qu1lines):
    name ="{}_{}".format(line.name.split()[0],int(line.cent*10))
    try:
        res[name] = line.measure(qu1)
    except Exception as err:
        print("Line {} failed to measure:".format(name))
        print(err)
        traceback.print_exc()

