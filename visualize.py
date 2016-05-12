import matplotlib.pyplot as pl
import scipy.stats as sta
import stats as st
import numpy as np

def splineevalplot(spline,spec,lmbd):
    pl.plot(lmbd,spec,'r')
    lmbds = np.linspace(lmbd[0],lmbd[-1],1e4)
    pl.plot(lmbds,spline(lmbds),'k')
    pl.show()

def kde(measure,axis=None,norm=False):
    rt = sta.gaussian_kde(measure)
    x  = np.linspace(measure.min(),measure.max(),121)
    if norm:
        norm = rt(x).sum()
    else:
        norm = 1
    if axis is not None:
        axis.plot(x,rt(x)/norm)
        return axis
    else:
        pl.plot(x,rt(x)/norm)
        pl.show()

def kde_multiplot(datlist,order=None):
    cols = 3
    rows = int(np.ceil(len(datlist)/cols))
    fig,ax = pl.subplots(rows,cols)
    if order is None:
        order = range(0,len(datlist))
    for el in list(order):       
        fig.axes[el] = kde(datlist[el],fig.axes[el])
    return fig

def make_linemap_lims(measure):    
    props = [ "bot", "vel", "fwhm", "as12", "fw13", "as13", "fw23", "as23", "err", "ew", "con"]
    proplims = {}
    lims = {}
    for i,name in enumerate(props):
        proplims[name] = []
        for line in measure:
            # Argsort returns the index, pick third lowest and third highest value as limits
            if name == "ew":
                mew = line[:,i].mean()
                proplims[name].append(line[ line[:,i].argsort()[[3,-2]],i ]/mew)
            else:
                proplims[name].append(line[ line[:,i].argsort()[[3,-2]],i ])

        most  = np.max(proplims[name]); least = np.min(proplims[name])
#        print(least,most)
        dx = (most-least)*5e-2
        if name == "bot":
            name = "rel"
            
        lims["{}lim".format(name)] = (least-dx,most+dx)
   
    return lims

def spline_linemap(measure,line,mesbin=None,lims=None,errs=None,regbins=73):
    bot  = 0; vel  = 1; fwhm = 2; as12 = 3; fw13 = 4; as13 = 5; fw23 = 6; as23 = 7; err  = 8; ew   = 9; con  = 10;
    if lims is None:
        ewlim   = ( 0.3 , 1.8  )
        vellim  = (-5.8 , 6.1  )
        rellim  = ( 0.2 , 1.3  )
        fw13lim = ( 0.0 , 1.0  )
        fwhmlim = ( 0.0 , 1.0  )
        fw23lim = ( 0.0 , 1.0  )
        as13lim = (-0.2 , 0.2  )
        as12lim = (-0.2 , 0.2  )
        as23lim = (-0.1 , 0.1  )
    else:
        ewlim   = lims["ewlim"]
        vellim  = lims["vellim"]
        rellim  = lims["rellim"]
        fw13lim = lims["fw13lim"]
        fwhmlim = lims["fwhmlim"]
        fw23lim = lims["fw23lim"]
        as13lim = lims["as13lim"]
        as12lim = lims["as12lim"]
        as23lim = lims["as23lim"]

    
    fig, axs  = pl.subplots(3,3,sharex=True)

    mew       =  measure[:,ew].mean()
    prop_plot(axs[0,0],measure[:,con],measure[:,ew]/mew,
        {"label" : "Average EW = {:.3f}".format(mew),
         "title" : "Relative equivalent width,\n " + str(line),
         "ylabel": "Relative equivalent width",
         "xlabel": "Continuum intensity",
         "ylim"  : ewlim},regbins)

    prop_plot(axs[0,1],measure[:,con],measure[:,vel],
        {"title" : "Line centre,\n " + str(line),
         "ylabel": "Line centre [km/s]",
         "xlabel": "Continuum intensity",
         "ylim"  : vellim},regbins)

    prop_plot(axs[0,2],measure[:,con],measure[:,bot],
        {"title" : "Relative line bottom,\n " + str(line),
         "ylabel": "Relative Line min intesity",
         "xlabel": "Continuum intensity",
         "ylim"  : rellim},regbins)

    prop_plot(axs[1,0],measure[:,con],measure[:,fw13],
        {"title" : "Full width 1/3 maximum,\n " + str(line),
         "ylabel": "Relative line width",
         "xlabel": "Continuum intensity",
         "ylim"  : fw13lim},regbins)

    prop_plot(axs[1,1],measure[:,con],measure[:,fwhm],
        {"title" : "Full width half maximum,\n " + str(line),
         "ylabel": "Relative line width",
         "xlabel": "Continuum intensity",
         "ylim"  : fwhmlim},regbins)

    prop_plot(axs[1,2],measure[:,con],measure[:,fw23],
        {"title" : "Full width 2/3 maximum,\n " + str(line),
         "ylabel": "Relative line width",
         "xlabel": "Continuum intensity",
         "ylim"  : fw23lim},regbins)

    prop_plot(axs[2,0],measure[:,con],measure[:,as13],
        {"title" : "Line asymmetry 1/3 maximum,\n " + str(line),
         "ylabel": "Line asymmetry",
         "xlabel": "Continuum intensity",
         "ylim"  : as13lim},regbins)
        
    prop_plot(axs[2,1],measure[:,con],measure[:,as12],
        {"title" : "Line asymmetry half maximum,\n " + str(line),
         "ylabel": "Line width [Å]",
         "xlabel": "Continuum intensity",
         "ylim"  : as12lim},regbins)

    prop_plot(axs[2,2],measure[:,con],measure[:,as23],
        {"title" : "Line asymmetry 2/3 maximum,\n " + str(line),
         "ylabel": "Line asymmetry",
         "xlabel": "Continuum intensity",
         "ylim"  : as23lim},regbins)

    if mesbin is not None:
        axs[0,0].plot(mesbin[:,con],mesbin[:,ew]/mew,'ko')
        axs[0,1].plot(mesbin[:,con],mesbin[:,vel],'ko')
        axs[0,2].plot(mesbin[:,con],mesbin[:,bot],'ko')
        axs[1,0].plot(mesbin[:,con],mesbin[:,fw13]/line.width,'ko')
        axs[1,1].plot(mesbin[:,con],mesbin[:,fwhm]/line.width,'ko')
        axs[1,2].plot(mesbin[:,con],mesbin[:,fw23]/line.width,'ko')
        axs[2,0].plot(mesbin[:,con],(mesbin[:,as13]-measure[:,as13].mean())/line.width,'ko')
        axs[2,1].plot(mesbin[:,con],(mesbin[:,as12]-measure[:,as12].mean())/line.width,'ko')
        axs[2,2].plot(mesbin[:,con],(mesbin[:,as23]-measure[:,as23].mean())/line.width,'ko')
    

    return fig

def add_errs_linemap(fig,errs,mew,xs=None,ys=None):
    porder = np.array([9,1,0,4,2,6,5,3,7])
    if xs is None:
        xs = np.ones(11)*1.17
    if ys is None:    
        ys = [1.2,4.5,0.93,0.85,0.87,0.9,0.0015,0.0015,0.0015,1,0] 
    perr = np.abs(errs[porder,:])
    perr[0,:] = perr[0,:]/mew

    for i in range(0,9):
        s2err = np.array([perr[i,3],perr[i,0]]).reshape(2,1)
        s1err = np.array([perr[i,2],perr[i,1]]).reshape(2,1)
        yl,yu = fig.axes[i].get_ylim(); 
        dy = (yu-yl)*0.05; y = yu - s2err[1] - dy
        fig.axes[i].errorbar( xs[i],y,yerr=s1err,fmt='b' )
        fig.axes[i].errorbar( xs[i],y,yerr=s2err,fmt='r' )
    return fig

def prop_plot(ax,x,y,conf,bins=73):
    regx,regy = st.kern_reg(x,y,bins=bins)
    ax.plot(x,y,'bo',alpha=0.1,label=conf.get("label",""))
    ax.plot(regx,regy,'w',linewidth=2.1)
    ax.plot(regx,regy,'r',linewidth=1.2)
    ax.set_title(conf["title"])
    ax.set_ylabel(conf["ylabel"])
    ax.set_xlabel(conf["xlabel"])
    if "label" in conf:
        ax.legend()
    if "ylim" in conf:
        ax.set_ylim(conf["ylim"])

def addreg(fig,measure,line,bins=73,colour="r",label=""):
    bot  = 0; vel  = 1; fwhm = 2; as12 = 3; fw13 = 4; as13 = 5; fw23 = 6; as23 = 7; err  = 8; ew   = 9; con  = 10;

    x = measure[:,con]

    ys =[
            measure[:,ew]/measure[:,ew].mean(),
            measure[:,vel],
            measure[:,bot],
            measure[:,fw13]/line.width,
            measure[:,fwhm]/line.width,
            measure[:,fw23]/line.width,
            (measure[:,as13]-measure[:,as13].mean())/line.width,
            (measure[:,as12]-measure[:,as13].mean())/line.width,
            (measure[:,as23]-measure[:,as13].mean())/line.width,
        ]

    for i,y in enumerate(ys):
        regx,regy = st.kern_reg(x,y,bins=bins)
        fig.axes[i].plot(regx,regy,color=colour)

    return fig
