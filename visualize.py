import numpy as np
import matplotlib.pyplot as pl
import matplotlib.cm as cm
import scipy.stats as sta
import lines as lin
import stats as st

def get_colours(num,colour_name="Oranges"):
    return iter(cm.__dict__[colour_name](np.linspace(0,1,num)))

def show_linegroup_on_curve(lines,xs,ys):
    pl.step(xs,ys)
    for line in lines:
        pl.plot(line.lmbd,line.ref,'ro-')

    pl.show()

def show_autoline_on_curve(autoline,xs,ys):
    line = list(autoline)
    pl.step(xs,ys)
    pl.plot(xs[line],ys[line],'ro-')
    pl.show()

def make_line_on_curve(line,xs,ys,ax,colour=""):
    ax.plot(xs[line.idx],ys[line.idx],'o',c=colour)
    return ax

def show_spec_lines(spec,lines):
    fig = pl.figure()
    ax  = fig.add_subplot(1,1,1)
    colour = get_colours(len(lines))

    ax.step(spec.lmbd,spec.ref)
    for i,line in enumerate(lines):
        ax = make_line_on_curve(line,spec.lmbd,spec.ref,ax,colour=next(colour))

    pl.show()

def show_dict_wins_on_curve(wins,xs,ys):
    colour = get_colours(len(wins))
    pl.step(xs,ys)
    for i,name in enumerate(wins):
        idx = wins[name]
        pl.plot(xs[idx],ys[idx],'o',c=next(colour))

    pl.show()

def show_fit_with_points(spec,idx,fit,line="mean"):
    if line == "mean":
        data = spec.mean
    else:
        data = spec.spec(line)
    pl.step(spec.lmbd,spec.mean)
    pl.plot(spec.lmbd[idx],spec.mean[idx],'ro')
    pl.plot(spec.lmbd[idx], fit[1]+ fit[0]*spec.lmbd[idx])
    pl.show()

def show_fit_on_curve(fit,xs,curve):
    pl.plot(xs,fit[0]*xs+fit[1],'r')
    pl.step(xs,curve,'b')
    pl.show()

def plot_fits_stats(data,ref,names,title="",bins=43,cols=2,yscale="linear",cutoff=0):
    fig = pl.figure()
    pl.suptitle(title,fontsize=14)
    rows = len(names)/cols
    if len(names)%cols > 0:
        rows+=1

    for i,key in enumerate(names):
        if cutoff > 0:
            idx = np.where(data[key] > np.percentile(data[key],q=cutoff))
        else:
            idx = list(range(0,len(data[key])))
        ax = fig.add_subplot(rows,cols,i+1)
        ax.hist(data[key][idx],bins=bins)
        ax.axvline(ref[key],linestyle="dashed",color="r")
        ax.set_yscale(yscale)
        ax.set_title(key)
        
    pl.tight_layout()
    pl.show()

def show_contfit(frame,line):
    xs  = frame.group.lmbd
    ys  = frame.spec(line)
    fit = (frame.cont.fit["k"][line],frame.cont.fit["m"][line])
    show_fit_on_curve(fit,xs,ys)

def show_line_and_corefit(line,frame,row,width=3,fast=True):
    spe     = frame.data[row,:]
    ref,lmd = frame.group.ref,frame.group.lmbd
    if fast:
        guess   = frame.group.ref[line.idx].argmin()
        bottom  = line.idx[slice(guess-4,guess+5)]
    else:
        guess   = spe[line.idx].argmin()
        bottom  = line.idx[guess+np.arange(-width,width+1)]
    a,b,c   = np.polyfit(lmd[bottom],spe[bottom],2)
    fit     = np.polyval((a,b,c),lmd[line.idx]); 
    idfit = line.idx[fit < 1.05]
    fit   = fit[fit < 1.05]

    pl.step(lmd[line.idx],spe[line.idx],'b')
    pl.plot(lmd[bottom],spe[bottom],'*k')
    pl.plot(lmd[idfit],fit,'r')

    pl.show()

def plot_linemap(measure,line,binned=()):
    vel = 0; bot = 1; con = 2; err = 3; ew  = 4; mn  = 5; var = 6; ske = 7; kur = 8
    cuts = measure[err] < np.percentile(measure[err],89)

    pl.subplot(3,2,1)
    regx,regy = st.kern_reg(measure[con][cuts],measure[ew][cuts],bins=73)
    pl.plot(measure[con][cuts],measure[ew][cuts],'bo',alpha=0.2)
    pl.plot(regx,regy,'r')
    pl.title("Equivalent width, " + str(line))
    pl.ylabel("Equivalent width")
    pl.xlabel("Continuum intensity")

    pl.subplot(3,2,2)
    regx,regy = st.kern_reg(measure[con][cuts],measure[vel][cuts],bins=73)
    pl.plot(measure[con][cuts],measure[vel][cuts],'bo',alpha=0.2)
    pl.plot(regx,regy,'r')
    pl.title("Line centre, " + str(line))
    pl.ylabel("Line centre")
    pl.xlabel("Continuum intensity")

    pl.subplot(3,2,3)
    regx,regy = st.kern_reg(measure[con][cuts],measure[mn][cuts],bins=73)
    pl.plot(measure[con][cuts],measure[mn][cuts],'bo',alpha=0.2)
    pl.plot(regx,regy,'r')
    pl.title("Line mean, " + str(line))
    pl.ylabel("Skewness")
    pl.xlabel("Continuum intensity")

    pl.subplot(3,2,4)
    regx,regy = st.kern_reg(measure[con][cuts],measure[bot][cuts]/measure[con][cuts],bins=73)
    pl.plot(measure[con][cuts],measure[bot][cuts]/measure[con][cuts],'bo',alpha=0.2)
    pl.plot(regx,regy,'r')
    pl.title("Relative line bottom, " + str(line))
    pl.ylabel("Relative Line min intesity")
    pl.xlabel("Continuum intensity")

    pl.subplot(3,2,5)
    regx,regy = st.kern_reg(measure[con][cuts],measure[kur][cuts],bins=73)
    pl.plot(measure[con][cuts],measure[kur][cuts],'bo',alpha=0.2)
    pl.plot(regx,regy,'r')
    pl.title("Line kurtosis, " + str(line))
    pl.ylabel("Kurtosis")
    pl.xlabel("Continuum intensity")

    pl.subplot(3,2,6)
    regx,regy = st.kern_reg(measure[con][cuts],measure[var][cuts],bins=73)
    pl.plot(measure[con][cuts],measure[var][cuts],'bo',alpha=0.2)
    pl.plot(regx,regy,'r')
    pl.title("Line variance, " + str(line))
    pl.ylabel("variance")
    pl.xlabel("Continuum intensity")
    pl.subplots_adjust(left=0.1, bottom=0.07, right=0.95, top=0.95,wspace=0.43, hspace=0.43)


    if len(binned):
        (vel,bot,err,ew ,mn ,var,ske,kur) = np.arange(0,8)
        mesbinn,cont = binned
        pl.subplot(3,2,1)
        pl.plot(cont,mesbinn[ew]/men,'ko')

        pl.subplot(3,2,2)
        pl.plot(cont,mesbinn[vel],'ko')

        pl.subplot(3,2,3)
        pl.plot(cont,mesbinn[mn],'ko')
   
        pl.subplot(3,2,4)
        pl.plot(cont,mesbinn[bot]/cont,'ko')

        pl.subplot(3,2,5)
        pl.plot(cont,mesbinn[kur],'ko')

        pl.subplot(3,2,6)
        pl.plot(cont,mesbinn[var],'ko')
    

    pl.show()

def kde(measure):
    rt = sta.gaussian_kde(measure)
    x  = np.linspace(measure.min(),measure.max(),121)
    pl.plot(x,rt(x))
    pl.show()

def spline_linemap(measure,line,mesbin=None,lims=None,errs=None):
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
        as12lim = lims["as12lim"]

    
    fig, axs  = pl.subplots(3,3,sharex=True)
    fig.subplots_adjust(wspace=0.3,hspace=0.3)

    mew       =  measure[:,ew].mean()
    prop_plot(axs[0,0],measure[:,con],measure[:,ew]/mew,
        {"label" : "Average EW = {:.3f}".format(mew),
         "title" : "Relative equivalent width,\n " + str(line),
         "ylabel": "Relative equivalent width",
         "xlabel": "Continuum intensity",
         "ylim"  : ewlim})

    prop_plot(axs[0,1],measure[:,con],measure[:,vel],
        {"title" : "Line centre,\n " + str(line),
         "ylabel": "Line centre [km/s]",
         "xlabel": "Continuum intensity",
         "ylim"  : vellim})

    prop_plot(axs[0,2],measure[:,con],measure[:,bot]/measure[:,con],
        {"title" : "Relative line bottom,\n " + str(line),
         "ylabel": "Relative Line min intesity",
         "xlabel": "Continuum intensity",
         "ylim"  : rellim})

    prop_plot(axs[1,0],measure[:,con],measure[:,fw13]/line.width,
        {"title" : "Full width 1/3 maximum,\n " + str(line),
         "ylabel": "Relative line width",
         "xlabel": "Continuum intensity",
         "ylim"  : fw13lim})

    prop_plot(axs[1,1],measure[:,con],measure[:,fwhm]/line.width,
        {"title" : "Full width half maximum,\n " + str(line),
         "ylabel": "Relative line width",
         "xlabel": "Continuum intensity",
         "ylim"  : fwhmlim})

    prop_plot(axs[1,2],measure[:,con],measure[:,fw23]/line.width,
        {"title" : "Full width 2/3 maximum,\n " + str(line),
         "ylabel": "Relative line width",
         "xlabel": "Continuum intensity",
         "ylim"  : fw23lim})

    prop_plot(axs[2,0],measure[:,con],(measure[:,as13]-measure[:,as13].mean())/line.width,
        {"title" : "Line asymmetry 1/3 maximum,\n " + str(line),
         "ylabel": "Line asymmetry",
         "xlabel": "Continuum intensity",
         "ylim"  : as13lim})
#        })
        
    prop_plot(axs[2,1],measure[:,con],(measure[:,as12]-measure[:,as12].mean())/line.width,
        {"title" : "Line asymmetry half maximum,\n " + str(line),
         "ylabel": "Line width [Å]",
         "xlabel": "Continuum intensity",
         "ylim"  : as12lim})
#         })
    prop_plot(axs[2,2],measure[:,con],(measure[:,as23]-measure[:,as23].mean())/line.width,
        {"title" : "Line asymmetry 2/3 maximum,\n " + str(line),
         "ylabel": "Line asymmetry",
         "xlabel": "Continuum intensity",
         "ylim"  : as23lim})
#         })

    if mesbin is not None:
        axs[0,0].plot(mesbin[:,con],mesbin[:,ew]/mew,'ko')
        axs[0,1].plot(mesbin[:,con],mesbin[:,vel],'ko')
        axs[0,2].plot(mesbin[:,con],mesbin[:,bot]/mesbin[:,con],'ko')
        axs[1,0].plot(mesbin[:,con],mesbin[:,fw13]/line.width,'ko')
        axs[1,1].plot(mesbin[:,con],mesbin[:,fwhm]/line.width,'ko')
        axs[1,2].plot(mesbin[:,con],mesbin[:,fw23]/line.width,'ko')
        axs[2,0].plot(mesbin[:,con],(mesbin[:,as13]-measure[:,as13].mean())/line.width,'ko')
        axs[2,1].plot(mesbin[:,con],(mesbin[:,as12]-measure[:,as12].mean())/line.width,'ko')
        axs[2,2].plot(mesbin[:,con],(mesbin[:,as23]-measure[:,as23].mean())/line.width,'ko')
    

    return fig

def dan_errplot(fig,errs,xs=None,ys=None):
    porder = np.array([9,1,0,4,2,6,5,3,7])
    if xs is None:
        xs = np.ones(11)*1.15
    if ys is None:    
        ys = [0.6,-4,1.0,0.2,0.2,0.1,-0.15,-0.15,-0.08,-0.08,0] 
    perr = np.abs(errs[porder,:])

    for i in range(0,9):
        s2err = np.array([perr[i,3],perr[i,0]]).reshape(2,1)
        s1err = np.array([perr[i,2],perr[i,1]]).reshape(2,1)
        fig.axes[i].errorbar( xs[i],ys[i],yerr=s1err,fmt='b' )
        fig.axes[i].errorbar( xs[i],ys[i],yerr=s2err,fmt='r' )
    return fig

def prop_plot(ax,x,y,conf):
    regx,regy = st.kern_reg(x,y,bins=73)
    ax.plot(x,y,'bo',alpha=0.1,label=conf.get("label",""))
    ax.plot(regx,regy,'r')
    ax.set_title(conf["title"])
    ax.set_ylabel(conf["ylabel"])
    ax.set_xlabel(conf["xlabel"])
    if "label" in conf:
        ax.legend()
    if "ylim" in conf:
        ax.set_ylim(conf["ylim"])
