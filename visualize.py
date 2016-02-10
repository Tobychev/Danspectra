import matplotlib.pyplot as pl

def show_fit_with_points(spec,idx,fit,line="mean"):
    if line == "mean":
        data = spec.mean
    else:
        data = spec.spec(line)
    pl.plot(spec.lmbd,spec.mean)
    pl.plot(spec.lmbd[idx],spec.mean[idx],'ro')
    pl.plot(spec.lmbd[idx], fit[1]+ fit[0]*spec.lmbd[idx])
    pl.show()

def show_fit_on_curve(fit,xs,curve):
    pl.plot(xs,fit[0]*xs+fit[1],'r')
    pl.plot(xs,curve,'b')
    pl.show()

def plot_fits_stats(data,ref,names,title="",bins=43,cols=2,yscale="linear"):
    fig = pl.figure()
    pl.suptitle(title,fontsize=14)
    rows = len(names)/cols
    if len(names)%cols > 0:
        rows+=1

    for i,key in enumerate(names):
        ax = fig.add_subplot(rows,cols,i+1)
        ax.hist(data[key],bins=bins)
        ax.axvline(ref[key],linestyle="dashed",color="r")
        ax.set_yscale(yscale)
        ax.set_title(key)
        

    pl.tight_layout()
    pl.show()
