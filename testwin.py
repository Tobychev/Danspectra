import stats
import numpy as np
import kontin as con
import danframe as dan
import visualize as vis
import matplotlib.pyplot as pl

def rolling_mean(data,size=10):
    return np.convolve(data, np.ones((size,))/size, mode='valid')

def gen_all_auto_wins(spec,wins={},line="mean"):
    metod = ["over 1","top 100","top 5%","top 20", "90-95 decile","ref top","segments"]
    for m in metod:
        wins[m] = con.select_bgwin_auto(spec,m,line)

    return wins
    
def gen_man_win(spec,wins={}):
    if len(spec.bgwindows) == 0:
        wins["manual"] = con.select_bgwin_manual(spec)
    else:
        wins["manual"] = con.make_idx_from_windows(spec.bgwindows)
    return wins

def fit_all_windows(spec,wins):
    fits = {}
    for key in wins.keys():
        idx = wins[key]
        fits[key] = con.fit_continuum(spec.lmbd[idx],spec.mean[idx])

    return fits

def perturb_window_fit(spec,window,itr=103):
    nfrac     = 10
    msg       = "{:.3f} %: mean {:.5e}, std {:.5e}, rel var {:.5e}".format
    fractions = (np.linspace(0.3,1.0,nfrac)*len(window)).round()
    fits      = np.zeros((itr,2))
    perbs     = np.zeros((nfrac,4))
    for i,fraction in enumerate(fractions):
        for j in range(0,itr):
            win = np.random.choice(window,fraction)
            fits[j,0],fits[j,1] = np.polyfit(spec.lmbd[win],spec.mean[win],1)
        
        c_cont = fits[:,1] + fits[:,0]*spec.lmbd.mean() -1 # Minus one only for meanspec 
        key = round(fraction/len(window),3)
        perbs[i,:] = [key, c_cont.mean(), c_cont.std() , c_cont.std()/c_cont.mean()]
        print msg(key,perbs[i,0],perbs[i,1],perbs[i,2])

    return perbs

def perturb_all_windows(spec,windows,itr=103):
    res = {}
    for key in windows.keys():
        res[key] = perturb_window_fit(spec,windows[key],itr)

    return res

def gen_frame_continuum(spec,win):
    fits = con.frame_fit(spec,win)
    return fits[1,:] + fits[0,:]*spec.lmbd.mean()

def all_wins_frame_continuum(spec,wins):
    lines = 800
    names = {}
    cont  = np.zeros((lines,
                      len(wins.keys()) ))
    for i,key in enumerate(wins.keys()):
        cont[:,i]  = gen_frame_continuum(spec,wins[key])
        names[key] = i

    return cont,names

def compare_win_continua(spec,wins,centre="mean",cols=2,bins=39,plot=True):
    conts,winame = all_wins_frame_continuum(spec,wins)

    if centre == "mean":
        conts = (conts.T - conts.mean(axis=1)).T #transpose to allow broadcasting
        title = "deviation from ensamble mean"
    elif centre == "median":
        centre = np.median(conts,axis=1)
        conts  = (conts.T - centre).T
        title = "deviation from ensamble median"
    elif centre == "refmax":
        centre = spec.data[:,spec.ref.argmax()]
        conts  = (conts.T - centre).T
        title = "deviation from highest pixel"
    elif centre == "smoothmax":
        centre = spec.data[:,(spec.ref.argmax()- 5):(spec.ref.argmax()+ 5)].mean(axis=1)
        conts  = (conts.T - centre).T
        title = "deviation from smoothed highest pixel"

    if plot:
        fig = pl.figure()
        pl.suptitle(title,fontsize=14)
        rows = len(winame)/cols
        if len(winame)%cols > 0:
            rows+=1
        for i,key in enumerate(winame):
            ax = fig.add_subplot(rows,cols,i+1)
            ax.hist(conts[:,i],bins=bins)
            ax.axvline(0.0,linestyle="dashed",color="r")
            ax.set_title(key)
        pl.tight_layout()
        pl.show()

    std  = conts.std(axis=0)
    mean = conts.mean(axis=0)
    print "\n"+ title
    print "\n{:>12}   {:<10} {}".format("","mean","std")
    for i,name in enumerate(winame):
        print "{:>12}: {:+.4e} {:6.4e}".format(name,mean[i],std[i])

    return mean,std 

def smooth_test(spec,win,size):
    rows = 800
    data = con.frame_subtract_cont(spec,win)
    lmd = rolling_mean(spec.lmbd,size)
    spe = np.zeros((rows,lmd.size))
    res = np.zeros((rows,2))

    for i in range(0,rows):
        spe[i,:] = rolling_mean(data[i,:],size)
        res[i,0] = spe[i,spe[i,:] > 0].sum()
        res[i,1] = spe[i,(spe[i,:] > -0.01) & (spe[i,:] < 0)].sum()

    return lmd,spe,res

def all_smooth_test(spec,wins,smooth=10,plot=False,bins=43,cols=2):
    winame = wins.keys()
    title  = "Integral of signal over continuum"
    if plot:
        fig = pl.figure()
        pl.suptitle(title,fontsize=14)
        rows = len(winame)/cols
        if len(winame)%cols > 0:
            rows+=1

    print "\n"+ title
    print "\n{:>12}  {:<7} {:<7} {:<8} {}".format("","mean","std","max","max_row")
    for i,key in enumerate(wins.keys()):
        lm,sp,re = smooth_test(spec,wins[key],smooth)
        print "{:>12}: {:7.3f} {:7.3f} {:8.3f} {}".format(key,re[:,0].mean(),re[:,0].std(),re[:,0].max(),re[:,0].argmax())
        if plot:
            ax = fig.add_subplot(rows,cols,i+1)
            ax.hist(re[:,0],bins=bins)
            ax.axvline(0.0,linestyle="dashed",color="r")
            ax.set_title(key)
    if plot:
        pl.tight_layout()
        pl.show()

def sim_poisson_noise(spec,rows):
    fil = "{}_{}__contmap.fits"
    cont_frame = dan.f.open(spec.Dir+fil.format(spec.wave,spec.series))
    mean_cont  = cont_frame[0].data.mean()

    return np.random.poisson( spec.ref*mean_cont,(rows,len(spec.ref)) )/mean_cont

def fit_on_poisson_noise(spec,rows):
    wins  = gen_all_auto_wins(spec,{},line="ref")

    ref_fit = {}
    for key in wins.keys():
        idx = wins[key]
        ref_fit[key] = np.polyfit(spec.lmbd[idx],spec.ref[idx],1)

    fits = {}
    for key in wins.keys():
        idx  = wins[key]
        data = sim_poisson_noise(spec,rows)
        fits[key] = np.polyfit(spec.lmbd[idx],data.T[idx,:],1)

    return wins.keys(),ref_fit,fits

def test_fit_with_noise(spec,rows,plot_excess=False,plot_cont_corr=False,plot_fit_par=False,
                        bins=43,cols=2,yscale="linear",cut=0):
    names,fit_ref,fits = fit_on_poisson_noise(spec,rows)

    title = "distribution of intercept"
    ms  = {key: fits[key][1,:] for key in names}
    ref = {key: fit_ref[key][1] for key in names}
    stats.print_dict_stats_with_ref(ms,ref,names,title)
    if plot_fit_par:
        vis.plot_fits_stats(ms,ref,names,title,bins,cols)

        
    title = "distribution of slope"
    ks  = {key: fits[key][0,:] for key in names}
    ref = {key: fit_ref[key][0] for key in names}
    stats.print_dict_stats_with_ref(ms,ref,names,title)
    if plot_fit_par:
        vis.plot_fits_stats(ks,ref,names,title,bins,cols)

    excess   = {}
    ref      = {}
    rmse_top = {}
    ref_rmse = {}
    for key in names:
        kx = np.outer(ks[key],spec.lmbd)
        m  = ms[key].reshape((rows,1))
        y  = kx + m
        extrema = np.unique(np.array([np.abs(y[:,0]).argmax(),np.abs(y[:,0]).argmin(),
                                      np.abs(y[:,-1]).argmax(),np.abs(y[:,-1]).argmin()])  )
        concorr = spec.ref - y
        ref_cor = spec.ref - (fit_ref[key][0]*spec.lmbd + fit_ref[key][1])
        zero    = np.zeros(concorr.shape)    

        excess[key]   = np.where(concorr > 0, concorr,zero).sum(axis=1) 
        rmse_top[key] = stats.rmse_at_zero_of_topq(concorr,q=50,ax=1)
        ref[key]      = ref_cor[ref_cor > 0].sum() 
        ref_rmse[key] = stats.rmse_at_zero_of_topq(ref_cor,q=50)
        
        if plot_cont_corr:
            fig = pl.figure()
            ax  = fig.add_subplot(111)
            for e in extrema:
                ax.plot(spec.lmbd,concorr[e,:])
            
            pl.show()
        
    stats.print_dict_stats_with_ref(excess,ref,names,"Excess")
    stats.print_dict_stats_with_ref(rmse_top,ref_rmse,names,"RMSE of top half")
    if plot_excess:
        title = "excess"
        vis.plot_fits_stats(excess,ref,names,title,bins,cols,yscale,cutoff=cut)
        
    return [excess,ms,ks,fit_ref,rmse_top,ref_rmse]
