#encoding: utf-8
import numpy as np
import scipy.signal as ss
import danframe as dan
import interactive as intr
import collections as col
import numpy.polynomial.polynomial as pol

def make_pkwin_from_linegroup(lines):
    pkwin = [] 
    for itm in lines:
        pkwin.append(list(itm.win))

    return pkwin

def fit_linecores(line,xs,ys):
    guess = line.ref.argmin()
    bottom = line.idx[guess-4:guess+5]
    a,b,c = np.polyfit(xs[bottom],ys.T[bottom,:],2)
    return get_linecore(a,b,c) 

def get_linecore(a,b,c):
    lam_min = -b/(2*a)
    lin_bot = np.polyval((a,b,c),lam_min)
    return lam_min,lin_bot

line_core_indices = col.namedtuple("Line_core_indices",["lcen","lbot","cont","Errs","EW","EWcont"])
lc = line_core_indices(0,1,2,3,0,1)

class line(object):
    def __init__(self,winbounds,group,weak=False):
        self.idx  = self.__trim_line_indices(winbounds,group.ref)
        self.win  = (self.idx[0], self.idx[-1])
        self.cent = self.__get_refcentre(group)
        self.name = "{:6.3f}".format(self.cent)
        self.weak = weak

    def __repr__(self):
        return "Line {} [{} to {}]".format(self.name,self.idx[0],self.idx[-1])

    def __get_refcentre(self,group):
        bot   = slice(group.ref[slice(self.win[0],self.win[1])].argmin()-3,group.ref[slice(self.win[0],self.win[1])].argmin()+4)
        a,b,c = np.polyfit(group.lmbd[self.idx[bot]],group.ref[self.idx[bot]],2)
        return -b/(2*a)

    def __trim_line_indices(self,winbounds,ref):
        # Trims out values above one
        idx = np.arange(winbounds[0],winbounds[1]+1)
        tmp = np.where(ref[idx] > 1)[0]
        cut = np.where(np.diff(tmp) > 1)[0] 
        if len(tmp) <  2:
            return idx
        elif len(tmp) == 2 and len(cut) == 1:
            return idx
        elif len(cut) > 1:
            return idx[tmp[cut]:tmp[cut+2]]
        elif tmp[0] == 0:
            return idx[tmp[-1]:]
        else:
            return idx[:tmp[1]]

    def __find_minima(self,curve):
        deriv = np.diff(curve)
        mid   = np.where(np.diff(np.sign(deriv)))[0]+1
        infl  = np.concatenate( ([0],mid,[len(deriv)]) )
        mins  = []
        for i,point in enumerate(infl):
            if self.__is_minima(deriv,point):                      
                if i!= 0 and i < len(infl)-1:
                    mins.append( (infl[i-1],point, infl[i+1]) )

        return mins
        
    def __is_minima(self,curve,point):
        if curve[point-3:point].mean() < 0 and curve[point:point+3].mean() > 0:
            return True
        else:
            return False

    def __select_bottom(self,spectra,cent):
        shift = np.arange(-3,4)
        bot = self.idx[cent+shift]
        for i in range(1,20):
            mn = spectra[bot].argmin()
            cent = cent + shift[mn]
            bot = self.idx[cent+shift]     
            if not shift[mn]:
                return bot

        raise Exception("Did not converge")        

    def subdivide_line(self, group):
        mins = self.__find_minima(group.ref[self.idx])
        if len(mins) > 1:
            lines = []    
            mins  = intr.manual_delete_minima(mins,group)
            for li in mins:
                lines.append(line( (line.idx[ li[0]] ,line.idx[li[2]]),group ))

            return lines
        else:
            return line
    
    def measure_linecores(self,group):
        cv = 299792.458
        nrows = group.frames[0].data.shape[0]
        out   = np.zeros((len(group.frames)*4,nrows))
        guess   = group.ref[self.idx].argmin()
        width = len(self.idx)*0.16 # Min fraction of points to be used in fit
        if width%2 == 0:
            width +=1
        dwn = int((width - 1)/2); up = dwn+1
        bottom  = self.idx[guess + np.arange(-dwn,up)]   
        test    = self.idx[guess + np.arange(-(dwn+1),(up+1))]
        i = 0
        for frame in group.frames:
            fit     = pol.polyfit(group.lmbd[bottom],frame.data[:,bottom].T,2)
            a,b,c   = fit[2,:],fit[1,:],fit[0,:]
            lam_min = -b/(2*a)
            lin_bot = pol.polyval(lam_min,fit,tensor=False)
            pred    = pol.polyval(group.lmbd[test],fit)

            out[i,:]   = cv*(lam_min-self.cent)/self.cent
            out[i+1,:] = lin_bot
            out[i+2,:] = frame.cont.val(self.cent) 
            # Extra error term to penalize fits that gets wildly off center, with extra weight so it *hurts*
            out[i+3,:] = np.sqrt( np.mean( (frame.data[:,test]-pred)**2,axis=1) + 2*(lam_min-self.cent)**2 )

            i +=4
        return out

    def measure_EW(self,group):
        nrows = group.frames[0].data.shape[0]
        out   = np.zeros((len(group.frames)*2,nrows))        
        i = 0
        for frame in group.frames:
            dlam = np.diff(group.lmbd[slice(self.idx[0]-1,self.idx[-1]+1)]).reshape((-1,1))*np.ones(nrows)
            out[i,:]   = ((frame.data[:,self.idx]-1)*dlam.T).sum(axis=1)*1e3 ## MiliÅngström
            out[i+1,:] = frame.cont.val(self.cent) 

            i +=2
        return out

    def measure(self,group):
        nrows = group.frames[0].data.shape[0]
        nfram = len(group.frames)
        guess = group.ref[self.idx].argmin()
        width = len(self.idx)*0.16 # Min fraction of points to be used in fit

        vel = np.zeros((nfram,nrows))
        bot = np.zeros((nfram,nrows))
        con = np.zeros((nfram,nrows))
        err = np.zeros((nfram,nrows))
        ew  = np.zeros((nfram,nrows))
        mn  = np.zeros((nfram,nrows))
        var = np.zeros((nfram,nrows))
        ske = np.zeros((nfram,nrows))
        kur = np.zeros((nfram,nrows))

        if width%2 == 0:
            width +=1
        dwn = int((width - 1)/2); up = dwn+1

        bottom  = self.idx[guess + np.arange(-dwn,up)]   
        test    = self.idx[guess + np.arange(-(dwn+1),(up+1))]
        for i,frame in enumerate(group.frames):
            con[i,:] = frame.cont.val(self.cent)
            vel[i,:],bot[i,:],err[i,:] = self.__linfit(frame,bottom,test)
            ew[i,:] = self.__equivalent_width(frame,nrows)
            mn[i,:],var[i,:],ske[i,:],kur[i,:] = self.__moments(frame)

        return (vel.reshape(1,-1),bot.reshape(1,-1),con.reshape(1,-1),err.reshape(1,-1),
                ew.reshape(1,-1) ,mn.reshape(1,-1) ,var.reshape(1,-1),ske.reshape(1,-1),
                kur.reshape(1,-1))

    def __linfit(self,frame,bottom,test):
        cv = 299792.458
        fit   = pol.polyfit(frame.group.lmbd[bottom],frame.data[:,bottom].T,2)
        a,b,c = fit[2,:],fit[1,:],fit[0,:]
        lmin  = -b/(2*a)
        bot   = pol.polyval(lmin,fit,tensor=False)
        pred  = pol.polyval(frame.group.lmbd[test],fit)
        vel   = cv*(lmin-self.cent)/self.cent
        # Error of fit with extra error term to penalize fits 
        # that gets wildly off center, with extra weight so it *hurts*
        err   = np.sqrt( np.mean( (frame.data[:,test]-pred)**2,axis=1) 
                                         + 2*(lmin-self.cent)**2 )
        return vel,bot,err

    def __equivalent_width(self,frame,nrows):
        dlam = np.diff(frame.group.lmbd[slice(self.idx[0]-1,self.idx[-1]+1)]
                       ).reshape((-1,1))*np.ones(nrows)
        return ((frame.data[:,self.idx]-1)*dlam.T).sum(axis=1)*1e3 ## MiliÅngström          

    def __moments(self,frame):
        x    = frame.group.lmbd[self.idx]
        dpdf = (1-frame.data[:,self.idx]/frame.data[:,self.idx].max(axis=1).reshape(-1,1))
        dpdf = dpdf/dpdf.sum(axis=1).reshape(-1,1)
        mu   = np.sum(dpdf*x,axis=1).reshape(-1,1) # Reshaping enables broadcasting
        mu2  = np.sum(dpdf*(x-mu)**2,axis=1)
        mu3  = np.sum(dpdf*(x-mu)**3,axis=1)
        mu4  = np.sum(dpdf*(x-mu)**4,axis=1)
        mu   = mu.reshape(-1) # Undoing reshape to allow assignment


        skew = mu3/mu2**(3/2) 
        kurt = (mu4/mu2**2 - 3)
        return mu,mu2,skew,kurt


class binned_framegroup(object):

    def __init__(self,line,group,bin_quant,cuts=None):
        self.idx   = line.idx
        self.cent  = line.cent
        self.group = group
        self.data,self.cont = self.__bin_by_quant(bin_quant,cuts)
    
    def __bin_by_quant(self,quant,cuts):
        try:             
            cuts, = np.where(cuts.reshape(-1))
        except AttributeError:
            cuts = np.ones(len(quant))

        datablock = self.group.frames[0].data[:,self.idx]        
        contblock = self.group.frames[0].cont.val(self.cent)
        for frm in self.group.frames[1:]:
            datablock = np.concatenate((datablock,frm.data[:,self.idx]),axis=0)
            contblock = np.concatenate((contblock,frm.cont.val(self.cent)),axis=0)
        datablock = datablock[cuts,:]
        contblock = contblock[cuts]

        counts, bins = ast.histogram(quant[cuts],bins='blocks')
        sorting = np.digitize(quant[cuts],bins,right=True)
        binned  =  datablock[sorting == 0,:]
        con     =  contblock[sorting == 0]
        for i in np.unique(sorting)[1:]:
            binned = np.vstack( (binned,np.mean(datablock[sorting == i,:],axis=0) ))
            con    = np.vstack( (con,np.mean(contblock[sorting == i]) ))
        return binned,con

    def measure(self):
        nfram = len(self.data)
        guess = self.group.ref[self.idx].argmin()
        width = len(self.idx)*0.16 # Min fraction of points to be used in fit

        vel = np.zeros(nfram)
        bot = np.zeros(nfram)
        err = np.zeros(nfram)
        ew  = np.zeros(nfram)
        mn  = np.zeros(nfram)
        var = np.zeros(nfram)
        ske = np.zeros(nfram)
        kur = np.zeros(nfram)

        if width%2 == 0:
            width +=1
        dwn = int((width - 1)/2); up = dwn+1

        bottom  = guess + np.arange(-dwn,up) 
        test    = guess + np.arange(-(dwn+1),(up+1))

        vel,bot,err = self.__linfit(bottom,test)
        ew = self.__equivalent_width(nfram)
        mn,var,ske,kur = self.__moments()

        return (vel,bot,err,ew ,mn ,var,ske,kur)

    def __linfit(self,bottom,test):
        cv = 299792.458
        fit   = pol.polyfit(self.group.lmbd[self.idx[bottom]],self.data[:,bottom].T,2)
        a,b,c = fit[2,:],fit[1,:],fit[0,:]
        lmin  = -b/(2*a)
        bot   = pol.polyval(lmin,fit,tensor=False)
        pred  = pol.polyval(self.group.lmbd[test],fit)
        vel   = cv*(lmin-self.cent)/self.cent
        # Error of fit with extra error term to penalize fits 
        # that gets wildly off center, with extra weight so it *hurts*
        err   = np.sqrt( np.mean( (self.data[:,test]-pred)**2,axis=1) 
                                         + 2*(lmin-self.cent)**2 )
        return vel,bot,err

    def __equivalent_width(self,nrows):
        dlam = np.diff(self.group.lmbd[slice(self.idx[0]-1,self.idx[-1]+1)]
                       ).reshape((-1,1))*np.ones(nrows)
        return ((self.data-1)*dlam.T).sum(axis=1)*1e3 ## MiliÅngström          

    def __momentss(self):
        x    = self.group.lmbd[self.idx]
        dpdf = (1-self.data/self.data.max(axis=1).reshape(-1,1))
        dpdf = dpdf/dpdf.sum(axis=1).reshape(-1,1)
        mu   = np.sum(dpdf*x,axis=1).reshape(-1,1) # Reshaping enables broadcasting
        mu2  = np.sum(dpdf*(x-mu)**2,axis=1)
        mu3  = np.sum(dpdf*(x-mu)**3,axis=1)
        mu4  = np.sum(dpdf*(x-mu)**4,axis=1)
        mu   = mu.reshape(-1) # Undoing reshape to allow assignment


        skew = mu3/mu2**(3/2) 
        kurt = (mu4/mu2**2 - 3)
        return mu,mu2,skew,kurt

    
    def __moments(self):
        x    = self.group.lmbd[self.idx]
        rows = 16
        mu, mu2, mu3, mu4 = np.zeros(rows),np.zeros(rows),np.zeros(rows),np.zeros(rows)
        for i in range(0,rows):
            dpdf = (1-self.data[i,:]/self.data[i,:].max())
            dpdf = dpdf/dpdf.sum()
            mu[i] = np.sum(dpdf*x) # Reshaping enables broadcasting
            mu2[i]= np.sum(dpdf*(x-mu[i])**2)
            mu3[i]= np.sum(dpdf*(x-mu[i])**3)
            mu4[i]= np.sum(dpdf*(x-mu[i])**4)

        skew = mu3/mu2**(3/2) 
        kurt = (mu4/mu2**2 - 3)
        return mu,mu2,skew,kurt
