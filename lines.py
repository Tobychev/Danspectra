#encoding: utf-8
import numpy as np
import scipy.signal as ss
import danframe as dan
import interactive as intr
import collections as col
import astropy.stats as ast
import numpy.polynomial.polynomial as pol
import scipy.interpolate as si

def make_lines_from_wins(frameseries,wins):
    lines = []
    for item in wins:
        lines.append(line(item,frameseries))

    return lines

def make_splines_from_wins(frameseries,wins):
    lines = []
    for item in wins:
        lines.append(spline_line(item,frameseries))

    return lines

line_core_indices = col.namedtuple("Line_core_indices",["lcen","lbot","cont","Errs","EW","EWcont"])
lc = line_core_indices(0,1,2,3,0,1)

class line(object):
    def __init__(self,winbounds,group,weak=False):
        self.idx   = self.__trim_line_indices(winbounds,group.ref)
        self.win   = (self.idx[0], self.idx[-1])
        self.cent  = self.__get_refcentre(group)
        self.name  = "{:6.3f}".format(self.cent)
        self.width = group.lmbd[self.win[0]] - group.lmbd[self.win[1]] 

    def __repr__(self):
        return "Line {} [$\Delta \lambda$ = {:.4f} Å]".format(self.name,self.width/10)

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

            out[i,:]   = lam_min
            out[i+1,:] = lin_bot
            out[i+2,:] = frame.cont.val(self.cent) 
            # Extra error term to penalize fits that gets wildly off center, with extra weight so it *hurts*
            out[i+3,:] = np.sqrt( np.mean( (frame.data[:,test]-pred)**2,axis=1) + 2*(lam_min-self.cent)**2 )

            i +=4
        return out

    def measure_EW(self,group,vel,bot):
        nrows = group.frames[0].data.shape[0]
        out   = np.zeros((len(group.frames)*5,nrows))        
        i = 0
        for frame in group.frames:
            dlam = np.diff(group.lmbd[slice(self.idx[0]-1,self.idx[-1]+1)]).reshape((-1,1))*np.ones(nrows)
            out[i,:]   = ((frame.data[:,self.idx]-1)*dlam.T).sum(axis=1)*1e3 ## MiliÅngström 
            out[i+1,:] = frame.cont.val(self.cent) 
            i += 2
        return out

    def measure(self,group):

        cv = 299792.458
        nrows = group.frames[0].data.shape[0]
        nfram = len(group.frames)
        guess = group.ref[self.idx].argmin()
        width = len(self.idx)*0.16 # Min fraction of points to be used in fit

        vel = np.zeros((nfram,nrows))
        bot = np.zeros((nfram,nrows))
        con = np.zeros((nfram,nrows))
        err = np.zeros((nfram,nrows))
        ew  = np.zeros((nfram,nrows))
        wvar = np.zeros((nfram,nrows))
        wske = np.zeros((nfram,nrows))
        wkur = np.zeros((nfram,nrows))
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
            ew[i,:]  = self._equivalent_width(frame,nrows)
            con[i,:] = frame.cont.val(self.cent)
            vel[i,:], bot[i,:], err[i,:]       = self.__linfit(frame,bottom,test)
            wvar[i,:], wske[i,:], wkur[i,:]    = self.__ew_moments(frame,nrows,vel[i,:],bot[i,:])
            mn[i,:],var[i,:],ske[i,:],kur[i,:] = self.__moments(frame)

        vel = vel.reshape(1,-1)
        vel = cv*(vel-self.cent)/self.cent

        return (vel,bot.reshape(1,-1),con.reshape(1,-1),err.reshape(1,-1),
                ew.reshape(1,-1) ,mn.reshape(1,-1) ,var.reshape(1,-1),ske.reshape(1,-1),
                kur.reshape(1,-1),wvar.reshape(1,-1),wske.reshape(1,-1),wkur.reshape(1,-1))

    def __linfit(self,frame,bottom,test):
        cv = 299792.458
        fit   = pol.polyfit(frame.group.lmbd[bottom],frame.data[:,bottom].T,2)
        a,b,c = fit[2,:],fit[1,:],fit[0,:]
        lmin  = -b/(2*a)
        bot   = pol.polyval(lmin,fit,tensor=False)
        pred  = pol.polyval(frame.group.lmbd[test],fit)
        vel   = lmin
        # Error of fit with extra error term to penalize fits 
        # that gets wildly off center, with extra weight so it *hurts*
        err   = np.sqrt( np.mean( (frame.data[:,test]-pred)**2,axis=1) 
                                         + 2*(lmin-self.cent)**2 )
        return vel,bot,err

    def __ew_moments(self,frame,nrows,vel,bot):
        lmbd = frame.group.lmbd[self.idx].reshape(-1,1)
        dlam = np.diff(frame.group.lmbd[slice(self.idx[0]-1,self.idx[-1]+1)]).reshape((-1,1))*np.ones(nrows)
        ew   = ((frame.data[:,self.idx]-1)*dlam.T).sum(axis=1) 

        # Variance by ratio between center and outer mass
        lsel = (lmbd > vel-ew/2) & (lmbd < vel+ew/2); 
        In   = ((frame.data[:,self.idx]-1)*dlam.T*lsel.T).sum(axis=1)
        var  =  In/ew

        # Skewness by ratio between left and right mass
        lsel = (lmbd < vel)
        lft  = ((frame.data[:,self.idx]-1)*dlam.T*lsel.T).sum(axis=1); rght = ((frame.data[:,self.idx]-1)*dlam.T*np.logical_not(lsel.T)).sum(axis=1);
        cut, = np.where(rght == 0); lft[cut] = 0; rght[cut] = 1
        ske  = lft/rght-1

        # Kurtosis 
        lsel = frame.data[:,self.idx] > (1 +   bot.reshape(-1,1))/2
        up   = ((frame.data[:,self.idx]-1)*dlam.T*lsel).sum(axis=1); dwn  = ((frame.data[:,self.idx]-1)*dlam.T*np.logical_not(lsel)).sum(axis=1);
        cut, = np.where(dwn == 0); dwn[cut] = 1; up[cut] = 0
        kur  = up/dwn

        return var, ske, kur

    def _equivalent_width(self,frame,nrows):
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

class spline_line(line):

    def measure_on_block(self,group,block,con):
        nrows = block.shape[0]
        lmbd  = group.lmbd[self.idx]
        dlam  = np.diff(group.lmbd[slice(self.idx[0]-1,self.idx[-1]+1)]).reshape((-1,1))*np.ones(nrows)
        ew    = ((block-1)*dlam.T).sum(axis=1)*1e3
        splmes = np.zeros((nrows,11))
        
        splmes[:,10] = con.reshape(-1)
        splmes[:, 9] = ew.reshape(-1)
        for i,row in enumerate(block):
            mf           = self.makespline(row,lmbd,9)
            splmes[i,:9] = self.measure_spline(mf,group) 

        return splmes

    def measure(self,group):
        nrows = group.frames[0].data.shape[0]
        lmbd  = group.lmbd[self.idx]
        nfram = len(group.frames)
        block = group.frames[0].data[:,self.idx]
        con   = group.frames[0].cont.val(self.cent)
        ew    = self._equivalent_width(group.frames[0],nrows)

        for i,frm in enumerate(group.frames[1:]):
            block = np.vstack((block,frm.data[:,self.idx]) )
            ew    = np.vstack((ew, self._equivalent_width(frm,nrows)) )
            con   = np.vstack((con,frm.cont.val(self.cent)) )
        splmes = np.zeros(block[:,1:12].shape)
        splmes[:,10] = con.reshape(-1)
        splmes[:, 9] = ew.reshape(-1)

        for i,row in enumerate(block):
            mf           = self.makespline(row,lmbd,9)
            splmes[i,:9] = self.measure_spline(mf,group) 

        return splmes

    def makespline(self,spec,lmbd,kns=6):
        _,kno = np.histogram(lmbd,kns+2)
        kno   = kno[1:-2]
        return si.LSQUnivariateSpline(lmbd[::-1],spec[::-1],kno)

    def measure_spline(self,spl,group,dl=2e-5):
        lmbd = group.lmbd[self.idx]
        lmbd = np.linspace(lmbd[0],lmbd[-1],int( (lmbd[0]-lmbd[-1])/dl )) 
        bot  = spl(lmbd).min()
        icnt = spl(lmbd).argmin()
        cnt  = lmbd[icnt]
        bo12 = (1 +   bot)/2
        bo23 = (1 + 2*bot)/3
        bo13 = (2 +   bot)/3
        fwhm,as12 = self.__width_assym(spl,lmbd,bo12,cnt)
        fw13,as13 = self.__width_assym(spl,lmbd,bo13,cnt)
        fw23,as23 = self.__width_assym(spl,lmbd,bo23,cnt)
        cnt = 299792.458*(cnt-self.cent)/self.cent

        return bot,cnt,fwhm,as12,fw13,as13,fw23,as23,spl.get_residual()

    def __width_assym(self,spl,lmbd,lev,cnt):
        ilev, = np.where(spl(lmbd) <= lev)
        
        # Check that we only got one interval
        spli, = np.where(np.diff(ilev) > 1)    # Either a number or empty
        if spli.sum() > 0:
            if   len(spli) == 1 :
                ilev  = ilev[slice(spli+1)]

        wdth  = lmbd[ilev[0]] - lmbd[ilev[-1]]
        assm  = cnt  - (lmbd[ilev[0]] + lmbd[ilev[-1]])/2
        return wdth,assm
