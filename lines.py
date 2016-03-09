#encoding: utf-8
import numpy as np
import scipy.signal as ss
import danspec as dan
import interactive as intr
import collections as col
import numpy.polynomial.polynomial as pol

cv = 299792.458

def make_pkwin_from_linegroup(lines):
    pkwin = [] 
    for itm in lines:
        pkwin.append(list(itm.win))

    return pkwin

def make_lines_from_wins(frameseries,wins):
    lines = []
    for item in wins:
        lines.append(line(item,frameseries))

    return lines

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

def smooth(data,method):
    if   method == "boxcar":
        return np.convolve(data,np.array([1,1,1,1])/4.,"valid")
    elif method == "savgol":
        return ss.savgol_filter(data,5,3)
    return data


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
        nrows = group.frames[0].data.shape[0]
        out   = np.zeros((len(group.frames)*4,nrows))
        guess   = group.ref[self.idx].argmin()
        width = len(self.idx)*0.16 # Min fraction of points to be used in fit
        if width%2 == 0:
            width +=1
        dwn = int(width - 1)/2; up = dwn+1
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

