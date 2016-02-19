import numpy as np
import scipy as sp
import danspec as dan
import interactive as intr
import collections as col

def make_pkwin_from_linegroup(lines):
    pkwin = [] 
    for itm in lines:
        pkwin.append(list(itm.win))

    return pkwin

def make_lines_from_wins(danspec_sac,wins):
    lines = []
    for item in wins:
        lines.append(line(item,danspec_sac))

    return lines

def fit_linecores(line,xs,ys):
    guess = line.ref.argmin()
    bottom = line.idx[guess-4:guess+5]
    a,b,c = np.polyfit(xs[bottom],ys.T[bottom,:],2)
    return get_linecore(a,b,c) 

def fit_linecore(line,xs,ys):
    guess = line.ref.argmin()
    bottom = line.idx[guess-4:guess+5]
    a,b,c = np.polyfit(xs[bottom],ys[bottom],2)
    return get_linecore(a,b,c)    

def get_linecore(a,b,c):
    lam_min = -b/(2*a)
    lin_bot = np.polyval((a,b,c),lam_min)
    return lam_min,lin_bot

line_core_indices = col.namedtuple("Line_core_indices",["lcen","lbot","cont"])
lc = line_core_indices(0,1,2)


class line(object):
    def __init__(self,winbounds,group):
        self.idx  = self.__trim_line_indices(winbounds,group.ref)
        self.win  = (self.idx[0], self.idx[-1])       
        self.name = str(group.lmbd[self.idx].mean())

    def __repr__(self):
        return "Line {} [{} to {}]".format(self.name,self.idx[0],self.idx[-1])
    
    def __trim_line_indices(self,winbounds,ref):
        # Trims out values above one
        idx = range(winbounds[0],winbounds[1]+1)
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
        out   = np.zeros((len(group.frames)*3,nrows))
        i = 0
        for frame in group.frames:
            guess   = group.ref[self.idx].argmin()
            bottom  = self.idx[guess-4:guess+5]   
            a,b,c   = np.polyfit(group.lmbd[bottom],frame.data.T[bottom,:],2)
            lam_min = -b/(2*a)
            lin_bot = np.polyval((a,b,c),lam_min)
            out[i,:]  = lam_min
            out[i+1:] = lin_bot
            out[i+2:] = frame.cont.val(lam_min) 
            i +=3
        return out

