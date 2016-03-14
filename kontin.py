import matplotlib.pyplot as pl
import numpy as np

def frame_fit(frame,idx):
    return np.polyfit(frame.lmbd[idx],frame.data.transpose()[idx,:],1)

def frame_subtract_cont(frame,idx):
    fits = frame_fit(frame,idx)
    cont = fits[1,:].reshape(len(fits[1,:]),1) + np.outer(fits[0,:],frame.lmbd)
    return frame.data - cont
    
def row_subtract_cont(frame,line,idx):
    b,a = row_fit_continuum(frame,line,idx)
    cont = a + b*frame.lmbd
    return frame.spec(line)-cont 

def row_fit_continuum(frame,line,idx):
    fit = np.polyfit(frame.lmbd[idx],frame.frame(line)[idx],1)
    return fit

def make_idx_from_windows(windows):
    idbg = []
    for win in windows:
        idbg = idbg + list(range(win[0],win[1]+1))
    return idbg

def make_windows_from_idx(idx):
    window  = []
    windows = []
 
    window.append(idx[0])
    # The subtraction allows for a cleaner
    # condition in the loop, by preventing an
    # out-of-bounds attempt. The add in the last
    # edge by hand.
    for i in range(1,len(idx)-1):
        if idx[i]+1 != idx[i+1]:
            print(idx[i]+1, idx[i+1])
            window.append(idx[i])
            windows.append(window)
            window = [idx[i+1]]
    window.append(idx[-1])
    windows.append(window)
    return windows

def select_bgwin_auto(danframe,metod,row="mean",
                        npoint=100,q=50):
    if row == "mean":
        data = danframe.mean
    elif row == "ref":
        data = danframe.ref
    else:
        data = danframe.frame(row)

    if   metod == "over 1":
        return np.flatnonzero(data > 1)
    elif metod == "top 100":
        return top_number(data,100)
    elif metod == "top 5%":
        return top_number(data,round(len(data)*0.05))
    elif metod == "top 20":
        return top_number(data,20)
    elif metod == "90-95 decile":
        t5  = top_number(data,round(len(data)*0.05))[-1]
        t10 = top_number(data,round(len(data)*0.1))
        return t10[ t10 >= t5 ]
    elif metod == "ref top":       
        return list(range( danframe.ref.argmax()- 5, danframe.ref.argmax()+ 5))
    elif metod == "segments":
        return top_of_segments(data,npoint,q)

def top_of_segments(data,npoint,q):
    ids = (data > np.percentile(data,q)).nonzero()[0]
    nregion = len(ids)/npoint
    perreg  = npoint/nregion
    regions = np.array_split(ids,nregion)
    idx = np.array( [])
    # Top perreg of data in each region
    # have global indices given by reg, 
    # top_number returns indiced local to data[reg]
    for reg in regions:
         idx = np.hstack( (idx, reg[top_number(data[reg],perreg)]) )
    return idx.astype("int")

def top_number(data,number):
    idx = np.argpartition(data,-number)[-number:]
    idx.sort()
    return idx

def save_bgwin_from_idx(danframe,idx):
    wins = make_windows_from_idx(idx)
    danframe.set_bgwindows(wins)


class continua(object):
    def __init__(self,danframe,method):
        self.frame = danframe
        self.fit   = self.auto_fit_frame(method,danframe.group.ref)

    def norm(self):
        dim = len(self.fit["m"])
        return np.outer(self.fit["k"],self.frame.group.lmbd)+ self.fit["m"].reshape((dim,1))

    def val(self,lam):
        return self.fit["k"]*lam+ self.fit["m"]

    def __top_number(self,data,number):
        idx = np.argpartition(data,-number)[-number:]
        idx.sort()
        return idx

    def __top_of_segments(self,data,npoint,q):
        ids = (data > np.percentile(data,q)).nonzero()[0]
        nregion = len(ids)/npoint
        perreg  = npoint/nregion
        regions = np.array_split(ids,nregion)
        idx = np.array( [])
        # Top perreg of data in each region
        # have global indices given by reg, 
        # top_number returns indiced local to data[reg]
        for reg in regions:
             idx = np.hstack( (idx, reg[top_number(data[reg],perreg)]) )
        return idx.astype("int")

    def auto_fit_frame(self,method,ref):
        idx = self.auto_select_bgwin(method,ref)
        k,m = np.polyfit(self.frame.group.lmbd[idx],self.frame.data.transpose()[idx,:],1)
        return {"k":k,"m":m}
            
    def auto_select_bgwin(self,method,data, npoint=100,q=50):
        if   method == "top 100":
            return top_number(data,100)
        elif method == "top 5%":
            return top_number(data,round(len(data)*0.05))
        elif method == "top 20":
            return top_number(data,20)
        elif method == "90-95 decile":
            t5  = top_number(data,round(len(data)*0.05))[-1]
            t10 = top_number(data,round(len(data)*0.1))
            return t10[ t10 >= t5 ]
        elif method == "segments":
            return top_of_segments(data,npoint,q)
  
  
class refcontinua(continua):
    def __init__(self,group,method):
        self.group = group
        self.fit   = self.auto_fit_frame(method)

    def cont(self):
        return self.fit["k"]*self.group.lmbd + self.fit["m"]

    def auto_fit_frame(self,method):
        idx = self.auto_select_bgwin(method,self.group.ref)
        k,m = np.polyfit(self.group.lmbd[idx],self.group.ref[idx],1)
        return {"k":k,"m":m}
