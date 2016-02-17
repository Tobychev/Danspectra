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
        idbg = idbg + range(win[0],win[1]+1)
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
            print idx[i]+1, idx[i+1]
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
        return top_number(data,round(len(data)*0.1))
    elif metod == "90-95 decile":
        t5  = top_number(data,round(len(data)*0.05))[-1]
        t10 = top_number(data,round(len(data)*0.1))
        return t10[ t10 >= t5 ]
    elif metod == "ref top":       
        return range( danframe.ref.argmax()- 5, danframe.ref.argmax()+ 5 )
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
    def __init__(self,danframe,extras=[]):
        self.frame = danframe
        self.top20 = self.auto_fit_frame("top 20")
        self.top5p = self.auto_fit_frame("top 5%")
        self.segmt = selt.auto_fit_frame("segments")
        if len(extras) > 0:
            self.extra = {}
            for method in extras:           
                self.extra[method] = self.auto_fit_frame(method)

    def top20cont(self,wavlen):
        return self.top20["k"]*wavlen + self.top20["m"]

    def top5pcont(self,wavlen):
        return self.top5p["k"]*wavlen + self.top5p["m"]

    def segmtcont(self,wavlen):
        return self.segmt["k"]*wavlen + self.segmt["m"]

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

    def auto_fit_frame(self,method):
        idx = self.auto_select_bgwin(method)
        k,m = np.polyfit(self.frame.lmbd[idx],self.frame.data.transpose()[idx,:],1)
        return {"k":k,"m":m}
            
    def auto_select_bgwin(self,method, npoint=100,q=50):
        data = self.frame.ref

        if   metod == "top 100":
            return top_number(data,100)
        elif metod == "top 5%":
            return top_number(data,round(len(data)*0.05))
        elif metod == "top 20":
            return top_number(data,round(len(data)*0.1))
        elif metod == "90-95 decile":
            t5  = top_number(data,round(len(data)*0.05))[-1]
            t10 = top_number(data,round(len(data)*0.1))
            return t10[ t10 >= t5 ]
        elif metod == "segments":
            return top_of_segments(data,npoint,q)

