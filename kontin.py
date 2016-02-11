from matplotlib.widgets import Cursor,AxesWidget
import matplotlib.pyplot as pl
import numpy as np

def frame_fit(spec,idx):
   return np.polyfit(spec.lmbd[idx],spec.data.transpose()[idx,:],1)

def frame_subtract_cont(spec,idx):
    fits = frame_fit(spec,idx)
    cont = fits[1,:].reshape(len(fits[1,:]),1) + np.outer(fits[0,:],spec.lmbd)
    return spec.data - cont
    
def row_subtract_cont(spec,line,idx):
    b,a = row_fit_continuum(spec,line,idx)
    cont = a + b*spec.lmbd
    return spec.spec(line)-cont 

def row_fit_continuum(spec,line,idx):
    fit = np.polyfit(spec.lmbd[idx],spec.spec(line)[idx],1)
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

def select_bgwin_auto(danspec,metod,row="mean",
                        npoint=100,q=50):
    if row == "mean":
        data = danspec.mean
    elif row == "ref":
        data = danspec.ref
    else:
        data = danspec.spec(row)

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
        return range( danspec.ref.argmax()- 5, danspec.ref.argmax()+ 5 )
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

def save_bgwin_from_idx(danspec,idx):
    wins = make_windows_from_idx(idx)
    danspec.set_bgwindows(wins)

