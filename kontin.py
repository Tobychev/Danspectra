from matplotlib.widgets import Cursor,AxesWidget
import matplotlib.pyplot as pl
import numpy as np

def make_idx_from_windows(windows):
    idbg = []
    for win in windows:
        idbg = idbg + range(win[0],win[1]+1)
    return idbg

def select_background_automatic(danspec,metod,line="mean"):
    if line == "mean":
        data = danspec.mean
    else:
        data = danspec.spec(line)

    if metod == "tophundred":
        return top_number(data,100)
    if metod == "top 5%":
        return top_number(data,round(len(data)*0.05))
    if metod == "top 10%":
        return top_number(data,round(len(data)*0.1))
    if metod == "90-95 decile":
        t5  = top_number(data,round(len(data)*0.05))[-1]
        t10 = top_number(data,round(len(data)*0.1))
        return t10[ t10 >= t5 ]

def top_number(data,number):
    idx = np.argpartition(data,-number)[-number:]
    idx.sort()
    return idx
    
window  = []
windows = []
def select_background_manual(danspec):
        print "Selecting background windows"

        if len(danspec.bgwindows) > 0:
            print "WARNING: background windows already defined"
            print "WARNING: all old values will be erased"
            if not raw_input("Continue Y/N? ").lower() == "y":
                return None
        windows = []
        fig = pl.figure()
        ax  = fig.add_subplot(111)
        ax.plot(danspec.mean)
        ax.plot([0,len(danspec.mean)],[1,1],'r')
        axis = AxesWidget(ax)
        axis.connect_event("key_press_event",__manageux)
        cursor = Cursor(ax,useblit=True, color='red', 
                        linewidth=1,horizOn=False)
        pl.show()
        if len(windows) > 0:
            danspec.set_bgwindows(windows,warn=False) #Already asked

        return make_idx_from_windows(windows)


def __manageux(event):
    global window, windows
    if event.key == "?":
        print """ 
a - add point to background window, after two points have been indicated a new window is defined
"""            
    elif event.key == "a":
        pos = int(round(event.xdata))
        print "adding ", pos
        if len(window) == 0:
            window.append(pos)
        elif len(window) == 1:
            if window[0] < event.xdata:
                window.append(pos)
            else:
                window.append(window[0])
                window[0] = pos
            windows.append(window)
            window=[]
