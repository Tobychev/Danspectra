import matplotlib.pyplot as pl
import matplotlib.widgets as wdg 
import lines as lin

window  = []
windows = []
def select_manual_windows(danspec):
        global windows
        print "Selecting windows"

        windows = []
        fig = pl.figure()
        ax  = fig.add_subplot(111)
        ax.plot(danspec.mean)
        ax.plot([0,len(danspec.mean)],[1,1],'r')
        axis = wdg.AxesWidget(ax)
        axis.connect_event("key_press_event",__bg_window_select)
        cursor = wdg.Cursor(ax,useblit=True, color='red', 
                        linewidth=1,horizOn=False)
        pl.show()
        return windows

def __bg_window_select(event):
    global window, windows
    if event.key == "?":
        print """ 
a - add point to window, after two points have been indicated a new window is defined
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

def manual_delete_minima(mins,group):
    outlist = list(mins)
    for itm in mins:
        idx = list(itm)
        fig = pl.figure()
        ax  = fig.add_subplot(111)
        ax.plot(group.lmbd,group.ref)
        ax.plot(group.lmbd[idx],group.ref[idx],'ro-')
        axis = wdg.AxesWidget(ax)
        axis.connect_event("key_press_event",lambda x: __minima_delete(outlist,itm,x))
        pl.show()

    return outlist


def __minima_delete(minlist,item,event):
    if event.key == "?":
        print """
d - delete this minima
k - keep this minima
"""
    elif event.key == "d":
        minlist.remove(item)
        print "Minima {} removed".format(item)
        pl.close()
    elif event.key == "k":
        print "Minima {} saved".format(item)
        pl.close()


def pick_extrema_scatteplot(xs,ys):
    fig = pl.figure()
    ax  = fig.add_subplot(1,1,1)
    ax.plot(xs,ys,'ro',alpha=0.1,picker=5)
    fig.canvas.mpl_connect("pick_event",__onpick)
    pl.show()

def __onpick(event):
    curve = event.artist
    xs, ys = curve.get_xdata(),curve.get_ydata()
    idx = event.ind
    print event.ind, xs[idx],ys[idx]
