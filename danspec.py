from matplotlib.widgets import Cursor,AxesWidget
import matplotlib.pyplot as pl
import numpy as np
import pyfits as f
import os

class danspectra(object):
    Dir = "data/"
    __meanname = "{}_{}__meanspec.fits"
    __lmbdname = "{}_{}__lambda.fits"
    bgwindows = []
    __window  = []
    def __init__(self,filename):
        self.fits    = f.open(self.Dir+filename,mode="update")
        self.data    = self.fits[0].data
        self.series  = filename.split("_")[1]
        self.wave    = filename.split("_")[0]
        self.lmbd    = f.open(self.Dir+self.__lmbdname.format(
                                        self.wave,self.series))
        self.lmbd    = self.lmbd[0].data
        self.mean     = f.open(self.Dir+self.__meanname.format(
                                        self.wave,self.series))
        self.mean    = self.mean[0].data
        self.isnorm  = False
        if self.fits[0].header.get("BGWINN"):
            self.__load_bgwindows()
            self.haswindows = True
        else:
            self.haswindows = False
            

    def residue(self,row,norm=1):
        return self.data[row,:] - self.mean/norm

    def specplot(self,yval):
        pl.plot(self.lmbd,yval)
        pl.show()
    
    def select_background(self):
        print "Defining new background windows"
        if self.haswindows:
            print "WARNING: background windows already defined"
        self.bgwindows = []
        fig = pl.figure()
        ax  = fig.add_subplot(111)
        ax.plot(self.lmbd,self.mean)
        axis = AxesWidget(ax)
        axis.connect_event("key_press_event",self.__manageux)
        cursor = Cursor(ax,useblit=True, color='red', 
                        linewidth=1,horizOn=False)
        pl.show()
        self.__set_bgwindows() 

    def __manageux(self,event):
        if event.key == "?":
            print """ 
a - add point to background window, after two points have been indicated a new window is defined
"""            
        elif event.key == "a":
            print event.xdata, self.__window
            if len(self.__window) == 0:
                self.__window.append(event.xdata)
            elif len(self.__window) == 1:
                if self.__window[0] < event.xdata:
                    self.__window.append(event.xdata)
                else:
                    self.__window.append(self.__window[0])
                    self.__window[0] = event.xdata
                self.bgwindows.append(self.__window)
                self.__window=[]

    def __set_bgwindows(self):
        keyword = "BGWIN{}".format
        LOW = 0; HIGH = 1 # Lower and upper window edge
    
        head = self.fits[0].header
        head.set(keyword("N"),len(self.bgwindows))

        # Wavlength increases in the negative x-direction
        # So switching things around here makes it look right
        for i in range(0,len(self.bgwindows)):
            idxLow  = self.__find_nearest(self.bgwindows[i][LOW])
            idxHigh = self.__find_nearest(self.bgwindows[i][HIGH])
            print i,idxLow,idxHigh, self.bgwindows[i]
            
            head.set(keyword(str(i+1)+"L"), idxLow)   #Want index 1,..,N
            head.set(keyword(str(i+1)+"H"), idxHigh)
    
        self.fits.flush() 
        self.__load_bgwindows()
        self.haswindows = True

    def __load_bgwindows(self):
        keyword = "BGWIN{}".format
        head = self.fits[0].header
        wins = head.get("BGWINN") 
        self.bgwindows = []
        print "load" 
        for i in range(1,wins+1):   #Want index 1,..,N
            self.__window = []
            self.__window.append(head.get(keyword(str(i)+"L")))
            self.__window.append(head.get(keyword(str(i)+"H")))
            self.bgwindows.append(self.__window)

    def __find_nearest(self,value):
        return (np.abs(self.lmbd-value)).argmin()

