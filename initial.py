from matplotlib.widgets import Cursor,AxesWidget
import matplotlib.pyplot as pl
import numpy as np
import pyfits as f
import os

class danspectra(object):
    Dir = "data/"
    meanname = "{}_{}__meanspec.fits"
    lmbdname = "{}_{}__lambda.fits"
    bgwindows = []
    __window  = []
    def __init__(self,filename):
        self.data    = f.open(self.Dir+filename)[0].data
        self.series  = filename.split("_")[1]
        self.wave    = filename.split("_")[0]
        self.lmbd    = f.open(self.Dir+self.lmbdname.format(
                                        self.wave,self.series))
        self.lmbd    = self.lmbd[0].data
        self.mean     = f.open(self.Dir+self.meanname.format(
                                        self.wave,self.series))
        self.mean    = self.mean[0].data
        self.isnorm  = False

    def residue(self,row,norm=1):
        return self.data[row,:] - self.mean/norm

    def specplot(self,yval):
        pl.plot(self.lmbd,yval)
        pl.show()
    
    def select_background(self):
        print "Defining new background windows"
        self.bgwindows = []
        fig = pl.figure()
        ax  = fig.add_subplot(111)
        ax.plot(self.lmbd,self.mean)
        axis = AxesWidget(ax)
        axis.connect_event("key_press_event",self.manageux)
        cursor = Cursor(ax,useblit=True, color='red', 
                        linewidth=1,horizOn=False)
        pl.show()        
        print self.bgwindows

    def manageux(self,event):
        if event.key == "?":
            print """ 
a - add point to background window, after two points have been indicated a new window is defined
"""            
        elif event.key == "a":
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
# lines - spectra
# rows  - dispersion/space axis
fil = "6405_aS1_388_cor.fits"
spec = danspectra(fil)



