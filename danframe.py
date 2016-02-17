import numpy as np
import pyfits as f
import os

class danframe(object):
    __meanname = "{}_{}__adjustspec.fits"
    __lmbdname = "{}_{}__lambda.fits"
    __refname  = "{}_{}__adjustfts.fits"

    def __init__(self,filename):
        self.filename = filename.split("/")[-1]
        if filename.find("/") > 0:
            self.Dir  = "/".join(filename.split("/")[:-1])+"/"  # everything before the last '/'
        else:
            self.Dir  = ""

        self.fits   = f.open(self.Dir+self.filename,mode="update")
        self.data   = self.fits[0].data
        self.header = self.fits[0].header
        self.series = self.filename.split("_")[1]
        self.wave   = self.filename.split("_")[0]
        self.lmbd   = self.__load_from_fits(self.__lmbdname)
        self.mean   = self.__load_from_fits(self.__meanname)
        self.ref    = self.__load_from_fits(self.__refname)
        if self.header.get("BGWINN") > 0:
            self.__load_bgwindows()
        else:
            self.bgwindows  = []
        if self.header.get("PEAKWINN") > 0:
            self.__load_pkwindows()
        else:
            self.pkwindows  = []

    def __load_from_fits(self,filename,hdu=0):
        fits = f.open(self.Dir+filename.format(self.wave,self.series))
        return fits[hdu].data

    def spec(self,row):
        return self.data[row,:]

    def col(self,col):
        return self.data[:,col]

    def specplot(self,yval):
        pl.plot(self.lmbd,yval)
        pl.show()

    def set_bgwindows(self,bgwindows,warn=True):
        if len(self.bgwindows) > 0 and warn :
            print "WARNING: Old values will be erased"
            if not raw_input("Continue Y/N? ").lower() == "y":
                return None
        self.__set_windows(bgwindows,"BGWIN{}".format)
        self.__load_bgwindows()

    def set_pkwindows(self,pkwindows,warn=True):
        if len(self.pkwindows) > 0 and warn :
            print "WARNING: Old values will be erased"
            if not raw_input("Continue Y/N? ").lower() == "y":
                return None
        self.__set_windows(pkwindows,"PEAKWIN{}".format)
        self.__load_pkwindows()

    def __set_windows(self,windows,keyword):
        LOW = 0; HIGH = 1 # Lower and upper window edge
    
        del self.header[keyword("?*")]
        self.header.set(keyword("N"),len(windows))
   
        windows.sort() # For nicer order in header
        for i in range(0,len(windows)):
            self.header.set(keyword(str(i+1)+"L"), windows[i][LOW])   #Want index 1,..,N
            self.header.set(keyword(str(i+1)+"H"), windows[i][HIGH])
    
        self.fits.flush() 

    def __load_bgwindows(self):
        self.bgwindows = self.__load_windows("BGWINN","BGWIN{}".format)

    def __load_pkwindows(self):
        self.pkwindows = self.__load_windows("PEAKWINN","PEAKWIN{}".format)

    def __load_windows(self,head,keyword):
        wins = self.header.get(head) 
        window  = []
        windows = []
        for i in range(1,wins+1):   #Want index 1,..,N
            window = []
            window.append(self.header.get(keyword(str(i)+"L")))
            window.append(self.header.get(keyword(str(i)+"H")))
            windows.append(window)
        return windows

