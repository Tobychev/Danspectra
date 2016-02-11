import numpy as np
import pyfits as f
import os

class danspectra(object):
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
        keyword = "BGWIN{}".format
        LOW = 0; HIGH = 1 # Lower and upper window edge
    
        del self.header[keyword("?*")]
        self.header.set(keyword("N"),len(bgwindows))
   
        bgwindows.sort() # For nicer order in header
        for i in range(0,len(bgwindows)):
            self.header.set(keyword(str(i+1)+"L"), bgwindows[i][LOW])   #Want index 1,..,N
            self.header.set(keyword(str(i+1)+"H"), bgwindows[i][HIGH])
    
        self.fits.flush() 
        self.__load_bgwindows()
        self.haswindows = True

    def __load_bgwindows(self):
        keyword = "BGWIN{}".format
        wins = self.header.get("BGWINN") 
        window = []
        self.bgwindows = []
        for i in range(1,wins+1):   #Want index 1,..,N
            window = []
            window.append(self.header.get(keyword(str(i)+"L")))
            window.append(self.header.get(keyword(str(i)+"H")))
            self.bgwindows.append(window)
