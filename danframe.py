import astropy.io.fits as fits
import kontin as con
import pickle as pic
import numpy as np
import glob as g
import os

class frameseries(object):
    __lmbdname = "{}_{}__lambda.fits"
    __refname  = "{}_{}__adjustfts.fits"
    __savename = "{}_{}__metadata"

    def __init__(self,fileglob,method,rows=800):
        self.__parse_fileglob(fileglob)
        self.meta  = self.__load_meta()
        self.ref   = self.__load_from_fits(self.Dir+self.__refname)
        self.lmbd  = self.__load_from_fits(self.Dir+self.__lmbdname)
        self.files = g.glob(self.glob) ; self.files.sort()
        self.rows  = list(range(0,rows))
        self.refcon = con.refcontinua(self,method)
        self.ref   = self.ref/self.refcon.cont()

        self.veto_rows([0,799])
        self.normed = False

        try:
            self.pkwindows = self.meta["peakwin"] 
        except KeyError:
            self.pkwindows = []

        self.frames = []
        for fil in self.files:
            self.frames.append(
                danframe(fil,self,method))

    def __parse_fileglob(self,fileglob):
        self.glob = fileglob +"_[1-9]*" 
        if fileglob.find("/") > 0:
            self.Dir  = "/".join(fileglob.split("/")[:-1])+"/"  # everything before the last '/'
        else:
            self.Dir  = ""
        filename    = self.glob.split("/")[-1]
        self.series = filename.split("_")[1]
        self.wave   = filename.split("_")[0]

    def __load_from_fits(self,filename,hdu=0):
        fit = fits.open(filename.format(self.wave,self.series))
        return fit[hdu].data

    def normalize(self):
        if not self.normed:
            for frame in self.frames:
                frame.data = frame.data/frame.cont.norm()
        self.normed = True

    def veto_rows(self,rows):
        if not isinstance(rows,list):
            rows = [rows]
        for itm in rows:
            try:
                self.rows.remove(itm)
            except ValueError:
                print("Row {} not found".format(itm))

    def veto_and_update_rows(rows):
        if len(rows) > 0:
            self.veto_rows(rows)
        for frm in self.frames:
            frm.update_veto()            

    def set_bgwindows(self,bgwindows,warn=True):
        if len(self.bgwindows) > 0 and warn :
            print("WARNING: Old values will be erased")
            if not input("Continue Y/N? ").lower() == "y":
                return None
        bgwindows.sort()
        self.__set_meta(bgwindows,"bgwin")
        self.__load_bgwindows()

    def set_pkwindows(self,pkwindows,warn=True):
        if len(self.pkwindows) > 0 and warn :
            print("WARNING: Old values will be erased")
            if not input("Continue Y/N? ").lower() == "y":
                return None
        pkwindows.sort()
        self.__set_meta(pkwindows,"peakwin")
        self.pkwindows =self.__load_meta("peakwin")

    def __set_meta(self,metadata,keyword=""):
        with open(self.Dir+self.__savename.format(self.wave,self.series),"r+b") as metafile:
            try:
                meta =  pic.load(metafile)
            except EOFError:
                meta = {}
        if keyword != "":
            meta[keyword] = metadata
        else:
            meta = metadata
        with open(self.Dir+self.__savename.format(self.wave,self.series),"wb") as metafile:
            pic.dump(meta,metafile,protocol=1)

    def __load_meta(self,keyword=""):
        try:        
            with open(self.Dir+self.__savename.format(self.wave,self.series),"r+b") as metafile:
                meta =  pic.load(metafile)
        except EOFError:
            return {}
        except IOError:
            fil = open(self.Dir+self.__savename.format(self.wave,self.series),"w")
            fil.close()
            return {}
        if keyword == "":
            return meta
        else:
            return meta[keyword]


class danframe(object):
    def __init__(self,filename,group,method):
        self.fits   = fits.open(filename,mode="update")
        self.name   = filename.split("/")[-1]
        self.group  = group
        self.rdata  = self.fits[0].data 
        self.data   = self.rdata[self.group.rows,:]
        self.header = self.fits[0].header
        self.cont   = con.continua(self,method)
        #print "Frame {} has dimensions {}".format(self.name,self.rdata.shape)

    def update_veto(self):
        self.data   = self.rdata[self.group.cols,:]

    def spec(self,row):
        return self.data[row,:]

    def col(self,col):
        return self.data[:,col]

    def specplot(self,yval):
        pl.plot(self.lmbd,yval)
        pl.show()

class danframe_sac(object):
    __meanname = "{}_{}__adjustspec.fits"
    __lmbdname = "{}_{}__lambda.fits"
    __refname  = "{}_{}__adjustfts.fits"

    def __init__(self,filename):
        self.filename = filename.split("/")[-1]
        if filename.find("/") > 0:
            self.Dir  = "/".join(filename.split("/")[:-1])+"/"  # everything before the last '/'
        else:
            self.Dir  = ""

        self.fits   = fits.open(self.Dir+self.filename,mode="update")
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
        fit = fits.open(self.Dir+filename.format(self.wave,self.series))
        return fit[hdu].data

    def spec(self,row):
        return self.data[row,:]

    def col(self,col):
        return self.data[:,col]

    def specplot(self,yval):
        pl.plot(self.lmbd,yval)
        pl.show()

    def set_bgwindows(self,bgwindows,warn=True):
        if len(self.bgwindows) > 0 and warn :
            print("WARNING: Old values will be erased")
            if not input("Continue Y/N? ").lower() == "y":
                return None
        self.__set_windows(bgwindows,"BGWIN{}".format)
        self.__load_bgwindows()

    def set_pkwindows(self,pkwindows,warn=True):
        if len(self.pkwindows) > 0 and warn :
            print("WARNING: Old values will be erased")
            if not input("Continue Y/N? ").lower() == "y":
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
