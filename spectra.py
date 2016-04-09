import numpy as np

class Spectra(object):

    def __init__(self,desc,lmbd,datablock,meta):
        if len(lmbd.shape) > 1:
            raise ValueError("Only 1-d lambda vectors please.")
        if lmbd.shape[0] != datablock.shape[1]:
            raise ValueError("Datablock must have same lenght in second dimension as lambda.")

        self.__data = np.copy(datablock)
        self.__data.flags.writeable = False
        self.__description = str(desc)
        self.meta = meta 
        self.lmbd = lmbd

    def __getitem__(self,key):
        return self.__data[key]

    def __repr__(self):
        return "[{},{}] Datablock: ".format(self.__data.shape[0],self.__data.shape[1])+self.__description

class MetaSpectra(object):
    pass

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
    
