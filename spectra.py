import astropy.io.fits as fits
import pickle as pic
import numpy as np
import glob as g
import os

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

class SpecMeta(object):
    pass

class SpectraFactory(object):
    __lmbdname = "{}_{}__lambda.fits"
    __refname  = "{}_{}__adjustfts.fits"
    __savename = "spec_{}_{}.metadata"
    __contrast = "{}_{}__concont.fits"

    def __init__(self,fileglob,framerows=800):
        self.__parse_fileglob(fileglob)
        self.rows     = list(range(0,framerows))
        self.files    = np.array( g.glob(self.glob)); self.files.sort()
        self.contrast = self.__load_from_fits(self.Dir+self.__contrast)
        self.__load_meta()     
        self.meta["state"] = "new"

    def __load_from_fits(self,filename,hdu=0):
        # Works for proper filenames too because
        # format does noting unless 'filename' contains '{}'
        fit = fits.open(filename.format(self.wave,self.series)) 
        return fit[hdu].data

    def __parse_fileglob(self,fileglob):
        self.glob = fileglob +"_[1-9]*" 
        if fileglob.find("/") > 0:
            self.Dir  = "/".join(fileglob.split("/")[:-1])+"/"  # everything before the last '/'
        else:
            self.Dir  = ""
        filename    = self.glob.split("/")[-1]
        self.series = filename.split("_")[1]
        self.wave   = filename.split("_")[0]

    def __load_meta(self,keyword=""):
        """
            Tricksy function, when called with keyword it will return a value,
            else it updates self.meta as a side effect. Maybe not a good idea?
        """
        try:        
            with open(self.Dir+self.__savename.format(self.wave,self.series),"r+b") as metafile:
                meta =  pic.load(metafile)
                if keyword == "":
                    self.meta = meta
                else:
                    return meta[keyword]
        except EOFError:
            self.meta = {}
        except IOError:
            fil  = self.Dir+self.__savename.format(self.wave,self.series)
            print("Creating metafile {}".format(fil))
            fil = open(fil,"w")
            fil.close()
            self.meta = {}

    def __set_meta(self,metavalue,keyword=""):
        if keyword != "":
            meta[keyword] = metavalue
        with open(self.Dir+self.__savename.format(self.wave,self.series),"wb") as metafile:
            pic.dump(self.meta,metafile,protocol=1)
            
    def __update_meta(self):
        self.__set_meta(None)
        self.__load_meta()

    def contrast_cut(self,cutval,mode="percentile"):
        """
            Discard frames based on continuum contrast
            cutval - integer or float respectively
            mode   - one of "ordinal" or "percentile" 
        """

        if mode == "percentile":
            idx, = np.where(self.contrast > np.percentile(self.contrast,cutval))
        elif mode == "ordinal":
            idx  = self.contrast.argsort()[-cutval:]

        print("Keeping frames")
        print(self.files[idx])
        print("With continuum contrast")
        print(self.contrast[idx])

        self.files = self.files[idx]
        if self.meta["state"] == "new":
            self.meta["state"]   = "continuum cut"
        else:
            self.meta["state"]   = self.meta["state"] + "+continuum cut"
        if "contcut" in self.meta:
            self.meta["contcut"] = self.meta["contcut"] + "and Mode:{}, cutval:{}".format(mode,cutval)
        else:
            self.meta["contcut"] = "Mode:{}, cutval:{}".format(mode,cutval)
        self.__update_meta()

    def frame_row_cut(self,cutrows):
        if not isinstance(cutrows,list):
            cutrows = list(cutrows)
        for itm in cutrows:
            try:
                self.rows.remove(itm)
            except ValueError:
                print("Row {} not found".format(itm))
        if self.meta["state"] == "new":
            self.meta["state"] = "row cut"
        else:
            self.meta["state"] = self.meta["state"] + "+row cut"

        if "row cut" in self.meta:
            self.meta["rowtcut"] = self.meta["rowtcut"]+ "and {}".format(cutrows)
        else:
            self.meta["rowtcut"] = "{}".format(cutrows)
        self.__update_meta()


    def set_continua(self,method,nump=100,q=50):
        if method in ["top 20","segments"]:
            self.meta["cont method"] = method
            self.meta["nump"]   = nump
            self.meta["q"]      = q
            self.meta["state"]  = "Continua defined"
            self.__update_meta()
        else:
            print("Unrecognized method:", method)
        
    def make_block(self,desc=""):
        lmbd = self.__load_from_fits(self.Dir+self.__lmbdname)
        ref  = self.__load_from_fits(self.Dir+self.__refname )

        data = self.__load_from_fits(self.files[0])
        block = data[self.rows,:]
        for fil in self.files[1:]:
            data  = self.__load_from_fits(fil)
            block = np.vstack((block,data[self.rows,:]))
        
        if "cont method" in self.meta:
            con = continua(ref,lmbd,self.meta["cont method"],self.meta["nump"],self.meta["q"])
            block = block/con(lmbd,block)        
        if desc == "":
            desc = ", ".join( ": ".join((str(k),str(v))) for k,v in self.meta.items())

        return Spectra(desc,lmbd,block,self.meta)

class continua(object):
    def __init__(self,refdata,lmbd,method,nump=30,q=80):
        self.idx      = self.__def_continua(refdata,method,nump,q)
        self.lmbd     = lmbd[self.idx]

    def __call__(self,lmbd,data):
        if len(data.shape) == 1:
            k,m = np.polyfit(self.lmbd,data[self.idx],1)
            return k*lmbd+m   
        elif len(data.shape) == 2:
            k,m = np.polyfit(self.lmbd,data[:,self.idx].T,1)
            return k.reshape(-1,1)*lmbd+m.reshape(-1,1) 
        else:
            raise ValueError("Data must be 1 or 2d")
        
    def __top_of_segments(self,data,npoint,q):
        ids, = np.where(data > np.percentile(data,q))
        nregion = int(len(ids)/npoint)
        perreg  = int(npoint/nregion)
        regions = np.array_split(ids,nregion)
        idx = []
        # Top perreg of data in each region
        # have global indices given by reg, 
        # selection returns indices local to data[reg]
        for reg in regions:
            idx.append( reg[ data[reg].argsort()[-perreg:] ])
        idx = np.array(idx).reshape(-1)
        return idx.astype("int")
            
    def __def_continua(self,data,method,nump,q):
        if   method == "top 20":
            return data.argsort()[-20:]
        elif method == "segments":
            return self.__top_of_segments(data,nump,q)
        elif   method == "top N":
            return data.argsort()[-q:]
