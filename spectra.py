import astropy.io.fits as fits
import numpy.polynomial.polynomial as pol
import scipy.interpolate as si
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
        return "[{},{}] ".format(self.__data.shape[0],self.__data.shape[1])+self.__description

    def modify(self,mutator):
        data = mutator(self.__data)
        assert(data.shape == self.__data.shape,"Data has been reshaped")
        self.__data.flags.writeable = True
        self.__data = data
        self.__data.flags.writeable = False
        print("Spectra updated")

class SpecMeta(object):
    __savename = "spec_{}_{}.metadata"    
    def __init__(self,filename,cont,lmbd,ref):
        self.__parse_filename(filename)
        self._load_meta()
        self.cont = cont
        self.lmbd = lmbd
        self.ref  = ref

    def __parse_filename(self,filename):
        self.glob = filename +"_[1-9]*" 
        if filename.find("/") > 0:
            self.Dir  = "/".join(filename.split("/")[:-1])+"/"  # everything before the last '/'
        else:
            self.Dir  = ""
        filname = filename.split("/")[-1]
        self.wave   = filename.split("_")[1]
        self.series = filename.split("_")[2].split(".")[0]

    def _load_meta(self,keyword=""):
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

    def _set_meta(self,metavalue,keyword=""):
        if keyword != "":
            meta[keyword] = metavalue
        with open(self.Dir+self.__savename.format(self.wave,self.series),"wb") as metafile:
            pic.dump(self.meta,metafile,protocol=1)
            
    def _update_meta(self):
        self._set_meta(None)
        self._load_meta()

class SpectraFactory(SpecMeta):
    __lmbdname = "{}_{}__lambda.fits"
    __refname  = "{}_{}__adjustfts.fits"
    __savename = "spec_{}_{}.metadata"
    __contrast = "{}_{}__concont.fits"

    def __init__(self,fileglob,framerows=800):
        self.__parse_fileglob(fileglob)
        self.rows     = list(range(0,framerows))
        self.files    = np.array( g.glob(self.glob)); self.files.sort()
        self.contrast = self.__load_from_fits(self.Dir+self.__contrast)
        self._load_meta()
        self.meta["state"] = "new"

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
        # Works for proper filenames too because
        # format does noting unless 'filename' contains '{}'
        fit = fits.open(filename.format(self.wave,self.series)) 
        return fit[hdu].data

    def __make_block(self):
        data = self.__load_from_fits(self.files[0])
        block = data[self.rows,:]
        for fil in self.files[1:]:
            data  = self.__load_from_fits(fil)
            block = np.vstack((block,data[self.rows,:]))
        return block
    
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
        self.meta["contcut"] = "Mode:{}, cutval:{}".format(mode,cutval)
        self._update_meta()

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

        self.meta["rowtcut"] = "{}".format(cutrows)
        self._update_meta()

    def set_continua(self,method,nump=30,q=80):
        if method in ["top 20","segments"]:
            self.meta["cont method"] = method
            self.meta["nump"]   = nump
            self.meta["q"]      = q
            self.meta["state"]  = "Continua defined"
            self._update_meta()
        else:
            print("Unrecognized method:", method)

    def set_continua(self,method,nump=30,q=80):
        if method in ["top 20","segments"]:
            self.meta["cont method"] = method
            self.meta["nump"]   = nump
            self.meta["q"]      = q
            self.meta["state"]  = "Continua defined"
            self._update_meta()
        else:
            print("Unrecognized method:", method)

    def make_spectra(self,desc=""):
        lmbd = self.__load_from_fits(self.Dir+self.__lmbdname)
        ref  = self.__load_from_fits(self.Dir+self.__refname )

        block = self.__make_block()
                
        if "cont method" in self.meta:
            con   = continua(ref,lmbd,self.meta["cont method"],self.meta["nump"],self.meta["q"])
            cont  = con.fit(block)
            block = block/con(lmbd,block)
        else: 
            con = None
        if desc == "":
            desc = ", ".join( ": ".join((str(k),str(v))) for k,v in self.meta.items())

        meta = SpecMeta(self.Dir+self.__savename.format(self.wave,self.series),cont,lmbd,ref)
        return Spectra(desc,lmbd,block,meta)

    def make_spectra_subset(self,spectra,rowsubset=None,colsubset=None,desc=""):
        shape = spectra[:,:].shape
        colid = range(0,shape[1]); rowid = range(0,shape[0])
        meta  = spectra.meta
        subset = None

        if desc == "":
            desc = "Subset of {}".format(spectra)
        if rowsubset is not None:
            if rowsubset.dtype == np.dtype('bool'):
                rowid, = np.where(rowsubset)
            else:
                if len(rowsubset) > shape[0]:
                    raise IndexError("Length of rowsubset grater than number of rows in spectra")
                rowid = rowsubset

        if colsubset is not None:
            if colsubset.dtype == np.dtype('bool'):
                colid, = np.where(colsubset)
            else:
                if len(colsubset) > shape[0]:
                    raise IndexError("Length of colsubset grater than number of cols in spectra")
                colid = colsubset

        if rowsubset is not None or colsubset is not None:
            subset = spectra[rowid,:]
            subset = subset[:,colid]
            return Spectra(desc,spectra.meta.lmbd[colid],subset,meta)
        else:
            print("No selection given")

class continua(object):
    def __init__(self,refdata,lmbd,method,nump=30,q=80):
        self.idx      = self.__def_continua(refdata,method,nump,q)
        self.lmbd     = lmbd[self.idx]

    def fit(self,data):
        if len(data.shape) == 1:
            k,m = np.polyfit(self.lmbd,data[self.idx],1)
            return k,m   
        elif len(data.shape) == 2:
            k,m = np.polyfit(self.lmbd,data[:,self.idx].T,1)
            return k,m 
        else:
            raise ValueError("Data must be 1 or 2d")

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

class line(object):
    def __init__(self,winbounds,cent,specmeta):
        self.idx   = np.arange(winbounds[0],winbounds[1]+1)
        self.cent  = cent
        self.name  = "{:6.3f}".format(self.cent)
        self.width = specmeta.lmbd[self.idx[0]] - specmeta.lmbd[self.idx[-1]] 
        self.spec  = specmeta

    def __repr__(self):
        return "Line {} [$\Delta \lambda$ = {:.4f} Å]".format(self.name,self.width/10)        

    def _equivalent_width(self,spec):
        dlam = np.diff(spec.meta.lmbd[slice(self.idx[0]-1,self.idx[-1]+1)]
                      ).reshape((-1,1))*np.ones(spec[:,0].shape)
        return ((spec[:,self.idx]-1)*dlam.T).sum(axis=1)*1e3 ## MiliÅngström

class statline(line):
    def measure(self,spectra):
        guess = spectra.meta.ref[self.idx].argmin()
        width = len(self.idx)*0.16 # Min fraction of points to be used in fit

        if width%2 == 0:
            width +=1
        dwn = int((width - 1)/2); up = dwn+1

        bottom  = self.idx[guess + np.arange(-dwn,up)]   
        test    = self.idx[guess + np.arange(-(dwn+1),(up+1))]

        ew = self._equivalent_width(spectra)
        vel,bot,err    = self.__linfit(spectra,bottom,test)
        mn,var,ske,kur = self.__moments(spectra)
        wvar,wske,wkur = self.__ew_moments(spectra,vel,bot,ew)        
        
        if spectra.meta.cont is not None:
            con = spectra.meta.cont[0].reshape(-1,1)*self.cent + spectra.meta.cont[1].reshape(-1,1)
        else:
            con = None

        return vel,bot,con,err,ew,mn,var,ske,kur,wvar,wske,wkur

    def __linfit(self,spec,bottom,test):
        cv = 299792.458
        fit   = pol.polyfit(spec.meta.lmbd[bottom],spec[:,bottom].T,2)
        a,b,c = fit[2,:],fit[1,:],fit[0,:]
        lmin  = -b/(2*a)
        bot   = pol.polyval(lmin,fit,tensor=False)
        pred  = pol.polyval(spec.meta.lmbd[test],fit)
        vel   = cv*(lmin-self.cent)/self.cent
        # Error of fit with extra error term to penalize fits 
        # that gets wildly off center, with extra weight so it *hurts*
        err   = np.sqrt( np.mean( (spec[:,test]-pred)**2,axis=1) 
                                         + 2*(lmin-self.cent)**2 )
        return vel,bot,err

    def __moments(self,spec):
        x    = spec.meta.lmbd[self.idx]
        dpdf = (1-spec[:,self.idx]/spec[:,self.idx].max(axis=1).reshape(-1,1))
        dpdf = dpdf/dpdf.sum(axis=1).reshape(-1,1)
        mu   = np.sum(dpdf*x,axis=1).reshape(-1,1) # Reshaping enables broadcasting
        mu2  = np.sum(dpdf*(x-mu)**2,axis=1)
        mu3  = np.sum(dpdf*(x-mu)**3,axis=1)
        mu4  = np.sum(dpdf*(x-mu)**4,axis=1)
        mu   = mu.reshape(-1) # Undoing reshape to allow assignment

        skew = mu3/mu2**(3/2) 
        kurt = (mu4/mu2**2 - 3)
        return mu,mu2,skew,kurt

    def __ew_moments(self,spec,vel,bot,ew):
        lmbd = spec.meta.lmbd[self.idx].reshape(-1,1)
        dlam = np.diff(spec.meta.lmbd[slice(self.idx[0]-1,self.idx[-1]+1)]).reshape((-1,1))*np.ones(spec[:,:].shape[0])
        ew   = ew*1e-3 # Cancels scaling

        # Variance by ratio between center and outer mass
        lsel = (lmbd > vel-ew/2) & (lmbd < vel+ew/2); 
        In   = ((spec[:,self.idx]-1)*dlam.T*lsel.T).sum(axis=1)
        var  =  In/ew

        # Skewness by ratio between left and right mass
        lsel = (lmbd < vel)
        lft  = ((spec[:,self.idx]-1)*dlam.T*lsel.T).sum(axis=1); rght = ((spec[:,self.idx]-1)*dlam.T*np.logical_not(lsel.T)).sum(axis=1);
        cut, = np.where(rght == 0); lft[cut] = 0; rght[cut] = 1
        ske  = lft/rght-1

        # Kurtosis 
        lsel = spec[:,self.idx] > (1 +   bot.reshape(-1,1))/2
        up   = ((spec[:,self.idx]-1)*dlam.T*lsel).sum(axis=1); dwn  = ((spec[:,self.idx]-1)*dlam.T*np.logical_not(lsel)).sum(axis=1);
        cut, = np.where(dwn == 0); dwn[cut] = 1; up[cut] = 0
        kur  = up/dwn

        return var, ske, kur

class splineline(line):
    def measure(self,spectra):
        nrows = spectra[:,:].shape[0]
        lmbd  = spectra.meta.lmbd[self.idx]
        ew = self._equivalent_width(spectra)

        splmes = np.zeros((nrows,11))
        splmes[:,10] = (spectra.meta.cont[0]*self.cent)+ (spectra.meta.cont[1])
        splmes[:, 9] = ew.reshape(-1)
        print("Making splines and measuring\n")
        for i,row in enumerate(spectra[:,self.idx]):
            mf           = self.makespline(row,lmbd,9)
            splmes[i,:9] = self.measure_spline(mf,lmbd)

        return splmes

    def makespline(self,spec,lmbd,kns=6):
        _,kno = np.histogram(lmbd,kns+2)
        kno   = kno[1:-2]
        return si.LSQUnivariateSpline(lmbd[::-1],spec[::-1],kno)

    def measure_spline(self,spl,lmbd,dl=2e-5):
        lmbd = np.linspace(lmbd[0],lmbd[-1],int( (lmbd[0]-lmbd[-1])/dl )) 
        bot  = spl(lmbd).min()
        icnt = spl(lmbd).argmin()
        cnt  = lmbd[icnt]
        bo12 = (1 +   bot)/2        
        bo13 = (1 + 2*bot)/3
        bo23 = (2 +   bot)/3
        fwhm,as12 = self.__width_assym(spl,lmbd,bo12,cnt)
        fw13,as13 = self.__width_assym(spl,lmbd,bo13,cnt)
        fw23,as23 = self.__width_assym(spl,lmbd,bo23,cnt)
        cnt = 299792.458*(cnt-self.cent)/self.cent

        return bot,cnt,fwhm,as12,fw13,as13,fw23,as23,spl.get_residual()

    def __width_assym(self,spl,lmbd,lev,cnt):
        ilev, = np.where(spl(lmbd) <= lev)
        # Check that we only got one interval
        spli, = np.where(np.diff(ilev) > 1)    # Either a number or empty
        if spli.sum() > 0:
            if   len(spli) == 1 :
                ilev  = ilev[slice(spli+1)]

        wdth  = lmbd[ilev[0]] - lmbd[ilev[-1]]
        assm  = cnt  - (lmbd[ilev[0]] + lmbd[ilev[-1]])/2
        return wdth,assm

