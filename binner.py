import numpy as np
import danframe as dan
import astropy.stats as ast
import numpy.polynomial.polynomial as pol
import scipy.interpolate as si

class binspec(object):
    def __init__(self,group,cont,cond,method='blocks'):
        """
        Averages spectra in [frame]group based on cond[ition]
            if cond has dtype bool it makes two bins
            else it histograms cond usin bayesian blocks finds average spectra per bin
        cont[inuum] - continuum level for each row in the framegroup 
        """
        self.group = group
        block      = group.frames[0].data
        for frm in group.frames[1:]:
            block = np.vstack((block,frm.data))
        if cond.dtype == np.dtype('bool'):
            idx,  = np.where(cond)
            nidx, = np.where(np.logical_not(cond))
            binned = np.mean( block[idx,:],axis=0)
            con    = np.mean( cont[idx])

            self.binned = np.vstack( (binned, np.mean( block[nidx,:],axis=0)) )
            self.con    = np.vstack( (con   , np.mean( cont[nidx])) )
            self.bins   = np.array([0,0,1])
            self.counts = np.array([cond.sum(),np.logical_not(cond).sum()])
        else:
            if len(cond) != block.shape[0]:
                raise IndexError("Length of cond does not match number of rows in framegroup")
            counts, bins = ast.histogram(cond,bins=method)
            sorting = np.digitize(cond,bins[:-1])  # :-1 to get the correct number of buckets from digitize
            binned  = np.mean(block[sorting == 1,:],axis=0)
            con     = np.mean( cont[sorting == 1])

            for i in np.unique(sorting)[1:]:
                binned = np.vstack( (binned, np.mean(block[sorting == i,:],axis=0) ))
                con    = np.vstack( (con   , np.mean( cont[sorting == i])          ))

            self.binned = binned
            self.con    = con
            self.counts = counts
            self.bins   = bins
            self.sorting = sorting

class binned_framegroup(object):

    def __init__(self,line,group,bin_quant,cuts=None):
        self.idx   = line.idx
        self.cent  = line.cent
        self.group = group
        self.data,self.cont,self.bins,self.counts = self.__bin_by_quant(bin_quant,cuts)
    
    def __bin_by_quant(self,quant,cuts):
        try:             
            cuts, = np.where(cuts.reshape(-1))
        except AttributeError:
            cuts = np.ones(len(quant))

        datablock = self.group.frames[0].data[:,self.idx]        
        contblock = self.group.frames[0].cont.val(self.cent)
        for frm in self.group.frames[1:]:
            datablock = np.concatenate((datablock,frm.data[:,self.idx]),axis=0)
            contblock = np.concatenate((contblock,frm.cont.val(self.cent)),axis=0)
        datablock = datablock[cuts,:]
        contblock = contblock[cuts]

        counts, bins = ast.histogram(quant[cuts],bins='blocks')
        sorting = np.digitize(quant[cuts],bins[:-1])  # :-1 to get the correct number of buckets from digitize
        binned  =  np.mean(datablock[sorting == 1,:],axis=0)
        con     =  np.zeros(len(counts)); con[0] = np.mean( contblock[sorting == 1] )
        for i in np.unique(sorting)[1:]:
            binned   = np.vstack( (binned,np.mean(datablock[sorting == i,:],axis=0) ))
            con[i-1] = np.mean(contblock[sorting == i])
        return binned,con,bins,counts

    def partition_data(self,data):
        return np.digitize(data,self.bins[:-1])  # :-1 to get the correct number of buckets from digitize

    def measure(self):
        nfram = len(self.data)
        guess = self.group.ref[self.idx].argmin()
        width = len(self.idx)*0.16 # Min fraction of points to be used in fit

        vel = np.zeros(nfram)
        bot = np.zeros(nfram)
        err = np.zeros(nfram)
        ew  = np.zeros(nfram)
        mn  = np.zeros(nfram)
        var = np.zeros(nfram)
        ske = np.zeros(nfram)
        kur = np.zeros(nfram)

        if width%2 == 0:
            width +=1
        dwn = int((width - 1)/2); up = dwn+1

        bottom  = guess + np.arange(-dwn,up) 
        test    = guess + np.arange(-(dwn+1),(up+1))

        vel,bot,err = self.__linfit(bottom,test)
        ew = self.__equivalent_width(nfram)
        mn,var,ske,kur = self.__moments()

        return (vel,bot,err,ew ,mn ,var,ske,kur)

    def __linfit(self,bottom,test):
        cv = 299792.458
        fit   = pol.polyfit(self.group.lmbd[self.idx[bottom]],self.data[:,bottom].T,2)
        a,b,c = fit[2,:],fit[1,:],fit[0,:]
        lmin  = -b/(2*a)
        bot   = pol.polyval(lmin,fit,tensor=False)
        pred  = pol.polyval(self.group.lmbd[test],fit)
        vel   = cv*(lmin-self.cent)/self.cent
        # Error of fit with extra error term to penalize fits 
        # that gets wildly off center, with extra weight so it *hurts*
        err   = np.sqrt( np.mean( (self.data[:,test]-pred)**2,axis=1) 
                                         + 2*(lmin-self.cent)**2 )
        return vel,bot,err

    def __equivalent_width(self,nrows):
        dlam = np.diff(self.group.lmbd[slice(self.idx[0]-1,self.idx[-1]+1)]
                       ).reshape((-1,1))*np.ones(nrows)
        return ((self.data-1)*dlam.T).sum(axis=1)*1e3 ## MiliÅngström

    def __moments(self):
        x    = self.group.lmbd[self.idx]
        dpdf = (1-self.data/self.data.max(axis=1).reshape(-1,1))
        dpdf = dpdf/dpdf.sum(axis=1).reshape(-1,1)
        mu   = np.sum(dpdf*x,axis=1).reshape(-1,1) # Reshaping enables broadcasting
        mu2  = np.sum(dpdf*(x-mu)**2,axis=1)
        mu3  = np.sum(dpdf*(x-mu)**3,axis=1)
        mu4  = np.sum(dpdf*(x-mu)**4,axis=1)
        mu   = mu.reshape(-1) # Undoing reshape to allow assignment

        skew = mu3/mu2**(3/2) 
        kurt = (mu4/mu2**2 - 3)
        return mu,mu2,skew,kurt
