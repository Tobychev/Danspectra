from wipfile import *
import astropy.stats as ast
import numpy.polynomial.polynomial as pol
s6405_t5p.normalize()

data = s6405_t5p.frames[0].data[:,myst.idx]
for frm in s6405_t5p.frames[1:]:
    data = np.concatenate((data,frm.data[:,myst.idx]),axis=0)

vel = 0; bot = 1; con = 2; err = 3; ew  = 4; mn  = 5; var = 6; ske = 7; kur = 8
mesFeI   = FeI.measure(s6405_t5p)
mesSiFe  = SiFe.measure(s6405_t5p)
mesmyst  = myst.measure(s6405_t5p)

cuts = mesmyst[err] < np.percentile(mesmyst[err],89)
if False:    
    # Note the fancy histogram function, mine for reference
    counts, bins = ast.histogram(mesmyst[con][cuts],bins='blocks')

    data = data[np.where(cuts),:]

    sorting = np.digitize(mesmyst[con][cuts],bins,right=True)
    binned = []
    for i in np.unique(sorting)[:-1]:
    #    binned.append( DATA[sorting == i,:].mean(axis=0) )
        binned.append( np.mean(data[sorting == i,:],axis=0) )
        print(counts[i])
        pl.step(s6405_t5p.lmbd[myst.idx],binned[-1]);
        pl.show()


if True:
    quant = mesmyst[con].reshape(-1)
    try:             
        cuts, = np.where(cuts.reshape(-1))
    except AttributeError:
        cuts = np.ones(len(quant))
        
    datablock = s6405_t5p.frames[0].data[:,myst.idx]
    contblock  = s6405_t5p.frames[0].cont.val(myst.cent)
    for frm in s6405_t5p.frames[1:]:
        datablock = np.concatenate((datablock,frm.data[:,myst.idx]),axis=0)
        contblock = np.concatenate((contblock,frm.cont.val(myst.cent)),axis=0)

    datablock = datablock[cuts,:]
    contblock = contblock[cuts]
    print(contblock.shape,datablock.shape)

    counts, bins = ast.histogram(quant[cuts],bins='blocks')
    sorting = np.digitize(quant[cuts],bins[:-1]) # :-1 to get the correct number of buckets from digitize
    binned  =  np.mean(datablock[sorting == 1,:],axis=0)
    con     =  np.zeros(len(counts))
    con[0]  = np.mean( contblock[sorting == 1] )

    print( counts.shape, np.unique(sorting).shape )
    for i in np.unique(sorting)[1:]:
        print(i)
#        print(binned.shape,con.shape)
        print(np.mean(datablock[sorting == i,:],axis=0).shape,(contblock[sorting == i]).shape)
        binned = np.vstack( (binned,np.mean(datablock[sorting == i,:],axis=0) ))
        con[i-1]    = np.mean(contblock[sorting == i])


class binned_framegroup(object):

    def __init__(self,line,group,bin_quant,cuts=None):
        self.idx   = line.idx
        self.cent  = line.cent
        self.group = group
        self.data,self.cont = self.__bin_by_quant(bin_quant,cuts)
    
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
        sorting = np.digitize(quant[cuts],bins,right=True)
        binned  =  datablock[sorting == 0,:]
        con     =  contblock[sorting == 0]
        for i in np.unique(sorting)[1:]:
            binned = np.vstack( (binned,np.mean(datablock[sorting == i,:],axis=0) ))
            con    = np.vstack( (con,np.mean(contblock[sorting == i]) ))
        return binned,con

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
        fit   = pol.polyfit(self.group.lmbd[bottom],self.data[:,bottom].T,2)
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



#quant = mesmyst[con].reshape(-1)
#test  = binned_framegroup(myst,s6405_t5p,quant,mesmyst[err] < np.percentile(mesmyst[err],89))
#mes   = test.measure()
