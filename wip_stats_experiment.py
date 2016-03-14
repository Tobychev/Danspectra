from wipfile import *
import statsmodels.distributions.empirical_distribution as empd

s6405_t5p.normalize()
frm = s6405_t5p.frames[0]
#frm.data = frm.data/frm.cont.norm()
line = FeI

if False:
    ecdf  = empd.ECDF(frm.data[331,line.idx])
    ecdf2 = empd.ECDF(frm.data[31,line.idx])
    ecdf3 = empd.ECDF(frm.data[421,line.idx])
    pl.step(ecdf.x, ecdf.y)
    pl.step(ecdf2.x, ecdf2.y)
    pl.step(ecdf3.x, ecdf3.y)
    pl.show()

if False:
    data = frm.data
    line = myst
    dpdf = (1-frm.data[34,line.idx])/np.sum(1-frm.data[35,line.idx])

    x    = lam[line.idx]

    mu = np.sum(dpdf*x) # First moment, mean
    mu2 = np.sum(dpdf*(x-mu)**2)
    mu3 = np.sum(dpdf*(x-mu)**3)
    mu4 = np.sum(dpdf*(x-mu)**4)

    print(mu,mu2,mu3,mu4)
    print("Vars: {}".format(mu2))
    print("Skew: {}".format(mu3/mu2**(3/2)))
    print("Kurt: {}".format(mu4/mu2**2  -3))

if False:
    frame = frm
    x    = frame.group.lmbd[line.idx]
    for row,slask in enumerate(frame.data[:,0]):
        dpdf = (1-frame.data[row,line.idx]); 
        if any(dpdf < 0.0):
           dpdf = dpdf -dpdf.min()
        dpdf = dpdf/dpdf.sum()
        mu   = np.sum(dpdf*x) # Reshaping enables broadcasting
        mu2  = np.sum(dpdf*(x-mu)**2)
        if mu2 < 0 :
            print("Row:",row)
            print("Central diff:",(x-mu)**2)
            print("Prod:",dpdf*(x-mu)**2)
            print("Sum: ",np.sum(dpdf*(x-mu)**2))

        mu3  = np.sum(dpdf*(x-mu)**3)
        mu4  = np.sum(dpdf*(x-mu)**4)
        mu   = mu.reshape(-1) # Undoing reshape to allow assignment

        skew = mu3/mu2**(3/2)
#        print("lam: ",x.mean())
#        print("ePdf: ",dpdf.max(),dpdf.min(),dpdf.mean())
#        print("Moments: ", mu.mean(),mu2.mean(),mu3.mean(),mu4.mean())
#        print("Skewness: ",skew)

x    = frm.group.lmbd[line.idx]
dpdf = (1-frm.data[:,line.idx]/frm.data[:,line.idx].max(axis=1).reshape(-1,1))
dpdf =dpdf/dpdf.sum(axis=1).reshape(-1,1)

#dpdf = ((1-frm.data[:,line.idx]).T/np.sum(1-frm.data[:,line.idx],axis=1)).T
mu = np.sum(dpdf*x,axis=1).reshape(-1,1) # First moment, mean
mu2 = np.sum(dpdf*((x-mu)**2),axis=1)
mu3 = np.sum(dpdf*(x-mu)**3,axis=1)
mu4 = np.sum(dpdf*(x-mu)**4,axis=1)
skew = mu3/mu2**(3/2)

pl.hist(mu3/mu2**(3/2),21)
pl.title("Skewness")
pl.show()

pl.hist(mu4/mu2**2 - 3,21)
pl.title("Kurtosis")
pl.show()


#print("Skew: {}".format(mu3/mu2**(3/2)))
#print("Kurt: {}".format(mu4/mu2**2  -3))
