from wipfile import *
import statsmodels.distributions.empirical_distribution as empd


ecdf  = empd.ECDF(frm.data[331,line.idx])
ecdf2 = empd.ECDF(frm.data[31,line.idx])
ecdf3 = empd.ECDF(frm.data[421,line.idx])

if False:
    pl.step(ecdf.x, ecdf.y)
    pl.step(ecdf2.x, ecdf2.y)
    pl.step(ecdf3.x, ecdf3.y)
    pl.show()

data = frm.data
line = lins[0]
dpdf = (1-frm.data[331,line.idx])/np.sum(1-frm.data[331,line.idx])

x    = lam[line.idx]

mu = np.sum(dpdf*x) # First moment, mean
mu2 = np.sum(dpdf*(x-mu)**2)
mu3 = np.sum(dpdf*(x-mu)**3)
mu4 = np.sum(dpdf*(x-mu)**4)

print(mu,mu2,mu3,mu4)
print("Skew: {}".format(mu3/mu2**(3/2)))
print("Kurt: {}".format(mu4/mu2**2  -3))

dpdf = ((1-frm.data[:,line.idx]).T/np.sum(1-frm.data[:,line.idx],axis=1)).T
mu = np.sum(dpdf*x,axis=1).reshape(-1,1) # First moment, mean
mu2 = np.sum(dpdf*((x-mu)**2),axis=1)
mu3 = np.sum(dpdf*(x-mu)**3,axis=1)
mu4 = np.sum(dpdf*(x-mu)**4,axis=1)

print(mu.mean(),mu2.mean(),mu3.mean(),mu4.mean())
pl.hist(mu3/mu2**(3/2),21,)
pl.title("Skewness")
pl.show()

pl.hist(mu4/mu2**2 - 3,21)
pl.title("Kurtosis")
pl.show()


#print("Skew: {}".format(mu3/mu2**(3/2)))
#print("Kurt: {}".format(mu4/mu2**2  -3))
