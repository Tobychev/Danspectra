from wipfile import *
import statsmodels.distributions.empirical_distribution as empd

spec  = frm.data[331,:]
spec2 = frm.data[31,:]
spec3 = frm.data[421,:]
ecdf  = empd.ECDF(spec[line.idx])
ecdf2 = empd.ECDF(spec2[line.idx])
ecdf3 = empd.ECDF(spec3[line.idx])

if False:
    pl.step(ecdf.x, ecdf.y)
    pl.step(ecdf2.x, ecdf2.y)
    pl.step(ecdf3.x, ecdf3.y)
    pl.show()

line = SiFe
dpdf1 = (1-spec[line.idx])/np.sum(1-spec[line.idx])
dpdf2 = (1-spec2[line.idx])/np.sum(1-spec2[line.idx])
dpdf3 = (1-spec3[line.idx])/np.sum(1-spec3[line.idx])


dpdf = dpdf1 #Choose which example to use
x    = lam[line.idx]

mu = np.sum(dpdf*x) # First moment, mean
mu2 = np.sum(dpdf*(x-mu)**2)
mu3 = np.sum(dpdf*(x-mu)**3)
mu4 = np.sum(dpdf*(x-mu)**4)

print mu,mu2,mu3,mu4
print "Skew: {}".format(mu3/mu2**(3/2))
print "Kurt: {}".format(mu4/mu2**2  -3)
