import numpy as np

def percentiles(measure):
    return ( np.percentile(measure,2.274),
             np.percentile(measure,15.87),
             measure.mean(),
             np.percentile(measure,84.13),
             np.percentile(measure,97.725))

def err_spline_mes(measure):
    # For spline measurement
    bot,cnt,fwhm,as12,fw13,as13,fw23,as23,err,ew,cont = np.arange(0,11)

    sigm = np.zeros((11,5))

    sigm[bot,:]  = percentiles(measure[:,bot]);  sigm[bot,[0,1,3,4]]  = (sigm[bot,[0,1,3,4]]-sigm[bot,2])/sigm[bot,2] # Relativa fel 
    sigm[ew,:]   = percentiles(measure[:,ew]);   sigm[ew, [0,1,3,4]]  = (sigm[ew, [0,1,3,4]]-sigm[ew,2] )/sigm[ew,2]  # Relativa fel 
    sigm[cnt,:]  = percentiles(measure[:,cnt]);  sigm[cnt,[0,1,3,4]]  = (sigm[cnt,[0,1,3,4]]-sigm[cnt,2])/sigm[cnt,2] # Relativa fel 
    sigm[fwhm,:] = percentiles(measure[:,fwhm]); sigm[fwhm,[0,1,3,4]] = sigm[fwhm,[0,1,3,4]] - sigm[fwhm,2] #Absolut fel
    sigm[as12,:] = percentiles(measure[:,as12]); sigm[as12,[0,1,3,4]] = sigm[as12,[0,1,3,4]] - sigm[as12,2] #Absolut fel
    sigm[fw13,:] = percentiles(measure[:,fw13]); sigm[fw13,[0,1,3,4]] = sigm[fw13,[0,1,3,4]] - sigm[fw13,2] #Absolut fel
    sigm[as13,:] = percentiles(measure[:,as13]); sigm[as13,[0,1,3,4]] = sigm[as13,[0,1,3,4]] - sigm[as13,2] #Absolut fel
    sigm[fw23,:] = percentiles(measure[:,fw23]); sigm[fw23,[0,1,3,4]] = sigm[fw23,[0,1,3,4]] - sigm[fw23,2] #Absolut fel
    sigm[as23,:] = percentiles(measure[:,as23]); sigm[as23,[0,1,3,4]] = sigm[as23,[0,1,3,4]] - sigm[as23,2] #Absolut fel
    sigm[err,:] = percentiles(measure[:,err]);
    sigm[cont,:] = percentiles(measure[:,cont]); sigm[cont,[0,1,3,4]] = sigm[cont,[0,1,3,4]] - sigm[cont,2] #Absolut fel

    return sigm

def scale_spline_err(line,measure,errs):
    # For spline measurement
    bot,cnt,fwhm,as12,fw13,as13,fw23,as23,err,ew,cont = np.arange(0,11)

    errs[bot,2] = measure[:,bot].mean()
    errs[bot,[0,1,3,4]] = errs[bot,2]*errs[bot,[0,1,3,4]]/measure[:,cont].mean()    
    errs[cnt,2] = measure[:,cnt].mean()
    errs[cnt,[0,1,3,4]] = errs[cnt,2]*errs[cnt,[0,1,3,4]]
    errs[ew,2] = measure[:,ew].mean()
    errs[ew,[0,1,3,4]] = errs[ew,2]*errs[ew,[0,1,3,4]]

    return errs

    
