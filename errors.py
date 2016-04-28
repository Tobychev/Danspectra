import numpy as np
import scipy.interpolate as intr
# For spline measurement
s_bot,s_cnt,s_fwhm,s_as12,s_fw13,s_as13,s_fw23,s_as23,s_err,s_ew,s_cont = np.arange(0,11)

def percentiles(measure):
    return ( np.percentile(measure,2.274),
             np.percentile(measure,15.87),
             measure.mean(),
             np.percentile(measure,84.13),
             np.percentile(measure,97.725))

def err_spline_mes(measure):
    sigm = np.zeros((11,5))

    sigm[s_bot,:]  = percentiles(measure[:,s_bot])
    sigm[s_bot,[0,1,3,4]]  = (sigm[s_bot,[0,1,3,4]]-sigm[s_bot,2])/sigm[s_bot,2] # Relativa fel 
    sigm[s_ew,:]   = percentiles(measure[:,s_ew])
    sigm[s_ew, [0,1,3,4]]  = (sigm[s_ew, [0,1,3,4]]-sigm[s_ew,2] )/sigm[s_ew,2]  # Relativa fel 
    sigm[s_cnt,:]  = percentiles(measure[:,s_cnt])
    sigm[s_cnt,[0,1,3,4]]  = (sigm[s_cnt,[0,1,3,4]]-sigm[s_cnt,2])/sigm[s_cnt,2] # Relativa fel 
    sigm[s_fwhm,:] = percentiles(measure[:,s_fwhm])
    sigm[s_fwhm,[0,1,3,4]] = sigm[s_fwhm,[0,1,3,4]] - sigm[s_fwhm,2] #Absolut fel
    sigm[s_as12,:] = percentiles(measure[:,s_as12])
    sigm[s_as12,[0,1,3,4]] = sigm[s_as12,[0,1,3,4]] - sigm[s_as12,2] #Absolut fel
    sigm[s_fw13,:] = percentiles(measure[:,s_fw13])
    sigm[s_fw13,[0,1,3,4]] = sigm[s_fw13,[0,1,3,4]] - sigm[s_fw13,2] #Absolut fel
    sigm[s_as13,:] = percentiles(measure[:,s_as13])
    sigm[s_as13,[0,1,3,4]] = sigm[s_as13,[0,1,3,4]] - sigm[s_as13,2] #Absolut fel
    sigm[s_fw23,:] = percentiles(measure[:,s_fw23])
    sigm[s_fw23,[0,1,3,4]] = sigm[s_fw23,[0,1,3,4]] - sigm[s_fw23,2] #Absolut fel
    sigm[s_as23,:] = percentiles(measure[:,s_as23])
    sigm[s_as23,[0,1,3,4]] = sigm[s_as23,[0,1,3,4]] - sigm[s_as23,2] #Absolut fel
    sigm[s_err,:] = percentiles(measure[:,s_err])
    sigm[s_cont,:] = percentiles(measure[:,s_cont])
    sigm[s_cont,[0,1,3,4]] = sigm[s_cont,[0,1,3,4]] - sigm[s_cont,2] #Absolut fel

    return sigm

def scale_spline_err(line,measure,errs):
    errs[s_bot,2] = measure[:,s_bot].mean()
    errs[s_bot,[0,1,3,4]] = errs[s_bot,2]*errs[s_bot,[0,1,3,4]]/measure[:,s_cont].mean()    
    errs[s_cnt,2] = measure[:,s_cnt].mean()
    errs[s_cnt,[0,1,3,4]] = errs[s_cnt,2]*errs[s_cnt,[0,1,3,4]]
    errs[s_ew,2] = measure[:,s_ew].mean()
    errs[s_ew,[0,1,3,4]] = errs[s_ew,2]*errs[s_ew,[0,1,3,4]]

    return errs

def error_inerpolate(vals,errors,quantity):
    return (intr.interp1d(vals,errors[:,quantity,0]),
            intr.interp1d(vals,errors[:,quantity,1]),
            intr.interp1d(vals,errors[:,quantity,3]),
            intr.interp1d(vals,errors[:,quantity,4]))

def make_intr_errs(errors,botrange):
    intrErrs = {}
    intrErrs["bot"]  = error_inerpolate(botrange,errors,s_bot)
    intrErrs["cnt"]  = error_inerpolate(botrange,errors,s_cnt)
    intrErrs["fwhm"] = error_inerpolate(botrange,errors,s_fwhm)
    intrErrs["fw13"] = error_inerpolate(botrange,errors,s_fw13)
    intrErrs["fw23"] = error_inerpolate(botrange,errors,s_fw23)
    intrErrs["as13"] = error_inerpolate(botrange,errors,s_as13)
    intrErrs["as12"] = error_inerpolate(botrange,errors,s_as12)
    intrErrs["as23"] = error_inerpolate(botrange,errors,s_as23)
    intrErrs["ew"]   = error_inerpolate(botrange,errors,s_ew)
    intrErrs["err"]   = error_inerpolate(botrange,errors,s_ew)
    return intrErrs


def scale_intr_spline_err(measure,line,intr_err):
    errs = np.zeros((10,4))
    bot = measure[:,s_bot ].mean()
    cnt = measure[:,s_cnt ].mean()
    ew  = measure[:,s_ew  ].mean()
    con = measure[:,s_cont].mean()

    errs[s_bot,:] = np.array([ intr_err["bot"][0](bot),
                               intr_err["bot"][1](bot),
                               intr_err["bot"][2](bot),
                               intr_err["bot"][3](bot)])
    errs[s_bot,:] = bot*errs[s_bot,:]/con

    errs[s_cnt,:] = np.array([ intr_err["cnt"][0](bot),
                               intr_err["cnt"][1](bot),
                               intr_err["cnt"][2](bot),
                               intr_err["cnt"][3](bot)])
    errs[s_cnt,:] = cnt*errs[s_cnt,:]

    errs[s_fwhm,:] = np.array([ intr_err["fwhm"][0](bot),
                               intr_err["fwhm"][1](bot),
                               intr_err["fwhm"][2](bot),
                               intr_err["fwhm"][3](bot)])
    errs[s_fwhm,:] = errs[s_fwhm,:]/line.width

    errs[s_as12,:] = np.array([ intr_err["as12"][0](bot),
                               intr_err["as12"][1](bot),
                               intr_err["as12"][2](bot),
                               intr_err["as12"][3](bot)])
    errs[s_as12,:] = errs[s_as12,:]/line.width

    errs[s_fw13,:] = np.array([ intr_err["fw13"][0](bot),
                               intr_err["fw13"][1](bot),
                               intr_err["fw13"][2](bot),
                               intr_err["fw13"][3](bot)])
    errs[s_fw13,:] = errs[s_fw13,:]/line.width

    errs[s_as13,:] = np.array([ intr_err["as13"][0](bot),
                               intr_err["as13"][1](bot),
                               intr_err["as13"][2](bot),
                               intr_err["as13"][3](bot)])
    errs[s_as13,:] = errs[s_as13,:]/line.width

    errs[s_fw23,:] = np.array([ intr_err["fw23"][0](bot),
                               intr_err["fw23"][1](bot),
                               intr_err["fw23"][2](bot),
                               intr_err["fw23"][3](bot)])
    errs[s_fw23,:] = errs[s_fw23,:]/line.width

    errs[s_as23,:] = np.array([ intr_err["as23"][0](bot),
                               intr_err["as23"][1](bot),
                               intr_err["as23"][2](bot),
                               intr_err["as23"][3](bot)])
    errs[s_as23,:] = errs[s_as23,:]/line.width

#   errs[s_err,:] = np.array([ intr_err["err"][0](bot),
#                              intr_err["err"][1](bot),
#                              intr_err["err"][2](bot),
#                              intr_err["err"][3](bot)])

    errs[s_ew,:] = np.array([ intr_err["ew"][0](bot),
                               intr_err["ew"][1](bot),
                               intr_err["ew"][2](bot),
                               intr_err["ew"][3](bot)])
    errs[s_ew,:] = ew*errs[s_ew,:]
    return errs
