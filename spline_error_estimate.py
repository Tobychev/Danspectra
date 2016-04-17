import matplotlib.pyplot as pl
import scipy.integrate as si
import scipy.stats as st
import spectra as spc
import numpy as np
import errors as er
import copy 

from imp import reload
reload(spc)

def mutant(x,line):
    x[:,wFlat[0]:wFlat[1]] = x[:,wFlat[0]:wFlat[1]]*line
    return x

def make_err_range(spline_measurement):
    sigm = np.zeros((11,5))
    for i in range(0,11):
        # The percentile ranges are one and two sigmas, From mathematics handbook
        sigm[i,:]= (
            np.percentile(spline_measurement[:,i],2.28),
            np.percentile(spline_measurement[:,i],15.87),
            spline_measurement[:,i].mean(),
            np.percentile(spline_measurement[:,i],84.13),
            np.percentile(spline_measurement[:,i],97.72))
        m2s,m1s,mu,p1s,p2s = sigm[i,:]
        print("-2s: {: 8.6f}, -1s: {: 8.6f}, mu: {: 8.6f}, +1s: {: 8.6f}, +2s: {: 8.6f}".format( (m2s-mu)/m2s,
                                                                                                 (m1s-mu)/m1s,
                                                                                                 mu,
                                                                                                 (p1s-mu)/p1s,
                                                                                                 (p2s-mu)/p2s ))

    return sigm 

def relative_errors(sigm):
    mu = sigm[:,2].reshape(-1,1).repeat(5,axis=1);
    sigm = (sigm - mu)/sigm
    sigm[:,2] = mu[:,2]
    return sigm

def kde(idx):
    rt =  st.gaussian_kde(mesFlat[:,idx])
    x  = np.linspace(mesFlat[:,idx].min(),mesFlat[:,idx].max(),121)
    pl.plot(x,rt(x))
    pl.show()

def printerr(err):
    for row in err:
        print("-2s: {: 8.6f}, -1s: {: 8.6f}, mu: {: 8.6f}, +1s: {: 8.6f}, +2s: {: 8.6f}".format(row[0],row[1],row[2],row[3],row[4]))

lineerror = "SiFe"

sf = spc.SpectraFactory("data/6405_aS1")
sf.frame_row_cut([0,799])
sf.contrast_cut(50)
sf.set_continua("segments")

sp = sf.make_spectra()

# Hand picked flat area
wFlat = [1179,1239]; cFlat = sp.lmbd[wFlat[0]:wFlat[1]].mean()


norm = st.norm.pdf(np.linspace(-6,6,60))
simFeI  = 1 - 0.62*norm/norm.max() # Line values eye estimated from atlas
simSiFe = 1 - 0.19*norm/norm.max() # Line values eye estimated from atlas
simMyst = 1 - 0.09*norm/norm.max() # Line values eye estimated from atlas

if lineerror == "SiFe":
    my = copy.deepcopy(sp) 
    my.modify(lambda x: mutant(x,simSiFe))
    ewtrue = si.simps(simSiFe-1,x=sp.lmbd[wFlat[0]:wFlat[1]])*1e3 

    Flat = spc.splineline(wFlat,cFlat,my.meta)
    mesFlat =  Flat.measure(my)
    err = er.err_spline_mes(mesFlat)
    print("Raw intervals")
    printerr(err)
    np.save("Estimate_errSiFe",err)

if lineerror == "Myst":
    my = copy.deepcopy(sp) 
    my.modify(lambda x: mutant(x,simMyst))
    ewtrue = si.simps(simMyst-1,x=sp.lmbd[wFlat[0]:wFlat[1]])*1e3 

    Flat = spc.splineline(wFlat,cFlat,my.meta)
    mesFlat =  Flat.measure(my)
    err = make_err_range(mesFlat)
    errMyst = er.err_spline_mes(mesFlat)

    np.save("Estimate_errMyst",errMyst)

if lineerror == "FeI":
    my = copy.deepcopy(sp) 
    my.modify(lambda x: mutant(x,simFeI))
    ewtrue = si.simps(simFeI-1,x=sp.lmbd[wFlat[0]:wFlat[1]])*1e3 

    Flat = spc.splineline(wFlat,cFlat,my.meta)
    mesFlat =  Flat.measure(my,smallstep=2e-9,numsmallstep=1e4)
    err = make_err_range(mesFlat)
    errFeI = er.err_spline_mes(mesFlat)

    np.save("Estimate_errFeI",errFeI)
