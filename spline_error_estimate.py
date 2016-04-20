import matplotlib.pyplot as pl
import scipy.integrate as si
import scipy.stats as st
import spectra as spc
import numpy as np
import errors as er
import pickle as pic
import copy 

# For spline measurement
s_bot,s_cnt,s_fwhm,s_as12,s_fw13,s_as13,s_fw23,s_as23,s_err,s_ew,s_cont = np.arange(0,11)

def mutant(x,line):
    x[:,wFlat[0]:wFlat[1]] = x[:,wFlat[0]:wFlat[1]]*line
    return x

def printerr(err):
    for row in err:
        print("-2s: {: 8.6f}, -1s: {: 8.6f}, mu: {: 8.6f}, +1s: {: 8.6f}, +2s: {: 8.6f}".format(row[0],row[1],row[2],row[3],row[4]))

sf = spc.SpectraFactory("data/6405_aS1")
sf.frame_row_cut([0,799])
sf.contrast_cut(50)
sf.set_continua("segments")

sp = sf.make_spectra()

# Hand picked flat area
wFlat = [1179,1239]; cFlat = sp.lmbd[wFlat[0]:wFlat[1]].mean()
botrange = np.linspace(0.04,0.9,22)
errors = []

#Chosen as it makes errors large
var = 1.5

for bot in botrange:
    norm = st.norm.pdf(np.linspace(-6,6,60),scale=var); norm = norm/norm.max()
    line = 1-bot*norm
    my = copy.deepcopy(sp) 
    my.modify(lambda x: mutant(x,line))
    print("Estimating errors for {}".format(bot))
    Flat = spc.splineline(wFlat,cFlat,my.meta)
    mesFlat =  Flat.measure(my)
    errors.append( er.err_spline_mes(mesFlat) )

errors = np.array(errors)

np.savez("SplineError.estimate",errors,botrange)
