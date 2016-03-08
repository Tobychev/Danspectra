#encoding: utf-8
from wipfile import *

line = myst
#dlam = np.diff(lam[slice(line.idx[0]-1,line.idx[-1]+1)]).reshape((-1,1))*np.ones(798)
#EW   = ((frm.data[:,line.idx]-1)*dlam.T).sum(axis=1)*1e3 ## MiliÅngström
s6405_t5p.normalize()

EW =  line.measure_EW(s6405_t5p)
