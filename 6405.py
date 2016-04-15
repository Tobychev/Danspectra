import matplotlib.pyplot as pl
import visualize as vis
import astropy.io.fits as fts
import binner as bn
import spectra as spc
import numpy as np

s1 = spc.SpectraFactory("data/6405_aS1")
s1.frame_row_cut([0,799])
s1.contrast_cut(50)
s1.set_continua("segments")
as1 = s1.make_spectra()

s2 = spc.SpectraFactory("data/6405_aS2",framerows=774)
s2.frame_row_cut([0,743])
s2.contrast_cut(85)
s2.set_continua("segments")
as2 = s2.make_spectra()
as2con = as2.meta.cont[0]*as2.lmbd.mean() + as2.meta.cont[1]

b1 = spc.SpectraFactory("data/6405_bS1")
b1.frame_row_cut([0,799])
b1.contrast_cut(50)
b1.set_continua("segments")
bS1 = b1.make_spectra()

# Sunspot measurement
# Need to cut first column, its trash
if True:
    umbra    = 0.35
    wall     = (0.35,0.75)
    penumbra = (0.75,0.89)
    quiet    = (0.89)
    as2um    = s2.make_spectra_subset(as2,rowsubset=(as2con < umbra),desc="Umbra subset")
    as2wl    = s2.make_spectra_subset(as2,rowsubset=((as2con >= wall[0]) & (as2con <= wall[1])),desc="Umbra/penumbra wall")
    as2pn    = s2.make_spectra_subset(as2,rowsubset=((as2con > penumbra[0]) & (as2con < penumbra[1])),desc="Penumbra")
    as2qu    = s2.make_spectra_subset(as2,rowsubset=(as2con >= quiet),desc="Quiet sun")
    
    pl.plot(as2um.lmbd[1:], as2qu[:,1:].mean(axis=0));pl.show()
