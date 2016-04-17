import matplotlib.pyplot as pl
import visualize as vis
import astropy.io.fits as fts
import binner as bn
import spectra as spc
import numpy as np
import errors as er

s1 = spc.SpectraFactory("data/6405_aS1")
s1.frame_row_cut([0,799])
s1.contrast_cut(50)
s1.set_continua("segments")
as1 = s1.make_spectra()

s2 = spc.SpectraFactory("data/6405_aS2",framerows=774,framecols=1446)
s2.frame_row_cut([0,743])
s2.frame_col_cut([0])
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
if False:
    umbra    = 0.35
    wall     = (0.35,0.75)
    penumbra = (0.75,0.89)
    quiet    = (0.89)
    as2um    = s2.make_spectra_subset(as2,rowsubset=(as2con < umbra),desc="Umbra subset")
    as2wl    = s2.make_spectra_subset(as2,rowsubset=((as2con >= wall[0]) & (as2con <= wall[1])),desc="Umbra/penumbra wall")
    as2pn    = s2.make_spectra_subset(as2,rowsubset=((as2con > penumbra[0]) & (as2con < penumbra[1])),desc="Penumbra")
    as2qu    = s2.make_spectra_subset(as2,rowsubset=(as2con >= quiet),desc="Quiet sun")
    
    pl.plot(as2um.lmbd, as2um[:,:].mean(axis=0));pl.show()

if True:
    # For spline measurement
    bot,cnt,fwhm,as12,fw13,as13,fw23,as23,err,ew,cont = np.arange(0,11)

    wH2O  = [289, 342];   cH2O  = 640.86681253796701
    wFeI  = [376, 441];   cFeI  = 640.80263785584464
    wSiFe = [467, 519];   cSiFe = 640.72889305314993
    wMyst = [646, 722];   cMyst = 640.57558287431198 #Hand placed limits from mean(axis=0)
    wCNq  = [1006,1047];  cCNq  = 640.31162454661717 #Hand placed limits from mean(axis=0)

    Myst = spc.splineline(wMyst,cMyst,as1.meta)
    FeI  = spc.splineline(wFeI ,cFeI ,as1.meta)
    SiFe = spc.splineline(wSiFe,cSiFe,as1.meta)
    H2O  = spc.splineline(wH2O ,cH2O ,as1.meta)

    print("Measuring unknown line")
    mesMyst = Myst.measure(as1)
    print("Measuring Atmo H2O line")
    mesH2O  = H2O.measure(as1) # Seems to depend on continuua, so probably not atmo line?
    print("Measuring Iron line")
    mesFeI  = FeI.measure(as1)
    print("Measuring Si + Fe line")
    mesSiFe = SiFe.measure(as1)
    linmap = vis.spline_linemap(mesSiFe,SiFe)
    err = np.load("Estimate_errSiFe.npy")
    errSiFe = er.scale_spline_err(SiFe,mesSiFe,err)
    vis.dan_errplot(linmap,errSiFe).show()
#    pl.plot(mesMyst[:,cont],mesMyst[:,bot],'o')
    
