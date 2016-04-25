import matplotlib.pyplot as pl
import visualize as vis
import astropy.io.fits as fts
import binner as bn
import spectra as spc
import pickle as pic
import numpy as np
import errors as er

#Hand placed limits from mean(axis=0)
wH2O  = [311, 343];   cH2O  = 640.8682
wFeI  = [378, 439];   cFeI  = 640.8026
wSiFe = [479, 519];   cSiFe = 640.7291
wMyst = [646, 722];   cMyst = 640.5763
wCN  = [1006,1047];   cCN   = 640.3127 

lims = {}
lims["ewlim"]   = ( -0.5 , 2.5  )
lims["vellim"]  = (-10 , 10  )
lims["rellim"]  = ( 0.1 , 1.3  )

lims["fw13lim"] = ( -0.2 , 1.2  )
lims["fwhmlim"] = ( -0.1 , 1.2  )
lims["fw23lim"] = (-0.2 , 1.2  )

lims["as13lim"] = (-0.7 , 0.7  )
lims["as12lim"] = (-0.8 , 0.7  )
lims["as23lim"] = (-1.1 , 0.6  )
as1lims = lims

s1 = spc.SpectraFactory("data/6405_aS1")
s1.frame_row_cut([0,799])
s1.contrast_cut(50)
s1.set_continua("segments")
as1 = s1.make_spectra()
as1con = as1.meta.cont[0]*as1.lmbd.mean() + as1.meta.cont[1]

as1H2O  = spc.splineline(wH2O,  cH2O,  as1.meta,"H2O ") # Define lines for as1 series...
as1FeI  = spc.splineline(wFeI,  cFeI,  as1.meta,"FeI ")
as1SiFe = spc.splineline(wSiFe, cSiFe, as1.meta,"SiFe ")
as1Myst = spc.splineline(wMyst, cMyst, as1.meta,"Myst ")
as1CN   = spc.splineline(wCN,   cCNq,  as1.meta,"CN ")
as1lines = [as1H2O,as1FeI,as1SiFe,as1Myst,as1CN] # ... and save in this list for as2

b1 = spc.SpectraFactory("data/6405_bS1")
b1.frame_row_cut([0,799])
b1.contrast_cut(50)
b1.set_continua("segments")
bs1 = b1.make_spectra()
bs1con = bs1.meta.cont[0]*bs1.lmbd.mean() + bs1.meta.cont[1]

bs1H2O  = spc.splineline(wH2O,  cH2O,  bs1.meta,"H2O ") # Define lines for bs1 series...
bs1FeI  = spc.splineline(wFeI,  cFeI,  bs1.meta,"FeI ")
bs1SiFe = spc.splineline(wSiFe, cSiFe, bs1.meta,"SiFe ")
bs1Myst = spc.splineline(wMyst, cMyst, bs1.meta,"Myst ")
bs1CN  = spc.splineline(wCN,  cCN,  bs1.meta,"CN ")
bs1lines = [bs1H2O,bs1FeI,bs1SiFe,bs1Myst,bs1CN] # ... and save in this list for as2

#Hand placed limits from mean(axis=0)
wH2O  = [442, 470];   cH2O  = 640.8682
wFeI  = [505, 571];   cFeI  = 640.8026
wSiFe = [611, 651];   cSiFe = 640.7291
wMyst = [798, 852];   cMyst = 640.5763
wCN   = [1139,1178];   cCN  = 640.3127 

s2 = spc.SpectraFactory("data/6405_aS2",framerows=774,framecols=1446)
s2.frame_row_cut([0,743])
s2.frame_col_cut([0])
s2.contrast_cut(85)
s2.set_continua("segments")
as2 = s2.make_spectra()
as2con = as2.meta.cont[0]*as2.lmbd.mean() + as2.meta.cont[1]

as2H2O  = spc.splineline(wH2O,  cH2O,  as2.meta,"H2O ") # Define lines for as2 series...
as2FeI  = spc.splineline(wFeI,  cFeI,  as2.meta,"FeI ")
as2SiFe = spc.splineline(wSiFe, cSiFe, as2.meta,"SiFe ")
as2Myst = spc.splineline(wMyst, cMyst, as2.meta,"Myst ")
as2CN  = spc.splineline(wCN,    cCN,   as2.meta,"CN ")

#OBS! H2O and CNq cause crash, removed
#as2lines = [as2H2O,as2FeI,as2SiFe,as2Myst,as2CN] # ... and save in this list for as2
as2lines = [as2FeI,as2SiFe,as2Myst] # ... and save in this list for as2


# Sunspot measurement
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

if False:
    # For spline measurement
    bot,cnt,fwhm,as12,fw13,as13,fw23,as23,err,ew,cont = np.arange(0,11)

    wFlat = [1179,1239];
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
#    mesMyst = Myst.measure(as1)
    print("Measuring Atmo H2O line")
#    mesH2O  = H2O.measure(as1) # Seems to depend on continuua, so probably not atmo line?
    print("Measuring Iron line")
#    mesFeI  = FeI.measure(as1)
    print("Measuring Si + Fe line")
    mesSiFe = SiFe.measure(as1)
    linmap = vis.spline_linemap(mesSiFe,SiFe)
    npz = np.load("SplineError.estimate.npz")
    errs,vals = npz["arr_0"],npz["arr_1"]
    intrErr = er.make_intr_errs(errs,vals) 
    errSiFe = er.scale_intr_spline_err(mesSiFe,SiFe,intrErr)
    vis.dan_errplot(linmap,errSiFe).show()
#    pl.plot(mesMyst[:,cont],mesMyst[:,bot],'o')
#    err = np.load("Estimate_errSiFe.npy")
#    errSiFe = er.scale_spline_err(SiFe,mesSiFe,err)
    
