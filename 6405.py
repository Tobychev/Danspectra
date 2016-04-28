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

s1 = spc.SpectraFactory("data/6405_aS1")
s1.frame_row_cut([0,799])
s1.contrast_cut(50)
s1.set_continua("segments")
qu1 = s1.make_spectra()
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1]

qu1H2O  = spc.splineline(wH2O,  cH2O,  qu1.meta,"H2O ") # Define lines for qu1 series...
qu1FeI  = spc.splineline(wFeI,  cFeI,  qu1.meta,"FeI ")
qu1SiFe = spc.splineline(wSiFe, cSiFe, qu1.meta,"SiFe ")
qu1Myst = spc.splineline(wMyst, cMyst, qu1.meta,"Myst ")
qu1CN   = spc.splineline(wCN,   cCN,   qu1.meta,"CN ")
qu1lines = [qu1H2O,qu1FeI,qu1SiFe,qu1Myst,qu1CN] # ... and save in this list for qu1 

qu1lims = {}
qu1lims["ewlim"]   = (-0.5 , 2.5 )
qu1lims["vellim"]  = (-10  , 10  )
qu1lims["rellim"]  = ( 0.2, 1.1)
qu1lims["fw13lim"] = (-0.1 , 1.4 )
qu1lims["fwhmlim"] = (-0.1 , 1.4 )
qu1lims["fw23lim"] = (-0.1 , 1.4 )
qu1lims["as13lim"] = (-0.7 , 0.7 )
qu1lims["as12lim"] = (-0.8 , 0.7 )
qu1lims["as23lim"] = (-1.1 , 0.6 )

b1 = spc.SpectraFactory("data/6405_bS1")
b1.frame_row_cut([0,799])
b1.contrast_cut(50)
b1.set_continua("segments")
qu2 = b1.make_spectra()
qu2con = qu2.meta.cont[0]*qu2.lmbd.mean() + qu2.meta.cont[1]

qu2H2O  = spc.splineline(wH2O,  cH2O,  qu2.meta,"H2O ") # Define lines for qu2 series...
qu2FeI  = spc.splineline(wFeI,  cFeI,  qu2.meta,"FeI ")
qu2SiFe = spc.splineline(wSiFe, cSiFe, qu2.meta,"SiFe ")
qu2Myst = spc.splineline(wMyst, cMyst, qu2.meta,"Myst ")
qu2CN   = spc.splineline(wCN,   cCN,   qu2.meta,"CN ")
qu2lines = [qu2H2O,qu2FeI,qu2SiFe,qu2Myst,qu2CN] # ... and save in this list for as2

qu2lims = {}
qu2lims["ewlim"]   = (-0.5 , 2.5 )
qu2lims["vellim"]  = (-10  , 10  )
qu2lims["rellim"]  = ( 0.1 , 1.3 )
qu2lims["fw13lim"] = (-0.2 , 1.2 )
qu2lims["fwhmlim"] = (-0.1 , 1.2 )
qu2lims["fw23lim"] = (-0.2 , 1.2 )
qu2lims["as13lim"] = (-0.7 , 0.7 )
qu2lims["as12lim"] = (-0.8 , 0.7 )
qu2lims["as23lim"] = (-1.1 , 0.6 )

#Hand placed limits from mean(axis=0)
wH2O  = [444, 465];   cH2O  = 640.8682
wFeI  = [505, 571];   cFeI  = 640.8026
wSiFe = [614, 647];   cSiFe = 640.7291
wMyst = [797, 853];   cMyst = 640.5763
wCN   = [1139,1178];   cCN  = 640.3127 

s2 = spc.SpectraFactory("data/6405_aS2",framerows=774,framecols=1446)
s2.frame_row_cut([0,743])
s2.frame_col_cut([0])
s2.contrast_cut(85)
s2.set_continua("segments")
spt = s2.make_spectra()
sptcon = spt.meta.cont[0]*spt.lmbd.mean() + spt.meta.cont[1]

sptH2O  = spc.splineline(wH2O,  cH2O,  spt.meta,"H2O ") # Define lines for spt series...
sptFeI  = spc.splineline(wFeI,  cFeI,  spt.meta,"FeI ")
sptSiFe = spc.splineline(wSiFe, cSiFe, spt.meta,"SiFe ")
sptMyst = spc.splineline(wMyst, cMyst, spt.meta,"Myst ")
sptlines = [sptH2O,sptFeI,sptSiFe,sptMyst] # ... and save in this list for spt
#OBS! H2O and CNq cause crash, removed

sptlims = {}
sptlims["ewlim"]   = ( 0.3 , 2.1 )
sptlims["vellim"]  = (-8   , 8   )
sptlims["rellim"]  = ( 0.2, 1.01)
sptlims["fw13lim"] = (-0.2 , 2   )
sptlims["fwhmlim"] = (-0.3 , 3   )
sptlims["fw23lim"] = (-0.4 , 3   )
sptlims["as13lim"] = (-1.5 , 0.74)
sptlims["as12lim"] = (-2.0 , 0.74)
sptlims["as23lim"] = (-2.0 , 0.84)

um = 0.38978
wl = 0.64927
pn = 0.83791

xspotlims = (640.38104838709671, 640.69435483870961)
yspotlims = (0.86458333333333337, 1.0437500000000002)


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
