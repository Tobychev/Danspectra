import spectra as spc

FeI  = {"El":0.9146,"gf":-4.318,"lam":640.0316,"dep":0.533,"name":"Fe I"}
Myst = {"El":    -1,"gf":     0,"lam":640.5763,"dep":0.058,"name":"Myst"}
SiI  = {"El":5.8709,"gf":-1.393,"lam":640.7291,"dep":0.155,"name":"Si I"}
FeI2 = {"El":3.6864,"gf":-1.018,"lam":640.8017,"dep":0.562,"name":"Fe I"}
SiI2 = {"El":5.9841,"gf":-1.554,"lam":640.8670,"dep":0.098,"name":"Si I"}

#Hand placed limits from mean(axis=0)
wFeI  = [1359,1388]
wMyst = [652, 721];
wSiI  = [478, 518];
wFeI2 = [377, 438];
wSiI2 = [310, 342];

sf_qu1 = spc.SpectraFactory("data/6405_aS1",framerows=800,framecols=1472)
sf_qu1.frame_col_cut([0,1471])
sf_qu1.frame_row_cut([0]+list(range(668,677))+[799])
sf_qu1.contrast_cut(50)
sf_qu1.set_continua("segments")

qu1 = sf_qu1.make_spectra()
qu1m = qu1[:,:].mean(axis=0) 
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1]

qu1FeI  = spc.splineline(wFeI,  FeI,  qu1.meta);qu1FeI.recenter(qu1m)
qu1Myst = spc.splineline(wMyst, Myst, qu1.meta);qu1Myst.recenter(qu1m)
qu1SiI  = spc.splineline(wSiI,  SiI,  qu1.meta);qu1SiI.recenter(qu1m)
qu1FeI2 = spc.splineline(wFeI2, FeI2, qu1.meta);qu1FeI2.recenter(qu1m)
qu1SiI2 = spc.splineline(wSiI2, SiI2, qu1.meta);qu1SiI2.recenter(qu1m)
qu1lines = [qu1FeI,qu1Myst,qu1SiI,qu1FeI2,qu1SiI2] # ... and save in this list for qu1 

sf_qu2 = spc.SpectraFactory("data/6405_bS1",framerows=756,framecols=1472)
sf_qu2.frame_col_cut([0,1,1471])
sf_qu2.frame_row_cut([0]+list(range(663,670))+[799])
sf_qu2.contrast_cut(50)
sf_qu2.set_continua("segments")
qu2 = sf_qu2.make_spectra()
qu2con = qu2.meta.cont[0]*qu2.lmbd.mean() + qu2.meta.cont[1]

qu2FeI  = spc.splineline(wFeI,  FeI,  qu2.meta)
qu2Myst = spc.splineline(wMyst, Myst, qu2.meta)
qu2SiI  = spc.splineline(wSiI,  SiI,  qu2.meta)
qu2FeI2 = spc.splineline(wFeI2, FeI2, qu2.meta)
qu2SiI2 = spc.splineline(wSiI2, SiI2, qu2.meta)
qu2lines = [qu2FeI,qu2Myst,qu2SiI,qu2FeI2,qu2SiI2] # ... and save in this list for qu2 


#Hand placed limits from mean(axis=0)
wMyst = [797, 850]
wSiI  = [611, 647]
wFeI2 = [505, 571]
wSiI2 = [442, 470]

sf_spt = spc.SpectraFactory("data/6405_aS2",framerows=774,framecols=1446)
sf_spt.frame_col_cut([0,1445])
sf_spt.frame_row_cut([0]+list(range(659,668))+[743])
sf_spt.contrast_cut(85)
sf_spt.set_continua("segments")
spt = sf_spt.make_spectra()
sptcon = spt.meta.cont[0]*spt.lmbd.mean() + spt.meta.cont[1]

sptMyst = spc.splineline(wMyst, Myst, spt.meta)
sptSiI  = spc.splineline(wSiI,  SiI,  spt.meta)
sptFeI2 = spc.splineline(wFeI2, FeI2, spt.meta)
sptSiI2 = spc.splineline(wSiI2, SiI2, spt.meta)
sptlines = [sptMyst,sptSiI,sptFeI2,sptSiI2] # ... and save in this list for spt 

um = 0.38978
wl = 0.64927
pn = 0.83791

xspotlims = (640.38104838709671, 640.69435483870961)
yspotlims = (0.86458333333333337, 1.0437500000000002)

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
