import spectra as spc

FeI   = {"El":4.4733,"gf":-2.000,"lam":565.1468,"dep":0.121,"name":"Fe I"}
Myst  = {"El":    -1,"gf":     0,"lam":565.4501,"dep":0.125,"name":"Myst"}
SiI   = {"El":5.6136,"gf":-1.887,"lam":565.4919,"dep":0.362,"name":"Si I"}
VI    = {"El":1.0636,"gf":-1.000,"lam":565.7440,"dep":0.082,"name":"V I"}
ScII  = {"El":1.5070,"gf":-0.603,"lam":565.7896,"dep":0.491,"name":"Sc II"}
FeI2  = {"El":4.2844,"gf":-1.736,"lam":566.1344,"dep":0.257,"name":"Fe I"}


#Hand selected limits
wFeI  = [1216,1248]
wMyst = [873,937]; 
wSiI  = [842,867]; 
wVI   = [567,595]; 
wScII = [512,548]; 
wFeI2 = [136,165]; 

sf_qu1 = spc.SpectraFactory("data/5654_aS2",framerows=756,framecols=1480)
sf_qu1.frame_row_cut([0]+list(range(653,660))+[755])
sf_qu1.frame_col_cut([0,1479])
sf_qu1.contrast_cut(20)
sf_qu1.set_continua("segments")

qu1 = sf_qu1.make_spectra() # Do this before defining lines
qu1m = qu1[:,:].mean(axis=0) 
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1] # Define continua for qu1 series

qu1FeI  = spc.splineline(wFeI , FeI , qu1.meta);qu1FeI.recenter(qu1m)# Define lines for qu1 series...
qu1Myst = spc.splineline(wMyst, Myst, qu1.meta);qu1Myst.recenter(qu1m)
qu1SiI  = spc.splineline(wSiI , SiI , qu1.meta);qu1SiI.recenter(qu1m)
qu1VI   = spc.splineline(wVI  , VI  , qu1.meta);qu1VI.recenter(qu1m)
qu1ScII = spc.splineline(wScII, ScII, qu1.meta);qu1ScII.recenter(qu1m)
qu1FeI2 = spc.splineline(wFeI2, FeI2, qu1.meta);qu1FeI2.recenter(qu1m)
qu1lines = [qu1FeI,qu1Myst,qu1SiI,qu1VI,qu1ScII,qu1FeI2] # ... and save in this list for qu1


sf_qu2 = spc.SpectraFactory("data/5654_cS2",framerows=756,framecols=1480)
sf_qu2.frame_row_cut([0]+list(range(657,664))+[755])
sf_qu2.frame_col_cut([0,1479])
sf_qu2.contrast_cut(20)
sf_qu2.set_continua("segments")

qu2 = sf_qu2.make_spectra()
qu2m = qu2[:,:].mean(axis=0) 
qu2con = qu2.meta.cont[0]*qu2.lmbd.mean() + qu2.meta.cont[1] # Define continua for qu2 series

qu2FeI = spc.splineline(wFeI  ,  FeI , qu2.meta);qu2FeI.recenter(qu2m)# # Define lines for qu2 series...
qu2Myst = spc.splineline(wMyst,  Myst, qu2.meta);qu2Myst.recenter(qu2m)
qu2SiI  = spc.splineline(wSiI ,  SiI , qu2.meta);qu2SiI.recenter(qu2m)
qu2VI   = spc.splineline(wVI  ,  VI  , qu2.meta);qu2VI.recenter(qu2m)
qu2ScII = spc.splineline(wScII,  ScII, qu2.meta);qu2ScII.recenter(qu2m)
qu2FeI2  = spc.splineline(wFeI2, FeI2, qu2.meta);qu2FeI2.recenter(qu2m)
qu2lines = [qu2FeI,qu2Myst,qu2SiI,qu2VI,qu2ScII,qu2FeI2] # ... and save in this list for qu2


# Updating the limits to better match the meanspectra of this series
wFeI  = [1217,1253]; cFeI  = 565.1477
wMyst = [873,937];   cMyst = 565.4501
wSiI  = [844,871];   cSiI  = 565.4937
wVI   = [564,596];   cVI   = 565.7456
wScII = [512,555];   cScII = 565.7880
wFeI2 = [139,167];   cFeI2 = 566.1354

sf_spt = spc.SpectraFactory("data/5654_bS2",framerows=756,framecols=1480)
sf_spt.frame_row_cut([0]+list(range(653,662))+[755])
sf_spt.frame_col_cut([0,1479])
sf_spt.contrast_cut(50)
sf_spt.set_continua("segments")

spt = sf_spt.make_spectra()
sptm = qu1m # spt[:,:].mean(axis=0) 
sptcon = spt.meta.cont[0]*spt.lmbd.mean() + spt.meta.cont[1] # Define continua for spt series

sptFeI  = spc.splineline(wFeI , FeI , spt.meta);sptFeI.recenter(sptm)# Define lines for spt series... 
sptMyst = spc.splineline(wMyst, Myst, spt.meta);sptMyst.recenter(sptm)
sptSiI  = spc.splineline(wSiI , SiI , spt.meta);sptSiI.recenter(sptm)
sptVI   = spc.splineline(wVI  , VI  , spt.meta);sptVI.recenter(sptm)
sptScII = spc.splineline(wScII, ScII, spt.meta);sptScII.recenter(sptm)
sptFeI2 = spc.splineline(wFeI2, FeI2, spt.meta);sptFeI2.recenter(sptm)
sptlines = [sptFeI,sptMyst,sptSiI,sptVI,sptScII,sptFeI2] # ... and save in this list for spt

um = 0.37468
wl = 0.60141
pn = 0.85094

xspotlims = (565.21287720331078, 565.7816432702152)
yspotlims = (0.75923171274393061, 1.0552637808945176) 

qu1lims = {}
qu1lims["ewlim"]   = ( 0.6 , 1.4 )
qu1lims["vellim"]  = (-4.0 , 4.1 )
qu1lims["rellim"]  = ( 0.3 , 1.01)
qu1lims["fw13lim"] = ( 0.0 , 0.6)
qu1lims["fwhmlim"] = ( 0.0 , 0.7 )
qu1lims["fw23lim"] = ( -0.2 , 1.1)
qu1lims["as13lim"] = (-0.08, 0.1 )
qu1lims["as12lim"] = (-0.08, 0.1 )
qu1lims["as23lim"] = (-0.7, 0.2 )

qu2lims = {}
qu2lims["ewlim"]   = ( 0.6 , 1.4 )
qu2lims["vellim"]  = (-4.0 , 4.1 )
qu2lims["rellim"]  = ( 0.3 , 1.01)
qu2lims["fw13lim"] = ( -0.1 , 0.6)
qu2lims["fwhmlim"] = ( -0.1 , 0.8)
qu2lims["fw23lim"] = (- 0.1 , 1.1)
qu2lims["as13lim"] = (-0.08, 0.1 )
qu2lims["as12lim"] = (-0.08, 0.1 )
qu2lims["as23lim"] = (-0.2, 0.25 )

sptlims = {}
sptlims["ewlim"]   = ( -0.5 , 4.2 )
sptlims["vellim"]  = (-4.0 , 5 )
sptlims["rellim"]  = ( 0.3 , 1.01)
sptlims["fw13lim"] = ( -0.1 , 1)
sptlims["fwhmlim"] = ( -0.1 , 1.1 )
sptlims["fw23lim"] = ( -0.2 , 1.1)
sptlims["as13lim"] = (-0.3, 0.2 )
sptlims["as12lim"] = (-0.7, 0.25 )
sptlims["as23lim"] = (-0.8, 0.3 )


