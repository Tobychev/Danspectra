import spectra as spc
FeI   = {"El":3.2368,"gf":-3.260,"lam":520.9884,"dep":0.121,"name":"Fe I"}
TiII  = {"El":2.5903,"gf":-1.410,"lam":521.1530,"dep":0.362,"name":"Ti II"}
#Blnd  = {"El":-1    ,"gf":0     ,"lam":5212.219,"dep":0.142,"name":"Blend"}
CoI   = {"El":3.5145,"gf":-0.110,"lam":521.2687,"dep":0.257,"name":"Co I"}
Myst  = {"El":    -1,"gf":     0,"lam":521.5571,"dep":0.181,"name":"Myst"}
CuI   = {"El":3.8167,"gf": 0.364,"lam":521.8198,"dep":0.529,"name":"Cu I"}
TiI   = {"El":0.0211,"gf":-2.220,"lam":521.9702,"dep":0.381,"name":"Fe I"}

wFeI  = [1184,1202]
wTiII = [1010,1036]
wBlnd = [941,964];
wCoI  = [888,921]; 
wMyst = [591,623]; 
wCuI  = [321,352]; 
wTiI  = [170,199]; 

sf_qu1 = spc.SpectraFactory("data/5215_aS1",framerows=802,framecols=1476)
sf_qu1.frame_row_cut([0]+list(range(671,679))+[801])
sf_qu1.frame_col_cut([0,1474,1475])
sf_qu1.contrast_cut(70)
sf_qu1.set_continua("segments")

qu1 = sf_qu1.make_spectra()
qu1m = qu1[:,:].mean(axis=0) 
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1] # Define continua for qu1 series

qu1FeI  = spc.splineline(wFeI,  FeI,  qu1.meta);qu1FeI.recenter(qu1m) # Define lines for qu1 series...
qu1TiII = spc.splineline(wTiII, TiII, qu1.meta);qu1TiII.recenter(qu1m)
#qu1Blnd  = spc.splineline(wBlnd,Blnd, qu1.meta);qu1Blnd.recenter(qu1m)
qu1CoI  = spc.splineline(wCoI,  CoI,  qu1.meta);qu1CoI.recenter(qu1m)
qu1Myst = spc.splineline(wMyst, Myst, qu1.meta);qu1Myst.recenter(qu1m)
qu1CuI  = spc.splineline(wCuI,  CuI,  qu1.meta);qu1CuI.recenter(qu1m)
qu1TiI  = spc.splineline(wTiI,  TiI,  qu1.meta);qu1TiI.recenter(qu1m)
#qu1lines = [qu1FeI,qu1TiII,qu1CoI,qu1Blnd,qu1Myst,qu1CuI,qu1TiI]
qu1lines = [qu1FeI,qu1TiII,qu1CoI,qu1Myst,qu1CuI,qu1TiI]


sf_qu2 = spc.SpectraFactory("data/5215_bS1",framerows=802,framecols=1476)
sf_qu2.frame_row_cut([0]+list(range(669,674))+[801])
sf_qu2.frame_col_cut([0,1475])
sf_qu2.contrast_cut(50)
sf_qu2.set_continua("segments")

qu2 = sf_qu2.make_spectra()
qu2con = qu2.meta.cont[0]*qu2.lmbd.mean() + qu2.meta.cont[1] # Define continua for qu2 series

qu2FeI  = spc.splineline(wFeI,  FeI,  qu2.meta) # Define lines for qu2 series...
qu2TiII = spc.splineline(wTiII, TiII, qu2.meta)
qu2CoI  = spc.splineline(wCoI,  CoI,  qu2.meta)
qu2Myst = spc.splineline(wMyst, Myst, qu2.meta)
qu2CuI  = spc.splineline(wCuI,  CuI,  qu2.meta)
qu2TiI  = spc.splineline(wTiI,  TiI,  qu2.meta)
qu2lines = [qu2FeI,qu2TiII,qu2CoI,qu2Myst,qu2CuI,qu2TiI]


wFeI  = [1253,1272]; cFeI  = 520.9892
wTiII = [1077,1110]; cTiII = 521.1535
wCoI  = [956,996];   cCoI  = 521.2691
wMyst = [659,691];   cMyst = 521.5571
wCuI  = [385,424];   cCuI  = 521.8209
wTiI  = [234,278];   cTiI  = 521.9706

sf_spt = spc.SpectraFactory("data/5215_aS2",framerows=774,framecols=1458)
sf_spt.frame_row_cut([0]+list(range(666,675))+[772,773])
sf_spt.frame_col_cut([0])
sf_spt.contrast_cut(50)
sf_spt.set_continua("segments")

spt = sf_spt.make_spectra()
sptm = qu1m # spt[:,:].mean(axis=0) 
sptcon = spt.meta.cont[0]*spt.lmbd.mean() + spt.meta.cont[1] # Define continua for spt series

sptFeI  = spc.splineline(wFeI,  FeI,  spt.meta);sptFeI.recenter(sptm) # Define lines for spt series... 
sptTiII = spc.splineline(wTiII, TiII, spt.meta);sptTiII.recenter(sptm)
sptCoI  = spc.splineline(wCoI,  CoI,  spt.meta);sptCoI.recenter(sptm)
sptMyst = spc.splineline(wMyst, Myst, spt.meta);sptMyst.recenter(sptm)
sptCuI  = spc.splineline(wCuI,  CuI,  spt.meta);sptCuI.recenter(sptm)
sptTiI  = spc.splineline(wTiI,  TiI,  spt.meta);sptTiI.recenter(sptm)
sptlines = [sptFeI,sptTiII,sptCoI,sptMyst,sptCuI,sptTiI]


um = 0.33427
wl = 0.55930
pn = 0.82905
xspotlims = (521.41521357925342, 521.68125806858586)
yspotlims = (0.63457553733797023, 1.1250485137224484)

qu1lims = {}
qu1lims["ewlim"]   = ( 0.5 , 1.6 )
qu1lims["vellim"]  = (-2.75, 4   )
qu1lims["rellim"]  = ( 0.32, 1.3 )
qu1lims["fw13lim"] = ( 0.05, 0.7 )
qu1lims["fwhmlim"] = ( 0.10, 0.90)
qu1lims["fw23lim"] = ( 0.10, 1.05)
qu1lims["as13lim"] = (-0.16, 0.16)
qu1lims["as12lim"] = (-0.16, 0.16)
qu1lims["as23lim"] = (-0.16, 0.16)

qu2lims = {}
qu2lims["ewlim"]   = ( 0.5 , 1.6 )
qu2lims["vellim"]  = (-2   , 3   )
qu2lims["rellim"]  = ( 0.4 , 1.01)
qu2lims["fw13lim"] = ( 0.10, 0.35)
qu2lims["fwhmlim"] = ( 0.15, 0.45)
qu2lims["fw23lim"] = ( 0.15, 0.55)
qu2lims["as13lim"] = (-0.03, 0.04)
qu2lims["as12lim"] = (-0.04, 0.04)
qu2lims["as23lim"] = (-0.04, 0.05)

sptlims = {}
sptlims["ewlim"]   = ( 0  , 3.3 )
sptlims["vellim"]  = (-6  , 9   )
sptlims["rellim"]  = ( 0.2, 1.  )
sptlims["fw13lim"] = ( 0.0, 1.0 )
sptlims["fwhmlim"] = ( 0.0, 1.0 )
sptlims["fw23lim"] = ( 0.0, 1.2 )
sptlims["as13lim"] = (-0.3, 0.2 )
sptlims["as12lim"] = (-0.3, 0.2 )
sptlims["as23lim"] = (-0.3, 0.2 )
