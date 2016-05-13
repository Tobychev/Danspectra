import spectra as spc
import numpy as np

NiI  = {"El":3.8474,"gf":-0.380,"lam":504.8847,"dep":0.580,"name":"Ni I"}
CI   = {"El":7.6848,"gf":-1.303,"lam":505.2144,"dep":0.141,"name":"C I"}
C2   = {"El":0.4891,"gf": 0.158,"lam":505.2620,"dep":0.020,"name":"C_2 "}
TiI  = {"El":2.1747,"gf":-0.270,"lam":505.2868,"dep":0.277,"name":"Ti I"}
Myst = {"El":    -1,"gf":     0,"lam":505.3577,"dep":0.111,"name":"Myst"}
FeI  = {"El":3.6398,"gf":-1.921,"lam":505.4642,"dep":0.476,"name":"Fe I"}
FeI2 = {"El":3.6417,"gf":-2.830,"lam":505.8496,"dep":0.130,"name":"Fe I"}

wNiI  = [1064,1094]
wCI   = [731,761]
wC2   = [690,701]
wTiI  = [666,676]
wMyst = [580,620]
wFeI  = [480,504]
wFeI2 = [87,114]

sf_qu1 = spc.SpectraFactory("data/5053_aS1",framerows=792,framecols=1466)
sf_qu1.frame_row_cut([0]+list(range(666,679))+[791])
sf_qu1.frame_col_cut([0,1,1465])
sf_qu1.contrast_cut(80)
sf_qu1.set_continua("segments")

## Select subset
lowNiIcut = -0.511135 # Read off from Ni I vel
qu1sel = np.arange(0,3107) # Givin pretty pure cut

qu1 = sf_qu1.make_spectra()
qu1.modify(lambda x : x[qu1sel,:]) # Cutting
qu1.meta.cont = (qu1.meta.cont[0][qu1sel], qu1.meta.cont[1][qu1sel]) # Updating continua

qu1m = qu1[:,:].mean(axis=0) 
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1] # Define continua for qu1 series

qu1NiI  = spc.splineline(wNiI,  NiI,  qu1.meta);qu1NiI.recenter(qu1m)
qu1CI   = spc.splineline(wCI,   CI,   qu1.meta);qu1CI.recenter(qu1m)
qu1C2   = spc.splineline(wC2,   C2,   qu1.meta);qu1C2.recenter(qu1m)
qu1TiI  = spc.splineline(wTiI,  TiI,  qu1.meta);qu1TiI.recenter(qu1m)
qu1Myst = spc.splineline(wMyst, Myst, qu1.meta);qu1Myst.recenter(qu1m)
qu1FeI  = spc.splineline(wFeI,  FeI,  qu1.meta);qu1FeI.recenter(qu1m)
qu1FeI2 = spc.splineline(wFeI2, FeI2, qu1.meta);qu1FeI2.recenter(qu1m)
qu1lines = [qu1NiI,qu1CI,qu1C2,qu1TiI,qu1Myst,qu1FeI,qu1FeI2]


wNiI  = [1065,1099]
wCI   = [732,764]; 
wC2   = [693,705]; 
wTiI  = [669,679]; 
wMyst = [588,622]; 
wFeI  = [485,505]; 
wFeI2 = [92,119];  

sf_qu2 = spc.SpectraFactory("data/5053_bS1",framerows=792,framecols=1466)
sf_qu2.frame_row_cut([0]+list(range(666,671))+[791])
sf_qu2.frame_col_cut([0,1465])
sf_qu2.contrast_cut(80)
sf_qu2.set_continua("segments")

qu2 = sf_qu2.make_spectra()
qu2con = qu2.meta.cont[0]*qu2.lmbd.mean() + qu2.meta.cont[1] # Define continua for qu2 series

qu2NiI  = spc.splineline(wNiI,  NiI,  qu1.meta)
qu2CI   = spc.splineline(wCI,   CI,   qu1.meta)
qu2C2   = spc.splineline(wC2,   C2,   qu1.meta)
qu2TiI  = spc.splineline(wTiI,  TiI,  qu1.meta)
qu2Myst = spc.splineline(wMyst, Myst, qu1.meta)
qu2FeI  = spc.splineline(wFeI,  FeI,  qu1.meta)
qu2FeI2 = spc.splineline(wFeI2, FeI2, qu1.meta)
qu2lines = [qu2NiI,qu2CI,qu2C2,qu2TiI,qu2Myst,qu2FeI,qu2FeI2]


wNiI  = [1228,1269]; cNiI  = 504.8853
wCI   = [901,927];   cCI   = 505.2151
wC2   = [859,870];   cC2   = 505.261
wTiI  = [834,844];   cTiI  = 505.2869
wMyst = [757,786];   cMyst = 505.3577
wFeI  = [649,674];   cFeI  = 505.4647
wFeI2 = [258,285];   cFeI2 = 505.8495

sf_spt = spc.SpectraFactory("data/5053_aS2",framerows=780,framecols=1486)
sf_spt.frame_row_cut([0,1,779])
sf_spt.frame_col_cut([0]+list(range(658,672))+[1485])
sf_spt.contrast_cut(80)
sf_spt.set_continua("segments")

spt = sf_spt.make_spectra()
sptcon = spt.meta.cont[0]*spt.lmbd.mean() + spt.meta.cont[1] # Define continua for spt series

sptNiI  = spc.splineline(wNiI,  NiI,  spt.meta)
sptCI   = spc.splineline(wCI,   CI,   spt.meta)
sptC2   = spc.splineline(wC2,   C2,   spt.meta)
sptTiI  = spc.splineline(wTiI,  TiI,  spt.meta)
sptMyst = spc.splineline(wMyst, Myst, spt.meta)
sptFeI  = spc.splineline(wFeI,  FeI,  spt.meta)
sptFeI2 = spc.splineline(wFeI2, FeI2, spt.meta)
sptlines = [sptNiI,sptCI,sptC2,sptTiI,sptMyst,sptFeI,sptFeI2]

um = 0.32048
wl = 0.58196
pn = 0.82306
xspotlims = (505.26553248218875, 505.44646047080431)
yspotlims = (0.63975783972631528, 1.0820081717200596)

qu1lims = {}
qu1lims["ewlim"]   = ( 0.5 ,1.5 )
qu1lims["vellim"]  = (-6   ,4   )
qu1lims["rellim"]  = ( 0.25,0.99)
qu1lims["fw13lim"] = (-0.05,1.06)
qu1lims["fwhmlim"] = (-0.05,1.06)
qu1lims["fw23lim"] = (-0.05,1.06)
qu1lims["as13lim"] = (-0.62,0.45)
qu1lims["as12lim"] = (-0.62,0.62)
qu1lims["as23lim"] = (-0.68,0.50)

qu2lims = {}
qu2lims["ewlim"]   = ( 0.5 ,1.5 )
qu2lims["vellim"]  = (-4   ,3   )
qu2lims["rellim"]  = ( 0.75,0.95)
qu2lims["fw13lim"] = ( 0.05,0.55)
qu2lims["fwhmlim"] = ( 0.05,0.7 )
qu2lims["fw23lim"] = ( 0.05,0.9 )
qu2lims["as13lim"] = (-0.15,0.1 )
qu2lims["as12lim"] = (-0.15,0.1 )
qu2lims["as23lim"] = (-0.15,0.15)

sptlims = {}
sptlims["ewlim"]   = ( 0.0 , 2.5 )
sptlims["vellim"]  = (-6   , 7   )
sptlims["rellim"]  = ( 0.2 , 1.3 )
sptlims["fw13lim"] = ( 0.05, 0.45)
sptlims["fwhmlim"] = ( 0.1 , 1.0 )
sptlims["fw23lim"] = ( 0.1 , 1.1 )
sptlims["as13lim"] = (-0.2 , 0.4 )
sptlims["as12lim"] = (-0.2 , 0.6 )
sptlims["as23lim"] = (-0.17, 0.2 )
