import matplotlib.pyplot as pl
import numpy as np
import spectra as spc
import visualize as vis

#Hand selected limits
wFeIp = [1216,1248]; cFeIp = 565.1477
wMyst = [873,937];   cMyst = 565.4501
wSiI  = [842,867];   cSiI  = 565.4937
wVI   = [567,595];   cVI   = 565.7456
wScII = [512,548];   cScII = 565.7880
wFeI  = [131,168];   cFeI  = 566.1354

sf56_as2 = spc.SpectraFactory("data/5654_aS2",framerows=756,framecols=1480)
sf56_as2.frame_row_cut([0,755])
sf56_as2.frame_col_cut([0])
sf56_as2.contrast_cut(20)
sf56_as2.set_continua("segments")

qu1 = sf56_as2.make_spectra() # Do this before defining lines
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1] # Define continua for qu1 series

qu1FeIp = spc.splineline(wFeIp, cFeIp, qu1.meta,"FeIp ") # Define lines for qu1 series...
qu1Myst = spc.splineline(wMyst, cMyst, qu1.meta,"Myst ")
qu1SiI  = spc.splineline(wSiI , cSiI , qu1.meta,"SiI " )
qu1VI   = spc.splineline(wVI  , cVI  , qu1.meta,"VI "  )
qu1ScII = spc.splineline(wScII, cScII, qu1.meta,"ScII ")
qu1FeI  = spc.splineline(wFeI, cFeI , qu1.meta,"FeI "  )
qu1lines = [qu1FeIp,qu1Myst,qu1SiI,qu1VI,qu1ScII,qu1FeI] # ... and save in this list for qu1

lims = {}
lims["ewlim"]   = ( 0.3 , 1.8 )
lims["vellim"]  = (-5.8 , 6.1 )
lims["rellim"]  = ( 0.2 , 1.5 )
lims["fw13lim"] = (-0.1 , 1.0 )
lims["fwhmlim"] = (-0.1 , 1.0 )
lims["fw23lim"] = (-0.2 , 1.5 )
lims["as13lim"] = (-0.2 , 0.2 )
lims["as12lim"] = (-0.3 , 0.2 )
lims["as23lim"] = (-0.7 , 0.5 )
qu1lims = lims

sf56_cs2 = spc.SpectraFactory("data/5654_cS2",framerows=756,framecols=1480)
sf56_cs2.frame_row_cut([0,755])
sf56_cs2.frame_col_cut([0])
sf56_cs2.contrast_cut(20)
sf56_cs2.set_continua("segments")

qu2 = sf56_cs2.make_spectra()
qu2con = qu2.meta.cont[0]*qu2.lmbd.mean() + qu2.meta.cont[1] # Define continua for qu2 series

qu2FeIp = spc.splineline(wFeIp, cFeIp, qu2.meta,"FeIp ") # Define lines for qu2 series...
qu2Myst = spc.splineline(wMyst, cMyst, qu2.meta,"Myst ")
qu2SiI  = spc.splineline(wSiI , cSiI , qu2.meta,"SiI " )
qu2VI   = spc.splineline(wVI  , cVI  , qu2.meta,"VI "  )
qu2ScII = spc.splineline(wScII, cScII, qu2.meta,"ScII ")
qu2FeI  = spc.splineline(wFeI, cFeI , qu2.meta,"FeI "  )
qu2lines = [qu2FeIp,qu2Myst,qu2SiI,qu2VI,qu2ScII,qu2FeI] # ... and save in this list for qu2

qu2lims=lims

# Updating the limits to better match the meanspectra of this series
wFeIp = [1217,1253]; cFeIp = 565.1477
wMyst = [873,937];   cMyst = 565.4501
wSiI  = [844,871];   cSiI  = 565.4937
wVI   = [564,596];   cVI   = 565.7456
wScII = [512,555];   cScII = 565.7880
wFeI  = [139,167];   cFeI  = 566.1354


sf56_bs2 = spc.SpectraFactory("data/5654_bS2",framerows=756,framecols=1480)
sf56_bs2.frame_row_cut([0,755])
sf56_bs2.frame_col_cut([0])
sf56_bs2.contrast_cut(50)
sf56_bs2.set_continua("segments")
spt = sf56_bs2.make_spectra()
sptcon = spt.meta.cont[0]*spt.lmbd.mean() + spt.meta.cont[1] # Define continua for spt series

sptFeIp = spc.splineline(wFeIp, cFeIp, spt.meta,"FeIp ") # Define lines for spt series...
sptMyst = spc.splineline(wMyst, cMyst, spt.meta,"Myst ")
sptSiI  = spc.splineline(wSiI , cSiI , spt.meta,"SiI " )
sptVI   = spc.splineline(wVI  , cVI  , spt.meta,"VI "  )
sptScII = spc.splineline(wScII, cScII, spt.meta,"ScII ")
sptFeI  = spc.splineline(wFeI, cFeI , spt.meta,"FeI "  )
sptlines = [sptFeIp,sptMyst,sptSiI,sptVI,sptScII,sptFeI] # ... and save in this list for spt

lims["ewlim"]   = ( 0   , 5   )
lims["vellim"]  = (-5.8 , 7   )
lims["rellim"]  = ( 0   , 6   )
lims["fw13lim"] = (-0.1 , 1.0 )
lims["fwhmlim"] = (-0.1 , 1.1 )
lims["fw23lim"] = (-0.2 , 1.5 )
lims["as13lim"] = (-0.4 , 0.3 )
lims["as12lim"] = (-0.5 , 0.4 )
lims["as23lim"] = (-0.9 , 0.5 )
sptlims = lims

umbra    = 0.28
wall     = (0.28,0.60)
penumbra = (0.60,0.89)
quiet    = (0.89)


