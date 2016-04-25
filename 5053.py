import matplotlib.pyplot as pl
import numpy as np
import spectra as spc
import visualize as vis
wC22  = [876,897];   cC22  = 505.0737

wVI   = [1216,1244]; cVI   = 504.7302
wNiI  = [1057,1100]; cNiI  = 504.8853
wC2   = [1015,1032]; cC2   = 504.9425
wCI   = [728,761];   cCI   = 505.2151
wMyst = [580,620];   cMyst = 505.3577
wFeMg = [252,280];   cFeMg = 505.6846
wFeI  = [87,117];    cFeI  = 505.8495

sf53_as1 = spc.SpectraFactory("data/5053_aS1",framerows=792,framecols=1466)
sf53_as1.frame_row_cut([0,791])
sf53_as1.frame_col_cut([0])
sf53_as1.contrast_cut(80)
sf53_as1.set_continua("segments")

qu1 = sf53_as1.make_spectra()
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1] # Define continua for qu1 series

qu1VI   = spc.splineline(wVI,   cVI,   qu1.meta,"VI ") # Define lines for qu1 series...
qu1NiI  = spc.splineline(wNiI,  cNiI,  qu1.meta,"NiI ")
qu1C2   = spc.splineline(wC2,   cC2,   qu1.meta,"C2 ")
qu1CI   = spc.splineline(wCI,   cCI,   qu1.meta,"CI ")
qu1Myst = spc.splineline(wMyst, cMyst, qu1.meta,"Myst ")
qu1FeMg = spc.splineline(wFeMg, cFeMg, qu1.meta,"FeMg ")
qu1FeI  = spc.splineline(wFeI,  cFeI,  qu1.meta,"FeI ")

qu1lines = [qu1VI,qu1NiI,qu1C2,qu1CI,qu1Myst,qu1FeMg,qu1FeI]

lims = {}
lims["ewlim"]   = ( 0.5 , 1.5)
lims["vellim"]  = (-6   , 7  )
lims["rellim"]  = ( 0.2 , 1.3)
lims["fw13lim"] = (-0.1 , 1.1)
lims["fwhmlim"] = (-0.1 , 1.1)
lims["fw23lim"] = (-0.2 , 1.5)
lims["qu13lim"] = (-0.8 , 0.8)
lims["qu12lim"] = (-0.7 , 0.8)
lims["as23lim"] = (-0.9 , 0.6)
qu1lims = lims

sf53_bs1 = spc.SpectraFactory("data/5053_bS1",framerows=792,framecols=1466)
sf53_bs1.frame_row_cut([0,791])
sf53_bs1.frame_col_cut([0])
sf53_bs1.contrast_cut(80)
sf53_bs1.set_continua("segments")

qu2 = sf53_bs1.make_spectra()
qu2con = qu2.meta.cont[0]*qu2.lmbd.mean() + qu2.meta.cont[1] # Define continua for qu2 series

qu2VI   = spc.splineline(wVI,   cVI,   qu2.meta,"VI ") # Define lines for qu2 series...
qu2NiI  = spc.splineline(wNiI,  cNiI,  qu2.meta,"NiI ")
qu2C2   = spc.splineline(wC2,   cC2,   qu2.meta,"C2 ")
qu2CI   = spc.splineline(wCI,   cCI,   qu2.meta,"CI ")
qu2Myst = spc.splineline(wMyst, cMyst, qu2.meta,"Myst ")
qu2FeMg = spc.splineline(wFeMg, cFeMg, qu2.meta,"FeMg ")
qu2FeI  = spc.splineline(wFeI,  cFeI,  qu2.meta,"FeI ")

qu2lines = [qu2VI,qu2NiI,qu2C2,qu2CI,qu2Myst,qu2FeMg,qu2FeI]
qu2lims = lims

wVI   = [1356,1413]; cVI   = 504.7302
wNiI  = [1229,1269]; cNiI  = 504.8853
wC2   = [1183,1200]; cC2   = 504.9425
wCI   = [898,928];   cCI   = 505.2151
wMyst = [754,787];   cMyst = 505.3577
wFeMg = [421,451];   cFeMg = 505.6846
wFeI  = [258,283];    cFeI  = 505.8495

sf53_as2 = spc.SpectraFactory("data/5053_aS2",framerows=780,framecols=1486)
sf53_as2.frame_row_cut([0,791])
sf53_as2.frame_col_cut([0])
sf53_as2.contrast_cut(80)
sf53_as2.set_continua("segments")

spt = sf53_as2.make_spectra()
sptcon = spt.meta.cont[0]*spt.lmbd.mean() + spt.meta.cont[1] # Define continua for spt series

sptVI   = spc.splineline(wVI,   cVI,   spt.meta,"VI ") # Define lines for spt series...
sptNiI  = spc.splineline(wNiI,  cNiI,  spt.meta,"NiI ")
sptC2   = spc.splineline(wC2,   cC2,   spt.meta,"C2 ")
sptCI   = spc.splineline(wCI,   cCI,   spt.meta,"CI ")
sptMyst = spc.splineline(wMyst, cMyst, spt.meta,"Myst ")
sptFeMg = spc.splineline(wFeMg, cFeMg, spt.meta,"FeMg ")
sptFeI  = spc.splineline(wFeI,  cFeI,  spt.meta,"FeI ")

sptlines = [sptVI,sptNiI,sptCI,sptMyst,sptFeI] #C2 and FeMg cause a crash
#sptlines = [sptVI,sptNiI,sptC2,sptCI,sptMyst,sptFeMg,sptFeI]

lims["ewlim"]   = ( 0.5 , 1.5)
lims["vellim"]  = (-6   , 7  )
lims["rellim"]  = ( 0.2 , 1.3)
lims["fw13lim"] = (-0.1 , 0.9)
lims["fwhmlim"] = (-0.1 , 1.3)
lims["fw23lim"] = (-0.2 , 1.6)
lims["as13lim"] = (-0.6 , 0.4)
lims["as12lim"] = (-0.7 , 0.6)
lims["spt3lim"] = (-0.9 , 0.6)
sptlims = lims
