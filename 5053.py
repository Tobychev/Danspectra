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

lims = {}
lims["ewlim"]   = ( 0.5 , 1.5  )
lims["vellim"]  = (-6 , 7  )
lims["rellim"]  = ( 0.2 , 1.3  )
lims["fw13lim"] = ( -0.1 , 1.1  )
lims["fwhmlim"] = ( -0.1 , 1.1  )
lims["fw23lim"] = (-0.2 , 1.5  )
lims["as13lim"] = (-0.8 , 0.8  )
lims["as12lim"] = (-0.7 , 0.8  )
lims["as23lim"] = (-0.9 , 0.6  )
as1lims = lims

sf53_as1 = spc.SpectraFactory("data/5053_aS1",framerows=792,framecols=1466)
sf53_as1.frame_row_cut([0,791])
sf53_as1.frame_col_cut([0])
sf53_as1.contrast_cut(80)
sf53_as1.set_continua("segments")

as1 = sf53_as1.make_spectra()
as1con = as1.meta.cont[0]*as1.lmbd.mean() + as1.meta.cont[1] # Define continua for as1 series

as1VI = spc.splineline(wVI, cVI, as1.meta,"VI ") # Define lines for as1 series...
as1NiI = spc.splineline(wNiI, cNiI, as1.meta,"NiI ")
as1C2  = spc.splineline(wC2 , cC2 , as1.meta,"C2 " )
as1CI   = spc.splineline(wCI  , cCI  , as1.meta,"CI "  )
as1Myst = spc.splineline(wMyst, cMyst, as1.meta,"Myst ")
as1FeMg  = spc.splineline(wFeMg, cFeMg , as1.meta,"FeMg " )
as1FeI  = spc.splineline(wFeI, cFeI , as1.meta,"FeI " )

as1lines = [as1VI,as1NiI,as1C2,as1CI,as1Myst,as1FeMg,as1FeI]

sf53_bs1 = spc.SpectraFactory("data/5053_bS1",framerows=792,framecols=1466)
sf53_bs1.frame_row_cut([0,791])
sf53_bs1.frame_col_cut([0])
sf53_bs1.contrast_cut(80)
sf53_bs1.set_continua("segments")

bs1 = sf53_bs1.make_spectra()
bs1con = bs1.meta.cont[0]*bs1.lmbd.mean() + bs1.meta.cont[1] # Define continua for bs1 series

bs1VI = spc.splineline(wVI, cVI, bs1.meta,"VI ") # Define lines for bs1 series...
bs1NiI = spc.splineline(wNiI, cNiI, bs1.meta,"NiI ")
bs1C2  = spc.splineline(wC2 , cC2 , bs1.meta,"C2 " )
bs1CI   = spc.splineline(wCI  , cCI  , bs1.meta,"CI "  )
bs1Myst = spc.splineline(wMyst, cMyst, bs1.meta,"Myst ")
bs1FeMg  = spc.splineline(wFeMg, cFeMg , bs1.meta,"FeMg " )
bs1FeI  = spc.splineline(wFeI, cFeI , bs1.meta,"FeI " )

bs1lines = [bs1VI,bs1NiI,bs1C2,bs1CI,bs1Myst,bs1FeMg,bs1FeI]

wVI   = [1356,1413]; cVI   = 504.7302
wNiI  = [1229,1269]; cNiI  = 504.8853
wC2   = [1183,1200]; cC2   = 504.9425
wCI   = [898,928];   cCI   = 505.2151
wMyst = [754,787];   cMyst = 505.3577
wFeMg = [421,451];   cFeMg = 505.6846
wFeI  = [258,283];    cFeI  = 505.8495

lims["ewlim"]   = ( 0.5 , 1.5  )
lims["vellim"]  = (-6 , 7  )
lims["rellim"]  = ( 0.2 , 1.3  )
lims["fw13lim"] = ( -0.1 , 0.9  )
lims["fwhmlim"] = ( -0.1 , 1.3  )
lims["fw23lim"] = (-0.2 , 1.6  )
lims["as13lim"] = (-0.6 , 0.4  )
lims["as12lim"] = (-0.7 , 0.6  )
lims["as23lim"] = (-0.9 , 0.6  )
as2lims = lims

sf53_as2 = spc.SpectraFactory("data/5053_aS2",framerows=780,framecols=1486)
sf53_as2.frame_row_cut([0,791])
sf53_as2.frame_col_cut([0])
sf53_as2.contrast_cut(80)
sf53_as2.set_continua("segments")

as2 = sf53_as2.make_spectra()
as2con = as2.meta.cont[0]*as2.lmbd.mean() + as2.meta.cont[1] # Define continua for as2 series

as2VI = spc.splineline(wVI, cVI, as2.meta,"VI ") # Define lines for as2 series...
as2NiI = spc.splineline(wNiI, cNiI, as2.meta,"NiI ")
as2C2  = spc.splineline(wC2 , cC2 , as2.meta,"C2 " )
as2CI   = spc.splineline(wCI  , cCI  , as2.meta,"CI "  )
as2Myst = spc.splineline(wMyst, cMyst, as2.meta,"Myst ")
as2FeMg  = spc.splineline(wFeMg, cFeMg , as2.meta,"FeMg " )
as2FeI  = spc.splineline(wFeI, cFeI , as2.meta,"FeI " )

as2lines = [as2VI,as2NiI,as2CI,as2Myst,as2FeI] #C2 and FeMg cause a crash
#as2lines = [as2VI,as2NiI,as2C2,as2CI,as2Myst,as2FeMg,as2FeI]
