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
as2 = sf56_as2.make_spectra()
as2con = as2.meta.cont[0]*as2.lmbd.mean() + as2.meta.cont[1]

as2FeIp = spc.splineline(wFeIp, cFeIp, as2.meta,"FeIp ")
as2Myst = spc.splineline(wMyst, cMyst, as2.meta,"Myst ")
as2SiI  = spc.splineline(wSiI , cSiI , as2.meta,"SiI " )
as2VI   = spc.splineline(wVI  , cVI  , as2.meta,"VI "  )
as2ScII = spc.splineline(wScII, cScII, as2.meta,"ScII ")
as2FeI  = spc.splineline(wFeIp, cFeI , as2.meta,"FeI " )

as2lines = [as2FeIp,as2Myst,as2SiI,as2VI,as2ScII,as2FeI]

sf56_cs2 = spc.SpectraFactory("data/5654_cS2",framerows=756,framecols=1480)
sf56_cs2.frame_row_cut([0,755])
sf56_cs2.frame_col_cut([0])
sf56_cs2.contrast_cut(20)
sf56_cs2.set_continua("segments")
cs2 = sf56_cs2.make_spectra()
cs2con = cs2.meta.cont[0]*cs2.lmbd.mean() + cs2.meta.cont[1]

cs2FeIp = spc.splineline(wFeIp, cFeIp, cs2.meta,"FeIp ")
cs2Myst = spc.splineline(wMyst, cMyst, cs2.meta,"Myst ")
cs2SiI  = spc.splineline(wSiI , cSiI , cs2.meta,"SiI " )
cs2VI   = spc.splineline(wVI  , cVI  , cs2.meta,"VI "  )
cs2ScII = spc.splineline(wScII, cScII, cs2.meta,"ScII ")
cs2FeI  = spc.splineline(wFeIp, cFeI , cs2.meta,"FeI " )

cs2lines = [cs2FeIp,cs2Myst,cs2SiI,cs2VI,cs2ScII,cs2FeI]

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
bs2 = sf56_bs2.make_spectra()
bs2con = bs2.meta.cont[0]*bs2.lmbd.mean() + bs2.meta.cont[1]

bs2FeIp = spc.splineline(wFeIp, cFeIp, bs2.meta,"FeIp ")
bs2Myst = spc.splineline(wMyst, cMyst, bs2.meta,"Myst ")
bs2SiI  = spc.splineline(wSiI , cSiI , bs2.meta,"SiI " )
bs2VI   = spc.splineline(wVI  , cVI  , bs2.meta,"VI "  )
bs2ScII = spc.splineline(wScII, cScII, bs2.meta,"ScII ")
bs2FeI  = spc.splineline(wFeIp, cFeI , bs2.meta,"FeI " )

bs2lines = [bs2FeIp,bs2Myst,bs2SiI,bs2VI,bs2ScII,bs2FeI]

umbra    = 0.28
wall     = (0.28,0.60)
penumbra = (0.60,0.89)
quiet    = (0.89)
as2um    = s2.make_spectra_subset(as2,rowsubset=(as2con < umbra),desc="Umbra subset")
as2wl    = s2.make_spectra_subset(as2,rowsubset=((as2con >= wall[0]) & (as2con <= wall[1])),desc="Umbra/penumbra wall")
as2pn    = s2.make_spectra_subset(as2,rowsubset=((as2con > penumbra[0]) & (as2con < penumbra[1])),desc="Penumbra")
as2qu    = s2.make_spectra_subset(as2,rowsubset=(as2con >= quiet),desc="Quiet sun")

