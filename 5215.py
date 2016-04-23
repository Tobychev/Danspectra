import matplotlib.pyplot as pl
import numpy as np
import spectra as spc
import visualize as vis

wFeI  = [1183,1202]; cFeI  = 520.9892
wFeTi = [1010,1036]; cFeTi = 521.1535
wCoI  = [888,924];   cCoI  = 521.2691
wMyst = [591,623];   cMyst = 521.5571
wCuI  = [321,352];   cCuI  = 521.8209
wTiI  = [170,199];   cTiI  = 521.9706

sf52_as1 = spc.SpectraFactory("data/5215_aS1",framerows=802,framecols=1476)
sf52_as1.frame_row_cut([0,801])
sf52_as1.frame_col_cut([0])
sf52_as1.contrast_cut(70)
sf52_as1.set_continua("segments")

as1 = sf52_as1.make_spectra()
as1con = as1.meta.cont[0]*as1.lmbd.mean() + as1.meta.cont[1] # Define continua for as1 series

as1FeI = spc.splineline(wFeI, cFeI, as1.meta,"FeI ") # Define lines for as1 series...
as1FeTi = spc.splineline(wFeTi, cFeTi, as1.meta,"FeTi ")
as1CoI  = spc.splineline(wCoI , cCoI , as1.meta,"CoI " )
as1Myst   = spc.splineline(wMyst  , cMyst  , as1.meta,"Myst "  )
as1CuI = spc.splineline(wCuI, cCuI, as1.meta,"CuI ")
as1TiI  = spc.splineline(wTiI, cTiI , as1.meta,"TiI " )

as1lines = [as1FeI,as1FeTi,as1CoI,as1Myst,as1CuI,as1TiI]

sf52_bs1 = spc.SpectraFactory("data/5215_bS1",framerows=802,framecols=1476)
sf52_bs1.frame_row_cut([0,801])
sf52_bs1.frame_col_cut([0])
sf52_bs1.contrast_cut(50)
sf52_bs1.set_continua("segments")

bs1 = sf52_bs1.make_spectra()
bs1con = bs1.meta.cont[0]*bs1.lmbd.mean() + bs1.meta.cont[1] # Define continua for bs1 series

bs1FeI = spc.splineline(wFeI, cFeI, bs1.meta,"FeI ") # Define lines for bs1 series...
bs1FeTi = spc.splineline(wFeTi, cFeTi, bs1.meta,"FeTi ")
bs1CoI  = spc.splineline(wCoI , cCoI , bs1.meta,"CoI " )
bs1Myst   = spc.splineline(wMyst  , cMyst  , bs1.meta,"Myst "  )
bs1CuI = spc.splineline(wCuI, cCuI, bs1.meta,"CuI ")
bs1TiI  = spc.splineline(wTiI, cTiI , bs1.meta,"TiI " )

bs1lines = [bs1FeI,bs1FeTi,bs1CoI,bs1Myst,bs1CuI,bs1TiI]

wFeI  = [1253,1272]; cFeI  = 520.9892
wFeTi = [1077,1110]; cFeTi = 521.1535
wCoI  = [956,996];   cCoI  = 521.2691
wMyst = [659,691];   cMyst = 521.5571
wCuI  = [385,424];   cCuI  = 521.8209
wTiI  = [234,278];   cTiI  = 521.9706

sf52_as2 = spc.SpectraFactory("data/5215_aS2",framerows=774,framecols=1458)
sf52_as2.frame_row_cut([0,773])
sf52_as2.frame_col_cut([0])
sf52_as2.contrast_cut(70)
sf52_as2.set_continua("segments")

as2 = sf52_as2.make_spectra()
as2con = as2.meta.cont[0]*as2.lmbd.mean() + as2.meta.cont[1] # Define continua for as2 series

as2FeI = spc.splineline(wFeI, cFeI, as2.meta,"FeI ") # Define lines for as2 series...
as2FeTi = spc.splineline(wFeTi, cFeTi, as2.meta,"FeTi ")
as2CoI  = spc.splineline(wCoI , cCoI , as2.meta,"CoI " )
as2Myst   = spc.splineline(wMyst  , cMyst  , as2.meta,"Myst "  )
as2CuI = spc.splineline(wCuI, cCuI, as2.meta,"CuI ")
as2TiI  = spc.splineline(wTiI, cTiI , as2.meta,"TiI " )

as2lines = [as2FeI,as2FeTi,as2CoI,as2Myst,as2CuI,as2TiI]
