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

sf52_bs1 = spc.SpectraFactory("data/5215_bS1",framerows=802,framecols=1476)
sf52_bs1.frame_row_cut([0,801])
sf52_bs1.frame_col_cut([0])
sf52_bs1.contrast_cut(50)
sf52_bs1.set_continua("segments")

bs1 = sf52_bs1.make_spectra()

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
