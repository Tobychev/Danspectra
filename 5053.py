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

as1 = sf53_as1.make_spectra()

sf53_bs1 = spc.SpectraFactory("data/5053_bS1",framerows=792,framecols=1466)
sf53_bs1.frame_row_cut([0,791])
sf53_bs1.frame_col_cut([0])
sf53_bs1.contrast_cut(80)
sf53_bs1.set_continua("segments")

bs1 = sf53_bs1.make_spectra()

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

as2 = sf53_as2.make_spectra()
