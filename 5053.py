import matplotlib.pyplot as pl
import numpy as np
import spectra as spc
import visualize as vis

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

sf53_as2 = spc.SpectraFactory("data/5053_aS2",framerows=780,framecols=1486)
sf53_as2.frame_row_cut([0,791])
sf53_as2.frame_col_cut([0])
sf53_as2.contrast_cut(80)
sf53_as2.set_continua("segments")

as2 = sf53_as2.make_spectra()
