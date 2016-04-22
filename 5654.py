import matplotlib.pyplot as pl
import numpy as np
import spectra as spc
import visualize as vis

sf56_as2 = spc.SpectraFactory("data/5654_aS2",framerows=756,framecols=1480)
sf56_as2.frame_row_cut([0,755])
sf56_as2.frame_col_cut([0])
sf56_as2.contrast_cut(20)
sf56_as2.set_continua("segments")

as2 = sf56_as2.make_spectra()

sf56_bs2 = spc.SpectraFactory("data/5654_bS2",framerows=756,framecols=1480)
sf56_bs2.frame_row_cut([0,755])
sf56_bs2.frame_col_cut([0])
sf56_bs2.contrast_cut(50)
sf56_bs2.set_continua("segments")

bs2 = sf56_bs2.make_spectra()

sf56_cs2 = spc.SpectraFactory("data/5654_cS2",framerows=756,framecols=1480)
sf56_cs2.frame_row_cut([0,755])
sf56_cs2.frame_col_cut([0])
sf56_cs2.contrast_cut(20)
sf56_cs2.set_continua("segments")

cs2 = sf56_cs2.make_spectra()
