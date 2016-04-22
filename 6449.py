import matplotlib.pyplot as pl
import numpy as np
import spectra as spc
import visualize as vis

sf64_as1 = spc.SpectraFactory("data/6449_aS1",framerows=802,framecols=1514)
sf64_as1.frame_row_cut([0,801])
sf64_as1.frame_col_cut([0])
sf64_as1.contrast_cut(20)
sf64_as1.set_continua("segments")

as1 = sf64_as1.make_spectra()

sf64_as2 = spc.SpectraFactory("data/6449_aS2",framerows=780,framecols=1438)
sf64_as2.frame_row_cut([0,779])
sf64_as2.frame_col_cut([0])
sf64_as2.contrast_cut(50)
sf64_as2.set_continua("segments")

as2 = sf64_as2.make_spectra()

sf64_bs1 = spc.SpectraFactory("data/6449_bS1",framerows=802,framecols=1514)
sf64_bs1.frame_row_cut([0,779])
sf64_bs1.frame_col_cut([0])
sf64_bs1.contrast_cut(20)
sf64_bs1.set_continua("segments")

bs1 = sf64_bs1.make_spectra()
