import matplotlib.pyplot as pl
import spectra as spc

sf64_as1 = spc.SpectraFactory("data/6405_aS1")
sf64_as1.frame_row_cut([0,799])
sf64_as1.contrast_cut(50)
wFlat = [1179,1239]
as1  = sf.make_spectra()
cont_as1 = as1[:,range(wFlat[0],wFlat[1])].mean(axis=1)

sf64_as2 = spc.SpectraFactory("data/6405_aS2",framerows=774,framecols=1446)
sf64_as2.frame_row_cut([0,743])
sf64_as2.frame_col_cut([0])
sf64_as2.contrast_cut(85)
as2 = sf64_as2.make_spectra()
cont_as1 = as1[:,range(wFlat[0],wFlat[1])].mean(axis=1)

sf64_as2 = spc.SpectraFactory("data/6405_bS1")
sf64_as2.frame_row_cut([0,799])
sf64_as2.contrast_cut(50)
bS1 = sf64_as2.make_spectra()
cont_as1 = as1[:,range(wFlat[0],wFlat[1])].mean(axis=1)


