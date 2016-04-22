import matplotlib.pyplot as pl
import numpy as np
import spectra as spc
import visualize as vis

wMyst = [1211,1327]; cMyst = 644.9127
wCa   = [1150,1198]; cCa   = 644.9820
wCoBl = [1088,1151]; cCoBl = 645.0179
wUnk2 = [834,889];   cUnk2 = 645.2315
wH2O  = [613,637];   cH2O  = 645.4139
wCo   = [449,534];   cCo   = 645.5001
wCa2  = [416,446];   cCa2  = 645.5602
wH2O2 = [8,34];      cH2O2 = 645.8892

sf64_as1 = spc.SpectraFactory("data/6449_aS1",framerows=802,framecols=1514)
sf64_as1.frame_row_cut([0,801])
sf64_as1.frame_col_cut([0])
sf64_as1.contrast_cut(50)
sf64_as1.set_continua("segments")

as1 = sf64_as1.make_spectra()

sf64_bs1 = spc.SpectraFactory("data/6449_bS1",framerows=802,framecols=1514)
sf64_bs1.frame_row_cut([0,779])
sf64_bs1.frame_col_cut([0])
sf64_bs1.contrast_cut(50)
sf64_bs1.set_continua("segments")

bs1 = sf64_bs1.make_spectra()

wMyst = [831,938]; cMyst = 644.9127
wCa   = [764,830]; cCa   = 644.9820
wCoBl = [695,763]; cCoBl = 645.0179
wUnk2 = [444,492]; cUnk2 = 645.2315
wH2O  = [227,254]; cH2O  = 645.4139
wCo   = [102,152]; cCo   = 645.5001
wCa2  = [18,89];   cCa2  = 645.5602

sf64_as2 = spc.SpectraFactory("data/6449_aS2",framerows=780,framecols=1438)
sf64_as2.frame_row_cut([0,779])
sf64_as2.frame_col_cut([0])
sf64_as2.contrast_cut(50)
sf64_as2.set_continua("segments")

as2 = sf64_as2.make_spectra()

