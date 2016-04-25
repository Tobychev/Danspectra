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

lims = {}
lims["ewlim"]   = ( 0 , 2  )
lims["vellim"]  = (-5 , 8  )
lims["rellim"]  = ( 0.1 , 1.3  )
lims["fw13lim"] = ( - 0.1, 1.1  )
lims["fwhmlim"] = (  -0.2 , 1.1 )
lims["fw23lim"] = (-0.2 , 1.1    )
lims["as13lim"] = (-1.4 , 0.4  )
lims["as12lim"] = (-1.5 , 0.5  )
lims["as23lim"] = (-1.3 , 0.4  )

as1lims=lims
###OBS! Co is crazy, not included in the interval

sf64_as1 = spc.SpectraFactory("data/6449_aS1",framerows=802,framecols=1514)
sf64_as1.frame_row_cut([0,801])
sf64_as1.frame_col_cut([0])
sf64_as1.contrast_cut(50)
sf64_as1.set_continua("segments")

as1 = sf64_as1.make_spectra()
as1con = as1.meta.cont[0]*as1.lmbd.mean() + as1.meta.cont[1] # Define continua for as1 series

as1Myst = spc.splineline(wMyst, cMyst, as1.meta,"Myst ") # Define lines for as1 series...
as1Ca = spc.splineline(wCa, cCa, as1.meta,"Ca ")
as1CoBl  = spc.splineline(wCoBl , cCoBl , as1.meta,"CoBl " )
as1Unk2   = spc.splineline(wUnk2  , cUnk2  , as1.meta,"Unk2 "  )
as1H2O = spc.splineline(wH2O, cH2O, as1.meta,"H2O ")
as1Co  = spc.splineline(wCo, cCo , as1.meta,"Co " )
as1Ca2 = spc.splineline(wCa2, cCa2, as1.meta,"Ca2 ")
as1H2O2  = spc.splineline(wH2O2, cH2O2 , as1.meta,"H2O2 " )

as1lines = [as1Myst,as1Ca,as1CoBl,as1Unk2,as1H2O,as1Co,as1Ca2,as1H2O2]

bs1lims = lims
sf64_bs1 = spc.SpectraFactory("data/6449_bS1",framerows=802,framecols=1514)
sf64_bs1.frame_row_cut([0,779])
sf64_bs1.frame_col_cut([0])
sf64_bs1.contrast_cut(50)
sf64_bs1.set_continua("segments")

bs1 = sf64_bs1.make_spectra()
bs1con = bs1.meta.cont[0]*bs1.lmbd.mean() + bs1.meta.cont[1] # Define continua for bs1 series

bs1Myst = spc.splineline(wMyst, cMyst, bs1.meta,"Myst ") # Define lines for bs1 series...
bs1Ca = spc.splineline(wCa, cCa, bs1.meta,"Ca ")
bs1CoBl  = spc.splineline(wCoBl , cCoBl , bs1.meta,"CoBl " )
bs1Unk2   = spc.splineline(wUnk2  , cUnk2  , bs1.meta,"Unk2 "  )
bs1H2O = spc.splineline(wH2O, cH2O, bs1.meta,"H2O ")
bs1Co  = spc.splineline(wCo, cCo , bs1.meta,"Co " )
bs1Ca2 = spc.splineline(wCa2, cCa2, bs1.meta,"Ca2 ")
bs1H2O2  = spc.splineline(wH2O2, cH2O2 , bs1.meta,"H2O2 " )

### OBS!!! Unk2 and H2O cause crash
bs1lines = [bs1Myst,bs1Ca,bs1CoBl,bs1Co,bs1Ca2,bs1H2O2]
#bs1lines = [bs1Myst,bs1Ca,bs1CoBl,bs1Unk2,bs1H2O,bs1Co,bs1Ca2,bs1H2O2]

wMyst = [831,938]; cMyst = 644.9127
wCa   = [764,830]; cCa   = 644.9820
wCoBl = [695,763]; cCoBl = 645.0179
wUnk2 = [444,492]; cUnk2 = 645.2315
wH2O  = [227,254]; cH2O  = 645.4139
wCo   = [102,152]; cCo   = 645.5001
wCa2  = [18,89];   cCa2  = 645.5602

lims["ewlim"]   = ( -0.7 , 3.6  )
lims["vellim"]  = (-8 , 10  )
lims["rellim"]  = ( -1 , 6  )
lims["fw13lim"] = ( - 0.1, 1.1  )
lims["fwhmlim"] = (  -0.2 , 1.1 )
lims["fw23lim"] = (-0.2 , 1.1    )
lims["as13lim"] = (-1.4 , 0.4  )
lims["as12lim"] = (-1.5 , 0.5  )
lims["as23lim"] = (-1.3 , 0.4  )
as2lims = lims

sf64_as2 = spc.SpectraFactory("data/6449_aS2",framerows=780,framecols=1438)
sf64_as2.frame_row_cut([0,779])
sf64_as2.frame_col_cut([0])
sf64_as2.contrast_cut(50)
sf64_as2.set_continua("segments")

as2 = sf64_as2.make_spectra()
as2con = as2.meta.cont[0]*as2.lmbd.mean() + as2.meta.cont[1] # Define continua for as2 series

as2Myst = spc.splineline(wMyst, cMyst, as2.meta,"Myst ") # Define lines for as2 series...
as2Ca = spc.splineline(wCa, cCa, as2.meta,"Ca ")
as2CoBl  = spc.splineline(wCoBl , cCoBl , as2.meta,"CoBl " )
as2Unk2   = spc.splineline(wUnk2  , cUnk2  , as2.meta,"Unk2 "  )
as2H2O = spc.splineline(wH2O, cH2O, as2.meta,"H2O ")
as2Co  = spc.splineline(wCo, cCo , as2.meta,"Co " )
as2Ca2 = spc.splineline(wCa2, cCa2, as2.meta,"Ca2 ")

#OBS!!! H2O causes crash
#as2lines = [as2Myst,as2Ca,as2CoBl,as2Unk2,as2H2O,as2Co,as2Ca2]
as2lines = [as2Myst,as2Ca,as2CoBl,as2Unk2,as2Co,as2Ca2]

