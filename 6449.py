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

qu1 = sf64_as1.make_spectra()
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1] # Define continua for qu1 series

qu1Myst = spc.splineline(wMyst, cMyst, qu1.meta,"Unkn"   ) # Define lines for qu1 series...
qu1Ca   = spc.splineline(wCa,   cCa,   qu1.meta,"CaI"    )
qu1CoBl = spc.splineline(wCoBl, cCoBl, qu1.meta,"CoI+Bl" )
qu1Unk2 = spc.splineline(wUnk2, cUnk2, qu1.meta,"SiI?VI" )
qu1H2O  = spc.splineline(wH2O,  cH2O,  qu1.meta,"tel H2O")
qu1Co   = spc.splineline(wCo,   cCo,   qu1.meta,"CoI"    )
qu1Ca2  = spc.splineline(wCa2,  cCa2,  qu1.meta,"CaI"    )
qu1H2O2 = spc.splineline(wH2O2, cH2O2, qu1.meta,"tel H2O")

qu1lines = [qu1Myst,qu1Ca,qu1CoBl,qu1Unk2,qu1H2O,qu1Co,qu1Ca2,qu1H2O2]

qu1lims = {}
qu1lims["ewlim"]   = ( 0   , 2  )
qu1lims["vellim"]  = (-3   , 8  )
qu1lims["rellim"]  = ( 0.4 , 1.0)
qu1lims["fw13lim"] = (-0.1 , 0.7)
qu1lims["fwhmlim"] = (-0.1 , 1.0)
qu1lims["fw23lim"] = (-0.1 , 1.2)
qu1lims["as13lim"] = (-0.3 , 0.2)
qu1lims["as12lim"] = (-0.8 , 0.4)
qu1lims["as23lim"] = (-0.8 , 0.4)

wMyst = [1210,1302]; cMyst = 644.9127
wCa   = [1144,1190]; cCa   = 644.9820
wCoBl = [1080,1144]; cCoBl = 645.0179
wUnk2 = [832,874];   cUnk2 = 645.2315
wH2O  = [608,628];   cH2O  = 645.4139
wCo   = [492,530];   cCo   = 645.5001
wCa2  = [412,452];   cCa2  = 645.5602
wH2O2 = [3,30];      cH2O2 = 645.8892

sf64_bs1 = spc.SpectraFactory("data/6449_bS1",framerows=802,framecols=1514)
sf64_bs1.frame_row_cut([0,1,2,799,800,801])
sf64_bs1.frame_col_cut([0])
sf64_bs1.contrast_cut(50)
sf64_bs1.set_continua("segments")

qu2 = sf64_bs1.make_spectra()
qu2con = qu2.meta.cont[0]*qu2.lmbd.mean() + qu2.meta.cont[1] # Define continua for qu2 series

qu2Myst = spc.splineline(wMyst, cMyst, qu2.meta,"Unkn"   )   # Define lines for qu2 series...
qu2Ca   = spc.splineline(wCa,   cCa,   qu2.meta,"CaI"    )
qu2CoBl = spc.splineline(wCoBl, cCoBl, qu2.meta,"CoI+Bl" )
qu2Unk2 = spc.splineline(wUnk2, cUnk2, qu2.meta,"SiI/VI" )
qu2H2O  = spc.splineline(wH2O,  cH2O,  qu2.meta,"tel-H2O")
qu2Co   = spc.splineline(wCo,   cCo,   qu2.meta,"CoI"    )
qu2Ca2  = spc.splineline(wCa2,  cCa2,  qu2.meta,"CaI"    )
qu2H2O2 = spc.splineline(wH2O2, cH2O2, qu2.meta,"tel_H2O")
qu2lines = [qu2Myst,qu2Ca,qu2CoBl,qu2Unk2,qu2H2O,qu2Co,qu2Ca2,qu2H2O2]

qu2lims = {}
qu2lims["ewlim"]   = (-0.8, 5.1 )
qu2lims["vellim"]  = (-6  , 4.8 )
qu2lims["rellim"]  = ( 0.90,1.02)
qu2lims["fw13lim"] = (-0.1, 1.1 )
qu2lims["fwhmlim"] = (-0.08,1.4 )
qu2lims["fw23lim"] = (-0.08,1.43)
qu2lims["as13lim"] = (-0.7, 0.4 )
qu2lims["as12lim"] = (-0.9, 0.5 )
qu2lims["as23lim"] = (-0.9, 0.5 )

###
# Spot section
###
wMyst = [831,938]; cMyst = 644.9127
wCa   = [764,830]; cCa   = 644.9820
wCoBl = [695,763]; cCoBl = 645.0179
wUnk2 = [444,492]; cUnk2 = 645.2315
wH2O  = [229,249]; cH2O  = 645.4139
wCo   = [102,152]; cCo   = 645.5001
wCa2  = [18,89];   cCa2  = 645.5602

sf64_as2 = spc.SpectraFactory("data/6449_aS2",framerows=780,framecols=1438)
sf64_as2.frame_row_cut([0,779])
sf64_as2.frame_col_cut([0])
sf64_as2.contrast_cut(50)
sf64_as2.set_continua("segments")

spt = sf64_as2.make_spectra()
sptcon = spt.meta.cont[0]*spt.lmbd.mean() + spt.meta.cont[1] # Define continua for spt series

sptMyst = spc.splineline(wMyst, cMyst, spt.meta,"Unkn"   ) # Define lines for spt series...
sptCa   = spc.splineline(wCa,   cCa,   spt.meta,"CaI"    )
sptCoBl = spc.splineline(wCoBl, cCoBl, spt.meta,"CoI+Bl" )
sptUnk2 = spc.splineline(wUnk2, cUnk2, spt.meta,"SiI/VI" )
sptH2O  = spc.splineline(wH2O,  cH2O,  spt.meta,"tel-H2O")
sptCo   = spc.splineline(wCo,   cCo,   spt.meta,"CoI"    )
sptCa2  = spc.splineline(wCa2,  cCa2,  spt.meta,"CaI."   )
sptlines = [sptMyst,sptCa,sptCoBl,sptUnk2,sptCo,sptCa2]

sptlims = {}
sptlims["ewlim"]   = ( 0.0 , 2.1 )
sptlims["vellim"]  = (-1.1 , 4.8 )
sptlims["rellim"]  = ( 0.2 , 0.8 )
sptlims["fw13lim"] = ( 0.05, 0.5 )
sptlims["fwhmlim"] = ( 0.05, 0.6 )
sptlims["fw23lim"] = ( 0.10, 0.7 )
sptlims["as13lim"] = (-0.1 , 0.11)
sptlims["as12lim"] = (-0.1 , 0.11)
sptlims["as23lim"] = (-0.1 , 0.11)

um = 0.38834
wl = 0.67615
pn = 0.84928

xspotlims = (644.71474436201663, 645.12763836853765)
yspotlims = (0.67393994279930325, 1.0507659350912264)

