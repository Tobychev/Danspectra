import spectra as spc

wMyst = [1211,1327]; cMyst = 644.9127
wCa   = [1150,1198]; cCa   = 644.9820
wBlnd = [1088,1151]; cBlnd = 645.0179
wSiI  = [834,889];   cSiI  = 645.2315
wH2O  = [613,637];   cH2O  = 645.4139
wCo   = [449,534];   cCo   = 645.5001
wH2O2 = [8,34];      cH2O2 = 645.8892

sf_qu1 = spc.SpectraFactory("data/6449_aS1",framerows=802,framecols=1514)
sf_qu1.frame_row_cut([0]+list(range(668,674))+[801])
sf_qu1.frame_col_cut([0,1513])
sf_qu1.contrast_cut(50)
sf_qu1.set_continua("segments")

qu1 = sf_qu1.make_spectra()
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1] # Define continua for qu1 series

qu1Myst = spc.splineline(wMyst, cMyst, qu1.meta,"Myst"   ) # Define lines for qu1 series...
qu1Ca   = spc.splineline(wCa,   cCa,   qu1.meta,"CaI"    )
qu1Blnd = spc.splineline(wBlnd, cBlnd, qu1.meta,"CoI+Bl" )
qu1SiI  = spc.splineline(wSiI,  cSiI,  qu1.meta,"SiI" )
qu1H2O  = spc.splineline(wH2O,  cH2O,  qu1.meta,"tel H2O")
qu1Co   = spc.splineline(wCo,   cCo,   qu1.meta,"CoI"    )
qu1H2O2 = spc.splineline(wH2O2, cH2O2, qu1.meta,"tel H2O")

qu1lines = [qu1Myst,qu1Ca,qu1CoBl,qu1Unk2,qu1H2O,qu1Co,qu1H2O2]

wMyst = [1210,1302]; cMyst = 644.9127
wCa   = [1144,1190]; cCa   = 644.9820
wBlnd = [1080,1144]; cBlnd = 645.0179
wSiI  = [832,874];   cSiI  = 645.2315
wH2O  = [608,628];   cH2O  = 645.4139
wCo   = [492,530];   cCo   = 645.5001
wH2O2 = [3,30];      cH2O2 = 645.8892

sf_qu2 = spc.SpectraFactory("data/6449_bS1",framerows=802,framecols=1514)
sf_qu2.frame_row_cut([0,1,2]+list(range(671,681))+[801])
sf_qu2.frame_col_cut([0,1513])
sf_qu2.contrast_cut(50)
sf_qu2.set_continua("segments")

qu2 = sf_qu2.make_spectra()
qu2con = qu2.meta.cont[0]*qu2.lmbd.mean() + qu2.meta.cont[1] # Define continua for qu2 series

qu2Myst = spc.splineline(wMyst, cMyst, qu2.meta,"Myst"   )   # Define lines for qu2 series...
qu2Ca   = spc.splineline(wCa,   cCa,   qu2.meta,"CaI"    )
qu2Blnd = spc.splineline(wBlnd, cBlnd, qu2.meta,"Blend" )
qu2SiI  = spc.splineline(wSiI,  cSiI,  qu2.meta,"SiI" )
qu2H2O  = spc.splineline(wH2O,  cH2O,  qu2.meta,"tel-H2O")
qu2Co   = spc.splineline(wCo,   cCo,   qu2.meta,"CoI"    )
qu2H2O2 = spc.splineline(wH2O2, cH2O2, qu2.meta,"tel_H2O")
qu2lines = [qu2Myst,qu2Ca,qu2CoBl,qu2Unk2,qu2H2O,qu2Co,qu2H2O2]


###
# Spot section
###
wMyst = [831,938]; cMyst = 644.9127
wCa   = [764,830]; cCa   = 644.9820
wBlnd = [695,763]; cBlnd = 645.0179
wSiI  = [444,492]; cSiI  = 645.2315
wH2O  = [229,249]; cH2O  = 645.4139
wCo   = [102,152]; cCo   = 645.5001

sf_spt = spc.SpectraFactory("data/6449_aS2",framerows=780,framecols=1438)
sf_spt.frame_row_cut([0]+list(range(661,672))+[779])
sf_spt.frame_col_cut([0,1437])
sf_spt.contrast_cut(50)
sf_spt.set_continua("segments")

spt = sf_spt.make_spectra()
sptcon = spt.meta.cont[0]*spt.lmbd.mean() + spt.meta.cont[1] # Define continua for spt series

sptMyst = spc.splineline(wMyst, cMyst, spt.meta,"Myst"   ) # Define lines for spt series...
sptCa   = spc.splineline(wCa,   cCa,   spt.meta,"CaI"    )
sptBlnd = spc.splineline(wBlnd, cBlnd, spt.meta,"Blend" )
sptSiI  = spc.splineline(wSiI,  cSiI,  spt.meta,"SiI" )
sptH2O  = spc.splineline(wH2O,  cH2O,  spt.meta,"tel-H2O")
sptCo   = spc.splineline(wCo,   cCo,   spt.meta,"CoI"    )
sptCa2  = spc.splineline(wCa2,  cCa2,  spt.meta,"CaI."   )
sptlines = [sptMyst,sptCa,sptCoBl,sptUnk2,sptCo]

um = 0.38834
wl = 0.67615
pn = 0.84928

xspotlims = (644.71474436201663, 645.12763836853765)
yspotlims = (0.67393994279930325, 1.0507659350912264)

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

