import spectra as spc

Myst = {"El":    -1,"gf":     0,"lam":644.9127,"dep":0.102,"name":"Myst"}
Ca   = {"El":2.5213,"gf":-0.502,"lam":644.9820,"dep":0.560,"name":"Ca I"}
Blnd = {"El":1.7104,"gf":-1.527,"lam":645.0179,"dep":0.374,"name":"Co I Si I"}
SiI  = {"El":5.6193,"gf":-1.433,"lam":645.2315,"dep":0.211,"name":"Si I"}
H2O  = {"El":    -1,"gf":     0,"lam":645.4139,"dep":0.026,"name":"tel H_2O"}
Co   = {"El":3.6321,"gf":-0.250,"lam":645.5001,"dep":0.141,"name":"Co I"}
H2O2 = {"El":    -1,"gf":     0,"lam":645.8892,"dep":0.046,"name":"tel H_2O"}

wMyst = [1211,1327];
wCa   = [1150,1198];
wBlnd = [1088,1151];
wSiI  = [834,889];  
wH2O  = [613,637];  
wCo   = [449,534];  
wH2O2 = [8,34];     

sf_qu1 = spc.SpectraFactory("data/6449_aS1",framerows=802,framecols=1514)
sf_qu1.frame_row_cut([0]+list(range(668,674))+[801])
sf_qu1.frame_col_cut([0,1513])
sf_qu1.contrast_cut(50)
sf_qu1.set_continua("segments")

qu1 = sf_qu1.make_spectra()
qu1con = qu1.meta.cont[0]*qu1.lmbd.mean() + qu1.meta.cont[1] # Define continua for qu1 series

qu1Myst = spc.splineline(wMyst, Myst, qu1.meta) # Define lines for qu1 series...
qu1Ca   = spc.splineline(wCa,   Ca,   qu1.meta)
qu1Blnd = spc.splineline(wBlnd, Blnd, qu1.meta)
qu1SiI  = spc.splineline(wSiI,  SiI,  qu1.meta)
qu1H2O  = spc.splineline(wH2O,  H2O,  qu1.meta)
qu1Co   = spc.splineline(wCo,   Co,   qu1.meta)
qu1H2O2 = spc.splineline(wH2O2, H2O2, qu1.meta)
qu1lines = [qu1Myst,qu1Ca,qu1Blnd,qu1SiI,qu1H2O,qu1Co,qu1H2O2]

wMyst = [1210,1302]
wCa   = [1144,1190]
wBlnd = [1080,1144]
wSiI  = [832,874]; 
wH2O  = [608,628]; 
wCo   = [492,530]; 
wH2O2 = [3,30];    

sf_qu2 = spc.SpectraFactory("data/6449_bS1",framerows=802,framecols=1514)
sf_qu2.frame_row_cut([0,1,2]+list(range(671,681))+[801])
sf_qu2.frame_col_cut([0,1513])
sf_qu2.contrast_cut(50)
sf_qu2.set_continua("segments")

qu2 = sf_qu2.make_spectra()
qu2con = qu2.meta.cont[0]*qu2.lmbd.mean() + qu2.meta.cont[1] # Define continua for qu2 series

qu2Myst = spc.splineline(wMyst, Myst, qu2.meta)   # Define lines for qu2 series...
qu2Ca   = spc.splineline(wCa,   Ca,   qu2.meta)
qu2Blnd = spc.splineline(wBlnd, Blnd, qu2.meta)
qu2SiI  = spc.splineline(wSiI,  SiI,  qu2.meta)
qu2H2O  = spc.splineline(wH2O,  H2O,  qu2.meta)
qu2Co   = spc.splineline(wCo,   Co,   qu2.meta)
qu2H2O2 = spc.splineline(wH2O2, H2O2, qu2.meta)
qu2lines = [qu2Myst,qu2Ca,qu2Blnd,qu2SiI,qu2H2O,qu2Co,qu2H2O2]

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

sptMyst = spc.splineline(wMyst, Myst, spt.meta) # Define lines for spt series...
sptCa   = spc.splineline(wCa,   Ca,   spt.meta)
sptBlnd = spc.splineline(wBlnd, Blnd, spt.meta)
sptSiI  = spc.splineline(wSiI,  SiI,  spt.meta)
sptH2O  = spc.splineline(wH2O,  H2O,  spt.meta)
sptCo   = spc.splineline(wCo,   Co,   spt.meta)
sptlines = [sptMyst,sptCa,sptBlnd,sptSiI,sptH2O,sptCo]

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

