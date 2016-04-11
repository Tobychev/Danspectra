import matplotlib.pyplot as pl
import spectra as spc
from imp import reload
reload(spc)

sf = spc.SpectraFactory("data/6405_aS1")
sf.frame_row_cut([0,799])
sf.contrast_cut(50)
sf.set_continua("segments")

bl = sf.make_spectra()

wCN   = [289, 342];   cCN   = 640.86681253796701 
wFeI  = [376, 441];   cFeI  = 640.80263785584464
wSiFe = [467, 519];   cSiFe = 640.72889305314993
wMyst = [647, 750];   cMyst = 640.57558287431198
wCNq  = [1010, 1039]; cCNq  = 640.31162454661717

FeI = spc.statline(wFeI,cFeI,bl.meta)
sFeI = spc.splineline(wFeI,cFeI,bl.meta)

Myst = spc.statline(wMyst,cMyst,bl.meta)
SiFe = spc.statline(wSiFe,cSiFe,bl.meta)

vel = 0; bot = 1; cont = 2; err = 3; ew  = 4; mn  = 5; var = 6; ske = 7; kur = 8; wvar = 9; wske = 10; wkur = 11;
mesFeI =  FeI.measure(bl)
smesFeI = sFeI.measure(bl,mesFeI[cont])


pl.plot(mesFeI[cont],mesFeI[ew],'o',alpha=0.1)
pl.show()

pl.plot(mesFeI[cont],smesFeI[:,3],'o',alpha=0.1)
pl.show()

