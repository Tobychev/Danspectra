import matplotlib.pyplot as pl
import spectra as spc
vel = 0; bot = 1; cont = 2; err = 3; ew  = 4; mn  = 5; var = 6; ske = 7; kur = 8; wvar = 9; wske = 10; wkur = 11;

wCN   = [289, 342];   cCN   = 640.86681253796701 
wFeI  = [376, 441];   cFeI  = 640.80263785584464
wSiFe = [467, 519];   cSiFe = 640.72889305314993
wMyst = [647, 750];   cMyst = 640.57558287431198
wCNq  = [1010, 1039]; cCNq  = 640.31162454661717

sf = spc.SpectraFactory("data/6405_aS1")
sf.frame_row_cut([0,799])
sf.contrast_cut(50)
sf.set_continua("segments")

sp  = sf.make_spectra()
FeI = spc.statline(wFeI,cFeI,sp.meta)
mesFeI =  FeI.measure(sp)

# Comparing line shape relative to a cut in continua 
# Don't actually need to measure on mystery line to make the plot
if True:
    cut  =  sf.make_spectra_subset(sp,(mesFeI[cont]>1.05).reshape(-1),desc="Continua greater than 1.05")
    acut =  sf.make_spectra_subset(sp,(mesFeI[cont]<=1.05).reshape(-1),desc="Continua less than 1.05")

    pl.plot(cut.meta.lmbd,sp[:,:].mean(axis=0),label="mean")
    pl.plot(cut.meta.lmbd,cut[:,:].mean(axis=0),label="> 1.05,{} spectra".format(cut[:,:].shape[0]))
    pl.plot(cut.meta.lmbd,acut[:,:].mean(axis=0),label="< 1.05, {} spectra".format(acut[:,:].shape[0]))
    pl.plot(cut.meta.lmbd,cut.meta.ref,label="ref")
    pl.xlim(640.538,640.596)
    pl.ylim(0.9337,0.9935)
    pl.legend()
    pl.show()
