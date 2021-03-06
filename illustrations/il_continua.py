import matplotlib.pyplot as pl
import numpy as np
import spectra as spc
import visualize as vis
import copy

sf64_as1 = spc.SpectraFactory("data/6405_aS1")
sf64_as1.frame_row_cut([0,799])
sf64_as1.frame_col_cut([0])
sf64_as1.contrast_cut(50)
# Hand picked range from meanspectra
as1flat = list(range(147,282))+list(range(1150,1272)) 
as1  = sf64_as1.make_spectra()
cont_as1 = as1[:,as1flat].mean(axis=1)

sf64_bs1 = spc.SpectraFactory("data/6405_bS1")
sf64_bs1.frame_row_cut([0,799])
sf64_bs1.contrast_cut(50)
bs1 = sf64_bs1.make_spectra()
bs1flat = list(range(148,199))+list(range(591,625))+list(range(1194,1240))
cont_bs1 = bs1[:,bs1flat].mean(axis=1)

sf64_as2 = spc.SpectraFactory("data/6405_aS2",framerows=774,framecols=1446)
sf64_as2.frame_row_cut([0,743])
sf64_as2.frame_col_cut([0])
sf64_as2.contrast_cut(85)
as2flat = list(range(197,229))+list(range(404,423))+list(range(972,995))
as2 = sf64_as2.make_spectra()
cont_as2 = as2[:,as2flat].mean(axis=1)

if False:
    pl.plot(as1.lmbd,as1[:,:].mean(axis=0),label="Quiet Sun, series a")
    pl.plot(bs1.lmbd,bs1[:,:].mean(axis=0),label="Quiet Sun, series b")
    pl.plot(as2.lmbd,as2[:,:].mean(axis=0),label="Sun spot")
    pl.xlabel("Wavelength [nm]")
    pl.ylabel("Relative intensity")
    pl.legend(loc="lower left")
    pl.show()

if False:
    ax = pl.subplot(111)
    ax = vis.kde(cont_as1,ax)
    ax = vis.kde(cont_bs1,ax)
    ax = vis.kde(cont_as2,ax)
    ax.lines[0].set_label("Quiet Sun, series a")
    ax.lines[1].set_label("Quiet Sun, series b")
    ax.lines[2].set_label("Sun spot")
    ax.legend(loc="upper left")
    ax.set_xlabel("Continuum value")
    ax.set_ylabel("Number of spectra")
    ax.set_title("Distribution of continuum value")
    ax.figure.show()

def selected(cont,spec):
    mnspec = spec[:,:].mean(axis=0)
    pl.plot(spec.lmbd,mnspec)
    pl.plot(cont.lmbd,mnspec[cont.idx],'ro')
    pl.show()

def normalize(block,lmbd,continua):
    return block/continua(lmbd,block)

# Defining
as1_const   = copy.deepcopy(as1)
as1_manual  = copy.deepcopy(as1); con_as1_const   = spc.continua(np.array(as1flat),as1.lmbd,"manual")
as1_top20   = copy.deepcopy(as1); con_as1_top20   = spc.continua(as1[:,:].mean(axis=0),as1.lmbd,"top 20")
as1_segment = copy.deepcopy(as1); con_as1_segment = spc.continua(as1[:,:].mean(axis=0),as1.lmbd,"segments")

# Normalize
as1_const.modify(lambda x: x/cont_as1.reshape(-1,1))
as1_manual.modify(lambda x:  normalize(x,as1.lmbd,con_as1_const))
as1_top20.modify(lambda x:   normalize(x,as1.lmbd,con_as1_top20))
as1_segment.modify(lambda x: normalize(x,as1.lmbd,con_as1_segment))

# Evaluate, special windows to not evaluate on training set
cas1_const = as1_const[:,list(range(584,622))].mean(axis=1)
cas1_manua = as1_manual[:,list(range(584,622))].mean(axis=1)
cas1_top20 = as1_top20[:,list(range(584,622))].mean(axis=1)
cas1_segme = as1_segment[:,list(range(584,622))].mean(axis=1)
cas1_orig  = as1[:,list(range(584,622))].mean(axis=1)

if True:
    fig,ax = pl.subplots(1)
    ax = vis.kde(cas1_const,ax)
    ax = vis.kde(cas1_manua,ax)
    ax = vis.kde(cas1_top20,ax)
    ax = vis.kde(cas1_segme,ax)
    ax = vis.kde(cas1_orig,ax)
    ax.lines[0].set_label("Constant")
    ax.lines[1].set_label("Manual")
    ax.lines[2].set_label("Top 20")
    ax.lines[3].set_label("Segment")
    ax.lines[4].set_label("Original")
    ax.legend(loc="upper right")
    ax.set_xlabel("Continuum value")
    ax.set_ylabel("Number of spectra")
    ax.set_title("Distribution of continuum value")
    fig.show()

# Defining
as2_const   = copy.deepcopy(as2)
as2_manual  = copy.deepcopy(as2); con_as2_const   = spc.continua(np.array(as2flat),as2.lmbd,"manual")
as2_top20   = copy.deepcopy(as2); con_as2_top20   = spc.continua(as2[:,:].mean(axis=0),as2.lmbd,"top 20")
as2_segment = copy.deepcopy(as2); con_as2_segment = spc.continua(as2[:,:].mean(axis=0),as2.lmbd,"segments")

# Normalize
as2_const.modify(lambda x: x/cont_as2.reshape(-1,1))
as2_manual.modify(lambda x:  normalize(x,as2.lmbd,con_as2_const))
as2_top20.modify(lambda x:   normalize(x,as2.lmbd,con_as2_top20))
as2_segment.modify(lambda x: normalize(x,as2.lmbd,con_as2_segment))

# Evaluate, special windows to not evaluate on training set
cas2_const = as2_const[:,list(range(738,758))].mean(axis=1)
cas2_manua = as2_manual[:,list(range(738,758))].mean(axis=1)
cas2_top20 = as2_top20[:,list(range(738,758))].mean(axis=1)
cas2_segme = as2_segment[:,list(range(738,758))].mean(axis=1)
cas2_orig  = as2[:,list(range(738,758))].mean(axis=1)

if True:
    fig,ax = pl.subplots(1)
    ax = vis.kde(cas2_const,ax)
    ax = vis.kde(cas2_manua,ax)
    ax = vis.kde(cas2_top20,ax)
    ax = vis.kde(cas2_segme,ax)
    ax = vis.kde(cas2_orig,ax)
    ax.lines[0].set_label("Constant")
    ax.lines[1].set_label("Manual")
    ax.lines[2].set_label("Top 20")
    ax.lines[3].set_label("Segment")
    ax.lines[4].set_label("Original")
    ax.legend(loc="upper left")
    ax.set_xlabel("Continuum value")
    ax.set_ylabel("Number of spectra")
    ax.set_title("Distribution of continuum value")
    fig.show()
