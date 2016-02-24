#encoding: utf8
import matplotlib.pyplot as pl
import interactive as intr
import visualize as vis
import danframe as dan
import kontin as con
import lines as lin
import numpy as np

reload(vis);reload(dan);reload(lin);reload(intr)

s6405_t5p = dan.frameseries("data/6405_aS1","top 5%")
s6405_seg = dan.frameseries("data/6405_aS1","segments")

lins = lin.make_lines_from_wins(s6405_t5p,s6405_t5p.pkwindows)

FeI  = lins[1]
SiFe = lins[2]
myst = lins[3]

mes = {}
mes["FeI top 5%"]  = FeI.measure_linecores(s6405_t5p)
mes["FeI segm"  ]  = FeI.measure_linecores(s6405_seg)
mes["Myst top 5%"] = myst.measure_linecores(s6405_t5p)
mes["Myst segm"  ] = myst.measure_linecores(s6405_seg)

mes["SiFe top 5%"]  = SiFe.measure_linecores(s6405_t5p)
mes["SiFe segm"  ]  = SiFe.measure_linecores(s6405_seg)

if False: #SiFe scatter
    pl.plot(mes["SiFe top 5%"][lin.lc.cont,:],mes["SiFe top 5%"][lin.lc.lbot,:],'bo',alpha=0.2);
    pl.title("Bottom, SiFe " + str(SiFe))
    pl.ylabel("Relative intensity of line bottom")
    pl.xlabel("Continuum value at line centre")
    pl.show()

    # Outliers
    # [332] [ 0.93543975] [ 640.71571066] Problem is that they have flat bottoms
    # [795] [ 1.10894355] [ 640.73377294] so that the narrow fit ends up poorly constrained in the
    #                                     spectral direction. Wonder how smoothing would work on that. 

#mes["SiFe top 5% by line"]  = SiFe.byline_measure_linecores(s6405_t5p)
#intr.select_linecore(mes["SiFe top 5% by line"])

vis.show_line_and_corefit(SiFe,s6405_t5p.frames[0],332)


