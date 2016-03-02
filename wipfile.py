import matplotlib.pyplot as pl
import interactive as intr
import visualize as vis
import danframe as dan
import kontin as con
import lines as lin
import numpy as np

reload(vis);reload(dan);reload(lin);reload(intr)

s6405_t5p = dan.frameseries("data/6405_aS1","top 5%")

lins = lin.make_lines_from_wins(s6405_t5p,s6405_t5p.pkwindows)

FeI  = lins[1]
SiFe = lins[2]
myst = lins[3]


frm = s6405_t5p.frames[0]

spec = frm.data[331,:]
lam  = frm.group.lmbd

line = myst
