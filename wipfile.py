import matplotlib.pyplot as pl
from . import interactive as intr
from . import visualize as vis
from . import danframe as dan
from . import kontin as con
from . import lines as lin
import numpy as np
import imp

imp.reload(vis);imp.reload(dan);imp.reload(lin);imp.reload(intr)

s6405_t5p = dan.frameseries("data/6405_aS1","top 5%")

lins = lin.make_lines_from_wins(s6405_t5p,s6405_t5p.pkwindows)

FeI  = lins[1]
SiFe = lins[2]
myst = lins[3]


frm = s6405_t5p.frames[0]

frm.data = frm.data/frm.cont.norm()

spec = frm.data[331,:]
lam  = frm.group.lmbd

line = myst
