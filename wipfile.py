import matplotlib.pyplot as pl
import interactive as intr
import visualize as vis
import danframe as dan
import kontin as con
import lines as lin
import numpy as np
import imp

imp.reload(vis);imp.reload(dan);imp.reload(lin);imp.reload(intr)

s6405_t5p = dan.frameseries("data/6405_aS1","top 5%")

lins = lin.make_lines_from_wins(s6405_t5p,s6405_t5p.pkwindows)

CN   = lins[0]
FeI  = lins[1]
SiFe = lins[2]
myst = lins[3]
CNq  = lins[4]
