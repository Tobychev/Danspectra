#encoding: utf8
import matplotlib.pyplot as pl
from . import interactive as intr
from . import visualize as vis
from . import danframe as dan
from . import kontin as con
from . import lines as lin
import numpy as np

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


#Konstiga mätningar
# Första bilden
# Myst segm   : 10,406,693, (695)
# Myst top 5% : 10,406,693, (695)
# FeI  segm   : 313,314,315,316,317,799

# 799 - helt klart skräp
# 406,693,695 - ett extremvärde vid randen


