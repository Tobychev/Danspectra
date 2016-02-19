import matplotlib.pyplot as pl
import danframe as dan
import kontin as con
import lines as lin
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
