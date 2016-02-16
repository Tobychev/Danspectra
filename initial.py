import matplotlib.pyplot as pl
import kontin as con
import numpy as np
import danspec as d
import testwin as tw
import visualize as vis
import lines as lin
import interactive as intr
reload(d)
reload(tw)
reload(con)
reload(lin)
reload(vis)
reload(intr)

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()

# rows     - spectra
# columns  - dispersion/space axis
fil2 = "local_data/6405_aS1_397_cor.fits"
spec2 = d.danspectra(fil2)

wn = tw.gen_all_auto_wins(spec2)
wn = tw.gen_man_win(spec2,wn)

lines = lin.make_lines_from_wins(spec2,spec2.pkwindows)

lines = lines + lin.subdivide_line(lines[0])
lines = lines[1:]
