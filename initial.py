import matplotlib.pyplot as pl
from . import kontin as con
import numpy as np
from . import danframe as dan
from . import testwin as tw
from . import visualize as vis
from . import lines as lin
from . import interactive as intr
import imp
imp.reload(dan)
imp.reload(tw)
imp.reload(con)
imp.reload(lin)
imp.reload(vis)
imp.reload(intr)

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()

# rows     - spectra
# columns  - dispersion/space axis
#wn = tw.gen_all_auto_wins(spec2)
#wn = tw.gen_man_win(spec2,wn)


fil1 = "local_data/6405_aS1_397_cor.fits"
spec1 = dan.danframe_sac(fil1)
lines = lin.make_lines_from_wins(spec1,spec1.pkwindows)

l1 = lines[2]

#lin.fit_linecores(l1,spec2.lmbd,spec2.data)

#lines = lines + lin.subdivide_line(lines[0])
#lines = lines[1:]

#con6405 = con.continua(spec2)

fil = "local_data/6405_aS1"
spec6405 = dan.frameseries(fil,"top 5%")

l6407 = l1.measure_linecores(spec6405)
