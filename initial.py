import matplotlib.pyplot as pl
import kontin as con
import numpy as np
import danframe as dan
import testwin as tw
import visualize as vis
import lines as lin
import interactive as intr
reload(dan)
reload(tw)
reload(con)
reload(lin)
reload(vis)
reload(intr)

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
spec2 = dan.frameseries(fil,"top 5%")

