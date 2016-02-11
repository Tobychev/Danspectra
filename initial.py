import matplotlib.pyplot as pl
import kontin as con
import numpy as np
import danspec as d
import testwin as tw
import visualize as vis

reload(d)
reload(tw)
reload(con)

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()



# rows     - spectra
# columns  - dispersion/space axis
fil2 = "local_data/6405_aS1_397_cor.fits"
spec2 = d.danspectra(fil2)

wn = tw.gen_all_auto_wins(spec2)
wn = tw.gen_man_win(spec2,wn)

#refmean,refstd = tw.compare_win_continua(spec2,wn,centre="median",plot=False,bins=44)
#refmean,refstd = tw.compare_win_continua(spec2,wn,centre="smoothmax",plot=False)
#menmean,menstd = tw.compare_win_continua(spec2,wn,centre="mean",plot=False)

result = tw.test_fit_with_noise(spec2,rows=1e4,plot_excess=True,bins=103,cut=25)

