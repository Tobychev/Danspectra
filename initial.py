import matplotlib.pyplot as pl
import kontin as con
import numpy as np
import danspec as d
import testwin as tw

reload(d)
reload(tw)
reload(con)

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()



# lines - spectra
# rows  - dispersion/space axis
fil = "6405_aS1_388_cor.fits"
spec = d.danspectra(fil)


fil2 = "6405_aS1_397_cor.fits"
spec2 = d.danspectra(fil2)

try:
    wn
except NameError:
    wn = tw.gen_all_auto_wins(spec2)
    wn = tw.gen_man_win(spec2,wn)

fits = tw.fit_all_windows(spec2,wn)
pt = tw.perturb_all_windows(spec2,wn,5000)

# Continium fit stability, meanspectra
if True:   
    pl.plot(pt["top 5%"][:,0], pt["top 5%"][:,3],label="top 5%")
    pl.plot(pt["manual"][:,0], pt["manual"][:,3],label="manual")
    pl.plot(pt["over 1"][:,0], pt["over 1"][:,3],label="over 1")
    pl.plot(pt["top 20"][:,0], pt["top 20"][:,3],label="top 20")
    pl.title("Continium fit stability on mean spectra, 5k iterations")
    pl.ylabel("std/mean")
    pl.xlabel("Fraction of points included")
    pl.legend()
    pl.show()
