import matplotlib.pyplot as pl
import kontin as con
import numpy as np
import danspec as d

reload(d)
reload(con)

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()



# lines - spectra
# rows  - dispersion/space axis
fil = "6405_aS1_388_cor.fits"
spec = d.danspectra(fil)


fil2 = "6405_aS1_397_cor.fits"
spec2 = d.danspectra(fil2)

