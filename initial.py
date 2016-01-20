import matplotlib.pyplot as pl
import numpy as np
import danspec as d

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()



# lines - spectra
# rows  - dispersion/space axis
fil = "6405_aS1_388_cor.fits"
spec = d.danspectra(fil)



