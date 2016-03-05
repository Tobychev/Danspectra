#encoding: utf8
from wipfile import * 

import statsmodels.api as sm
spec = frm.data[331,:]
kde  = sm.nonparametric.KDEUnivariate(1-spec[line.idx])
kdeM = sm.nonparametric.KDEMultivariate(1-spec[line.idx],bw="cv_ml",var_type="c")
kde.fit(kernel="epa",bw="scott",fft=False)
pl.plot(kde.support,kde.density);pl.show()

