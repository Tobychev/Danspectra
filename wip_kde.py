from wipfile import *

np.seterr(invalid="raise")
s6405_t5p.normalize()
mes = {}
FeImes   = FeI.measure(s6405_t5p)

vel = 0; bot = 1; con = 2; err = 3; ew  = 4; mn  = 5; var = 6; ske = 7; kur = 8
cuts = FeImes[err] < np.percentile(FeImes[err],89)

import pyqt_fit as fitter
import pyqt_fit.nonparam_regression as smooth
#k0 = smooth.NonParamRegression(FeImes[con][cuts],FeImes[ew][cuts],
#                                method=fitter.npr_methods.SpatialAverage())
#k1 = smooth.NonParamRegression(FeImes[con][cuts],FeImes[ew][cuts],
#                                method=fitter.npr_methods.LocalPolynomialKernel(q=1))
#k0.fit()
#k1.fit()

import statsmodels.nonparametric as npar
import statsmodels.nonparametric.kde as kde


kcon = kde.KDEUnivariate(FeImes[con][cuts])
kcon.fit()
kew = kde.KDEUnivariate(FeImes[ew][cuts])
kew.fit()

a, xs = np.histogram(FeImes[con][cuts],63)
a, ys = np.histogram(FeImes[ew][cuts],63)

pl.figure(1)
pl.plot(FeImes[con][cuts],FeImes[ew][cuts],'bo',alpha=0.02);
#pl.plot(xs,kcon.evaluate(xs),label="Continuum density")
#pl.plot(kew.evaluate(ys),ys,label="EW density")
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("Equivalent width")
pl.xlabel("Continuum intensity")
pl.legend(loc='best')
pl.show()

#yopts = k1(FeImes[con][cuts])
#res = FeImes[ew][cuts] - yopts
#fitter.plot_fit.plot_residual_tests(FeImes[con][cuts], yopts, res, 'Spatial Average')
