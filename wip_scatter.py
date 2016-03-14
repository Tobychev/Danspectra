from wipfile import *

np.seterr(invalid="raise")
s6405_t5p.normalize()
mes = {}
FeImes   = FeI.measure(s6405_t5p)
#SiFe_t5_lcen, SiFe_t5_lbot, SiFe_t5_cont, SiFe_t5_err = SiFe.measure_linecores(s6405_t5p)
#myst_t5_lcen, myst_t5_lbot, myst_t5_cont, myst_t5_err = myst.measure_linecores(s6405_t5p)
#CN_t5_lcen  , CN_t5_lbot  , CN_t5_cont  , CN_t5_err   = CN.measure_linecores(s6405_t5p)
#CNq_t5_lcen , CNq_t5_lbot , CNq_t5_cont , CNq_t5_err  = CNq.measure_linecores(s6405_t5p)

vel = 0; bot = 1; con = 2; err = 3; ew  = 4; mn  = 5; var = 6; ske = 7; kur = 8
cuts = FeImes[err] < np.percentile(FeImes[err],89)

pl.plot(FeImes[con][cuts],FeImes[bot][cuts],'bo',alpha=0.2);
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("Line min intesity")
pl.xlabel("Continuum intensity")
pl.show()


pl.plot(FeImes[con][cuts],FeImes[var][cuts],'bo',alpha=0.2);
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("Variance")
pl.xlabel("Continuum intensity")
pl.show()

pl.plot(FeImes[con][cuts],FeImes[kur][cuts],'bo',alpha=0.2);
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("Kurtosis")
pl.xlabel("Continuum intensity")
pl.show()

pl.plot(FeImes[con][cuts],FeImes[ske][cuts],'bo',alpha=0.2);
pl.title("Line bottom , FeI " + str(FeI))
pl.ylabel("Skewness")
pl.xlabel("Continuum intensity")
pl.show()




EW = {}
EW["FeI top 5%"]  = FeI.measure_EW(s6405_t5p)
EW["Myst top 5%"] = myst.measure_EW(s6405_t5p)
EW["SiFe top 5%"] = SiFe.measure_EW(s6405_t5p)
EW["CN top 5%"]   = CN.measure_EW(s6405_t5p)
EW["CN? top 5%"]   = CNq.measure_EW(s6405_t5p)
